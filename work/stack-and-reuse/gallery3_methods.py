"""Gallery figures 27-30: methods, chemistry and file formats (descriptive only).

27-isobaric-plex   isobaric labelling plex stated in metadata text, 2013-2025 (100% stacked among datasets
                   stating a plex), mean channels per labelled dataset (text and table headers), and the
                   isobaric share of all datasets.
28-ptm-landscape   which modifications appear in deposited MaxQuant site tables, and how often each
                   metadata source reflects them.
29-species         organisms named in PRIDE species metadata, 2013-2025.
30-file-formats    file types present in the top-level file listing of each dataset, 2013-2025, and the
                   open-spectra decline split by PRIDE submission type.

Rules: `standards` is the corrected nongenerated_expanded flag (the legacy inclusive flag is never used).
Contributor groups are hashed `pi_group` only. Uncertainty: 95% percentile interval from a group
bootstrap resampling `pi_group` (datasets without a group are their own cluster), seed 20261001,
2,000 replicates (ratio of sums). Absence of a signal means "not detected". Claims are descriptive.
"""
import re, sys
import numpy as np, pandas as pd
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3
from common import O, SEED, cluster_boot, R

D = G / 'data'; D.mkdir(exist_ok=True)
H = A / 'H5-H6-H8-strategies-software-docs'
HALO = [pe.withStroke(linewidth=2.0, foreground='white')]


def headmm(fig, title, sub, h_mm):
    """House headline with offsets fixed in millimetres (as in gallery2_stack.py)."""
    if record_head(title, sub):
        return
    fig.text(.012, 1 - 1.6 / h_mm, title, fontsize=10.5, fontweight='bold', va='top', color=INK)
    fig.text(.012, 1 - 6.2 / h_mm, sub, fontsize=7, va='top', color=MUTED, linespacing=1.25)
    r = fig.canvas.get_renderer()
    for tx in fig.texts[-2:]:
        wmm = tx.get_window_extent(r).width / fig.dpi * 25.4
        if wmm > 181:
            print(f'WARNING headline text {wmm:.0f} mm wide: {tx.get_text()[:60]}')


def gci(values, groups):
    lo, hi = cluster_boot(np.asarray(values, float), np.asarray(groups), reps=2000, seed=SEED)
    return lo, hi


def strat():
    d = pd.read_csv(H / 'Dataset-strategy-software-docs.csv.gz')
    d['cluster'] = d.pi_group.fillna('acc:' + d.accession)
    return d


# ------------------------------------------------------------------------------------------------
# 27. Isobaric plex
# ------------------------------------------------------------------------------------------------
PLEX_GROUPS = [  # (label, text_plex values, colour)
    ('iTRAQ 4', ['iTRAQ 4'], '#f2c46d'),
    ('iTRAQ 8', ['iTRAQ 8'], '#d99a2b'),
    ('TMT 2/6', ['TMT 2', 'TMT 6'], '#bcd3e6'),
    ('TMT 10/11', ['TMT 10', 'TMT 11'], '#4f86b8'),
    ('TMTpro 16', ['TMT 16'], '#b23a48'),
    ('TMTpro 18/35', ['TMT 18', 'TMT 35'], '#6e1423'),
    ('TMTpro, plex not stated', ['TMTpro, plex not stated'], '#e3a4ab'),
]
PLEX_N = {'iTRAQ 4': 4, 'iTRAQ 8': 8, 'TMT 2': 2, 'TMT 6': 6, 'TMT 10': 10, 'TMT 11': 11, 'TMT 16': 16, 'TMT 18': 18, 'TMT 35': 35}


def fig_isobaric(d):
    iso = d[d.isobaric_meta & d.release_year.between(2013, 2025)].copy()
    # cross-check against the frozen per-year table
    fz = pd.read_csv(H / 'H5-isobaric-plex-text-by-year.csv').set_index('release_year')
    assert (iso.groupby('release_year').size().reindex(YEARS) == fz.isobaric_datasets.reindex(YEARS)).all()
    st = iso[iso.text_plex.notna() & (iso.text_plex != 'TMT and iTRAQ both stated')]
    cnt = pd.DataFrame({lab: st[st.text_plex.isin(v)].groupby('release_year').size().reindex(YEARS).fillna(0)
                        for lab, v, _ in PLEX_GROUPS})
    n_state = cnt.sum(axis=1)
    share = 100 * cnt.div(n_state, axis=0)
    # period numbers for the title (ledger H5.5 periods)
    per = lambda y0, y1: 100 * cnt.loc[y0:y1].sum() / n_state.loc[y0:y1].sum()
    p0, p1 = per(2013, 2016), per(2023, 2025)
    itraq0, itraq1 = p0['iTRAQ 4'] + p0['iTRAQ 8'], p1['iTRAQ 4'] + p1['iTRAQ 8']
    pro1 = p1['TMTpro 16'] + p1['TMTpro 18/35'] + p1['TMTpro, plex not stated']
    # mean channels per dataset (text-stated plex; TMTpro without plex excluded)
    iso['plex_n'] = iso.text_plex.map(PLEX_N)
    rows = []
    for y in YEARS:
        s = iso[(iso.release_year == y) & iso.plex_n.notna()]
        lo, hi = gci(s.plex_n, s.cluster)
        hd = d[(d.release_year == y) & (d.h_iso_channels > 0)]
        hlo, hhi = gci(hd.h_iso_channels, hd.cluster) if len(hd) >= 20 else (np.nan, np.nan)
        a = d[d.release_year == y]
        alo, ahi = gci(a.isobaric_meta, a.cluster)
        rows.append(dict(release_year=y, isobaric_datasets=int(a.isobaric_meta.sum()), all_datasets=len(a),
                         isobaric_pct=100 * a.isobaric_meta.mean(), isobaric_lo=100 * alo, isobaric_hi=100 * ahi,
                         plex_stating=int(n_state[y]), mean_text_plex=s.plex_n.mean(), text_lo=lo, text_hi=hi,
                         n_text_plex_numeric=len(s), header_datasets=len(hd),
                         mean_header_channels=hd.h_iso_channels.mean() if len(hd) else np.nan, header_lo=hlo, header_hi=hhi))
    m = pd.DataFrame(rows).set_index('release_year')
    out = share.add_prefix('pct_').join(cnt.add_prefix('n_')).join(m)
    out.round(3).to_csv(D / '27-isobaric-plex.csv')

    y16 = iso[iso.release_year.between(2013, 2016) & iso.plex_n.notna()].plex_n.mean()
    y25 = m.mean_text_plex[2025]
    HM = 100; fig = plt.figure(figsize=(W, HM * MM))
    headmm(fig, f'iTRAQ gave way to TMTpro, and the mean isobaric experiment grew from {m.mean_text_plex[2014]:.0f} to {y25:.0f} channels',
           f'PRIDE datasets with isobaric labelling in metadata, 2013–2025 (n = {len(iso):,}). a, plex stated in text, among the {int(n_state.sum()):,} that state one\n'
           f'(highest TMT plex stated is used); iTRAQ {itraq0:.0f}% → {itraq1:.0f}% and TMTpro 0% → {pro1:.0f}% from 2013–16 to 2023–25. b, mean channels per dataset:\n'
           f'plex stated in text vs reporter channels counted in table headers (independent; channels used ≤ kit size; years with ≥20 datasets).\n'
           f'c, isobaric datasets among all datasets. Shading: 95% CI, group bootstrap over contributor groups.', HM)
    ax = fig.add_axes([.068, .09, .452, .655])
    x = np.array(YEARS, float)
    base = np.zeros(len(x)); right = {}
    inside = {'iTRAQ 4': 2014.4, 'iTRAQ 8': 2015.6, 'TMT 2/6': 2017.2, 'TMT 10/11': 2020.6}
    for lab, _, col in PLEX_GROUPS:
        v = share[lab].to_numpy()
        ax.fill_between(x, base, base + v, color=col, lw=.4, edgecolor='white')
        dark = col in ('#4f86b8', '#b23a48', '#6e1423', '#d99a2b')
        if lab in inside:
            xl = inside[lab]; yc = np.interp(xl, x, base + v / 2)
            ax.text(xl, yc, lab, ha='center', va='center', fontsize=6.5, fontweight='bold', color='white' if dark else INK)
        right[lab] = (base[-1] + v[-1] / 2, v[-1], col)
        base = base + v
    pos = repel({k: v[0] for k, v in right.items()}, 7.5)
    lo = min(pos.values())
    if lo < 2:  # keep the bottom label inside the axis
        pos = {k: v + (2 - lo) for k, v in pos.items()}
    for lab, (yc, v, col) in right.items():
        txt = {'TMTpro, plex not stated': 'TMTpro, plex\nnot stated', 'TMTpro 18/35': 'TMTpro 18+'}.get(lab, lab)
        ax.annotate(f'{txt} {v:.0f}%', xy=(2025, yc), xytext=(2025.45, pos[lab]), fontsize=6.2, va='center', ha='left',
                    color=col if col not in ('#bcd3e6', '#f2c46d', '#e3a4ab') else {'#bcd3e6': '#6f93b3', '#f2c46d': '#b07d1c', '#e3a4ab': '#c0606b'}[col],
                    fontweight='bold', linespacing=.95, annotation_clip=False,
                    arrowprops=dict(arrowstyle='-', color=HAIR, lw=.5, shrinkA=0, shrinkB=0))
    ax.set_xlim(2013, 2025); ax.set_ylim(0, 100); yfmt(ax)
    ax.set_xticks([2013, 2015, 2017, 2019, 2021, 2023, 2025])
    ax.set_ylabel('Share of plex-stating isobaric datasets')
    for y in YEARS:
        ax.text(y, 101.5, f'{int(n_state[y])}', ha='center', va='bottom', fontsize=6, color=MUTED)
    ax.text(2012.75, 101.5, 'n', ha='right', va='bottom', fontsize=6, color=MUTED, style='italic')
    panel(ax, 'a', x=-.085, y=1.05)

    # b: channels per dataset
    bx = fig.add_axes([.667, .09, .135, .655])
    ok = m.n_text_plex_numeric >= 20
    yy = m.index[ok]
    bx.fill_between(yy, m.text_lo[ok], m.text_hi[ok], color=ACC, alpha=.18, lw=0)
    bx.plot(yy, m.mean_text_plex[ok], color=ACC, lw=1.6)
    hk = m.header_datasets >= 20
    bx.fill_between(m.index[hk], m.header_lo[hk], m.header_hi[hk], color=ACC2, alpha=.18, lw=0)
    bx.plot(m.index[hk], m.mean_header_channels[hk], color=ACC2, lw=1.4, ls=(0, (3, 1.2)))
    bx.text(2025, m.mean_text_plex[2025] + 1.2, f'{m.mean_text_plex[2025]:.0f}', color=INK, fontsize=7.5, fontweight='bold', ha='center')
    y0 = int(yy.min())
    bx.text(y0 + .45, m.mean_text_plex[y0] - 2.0, f'{m.mean_text_plex[y0]:.0f}', color=ACC, fontsize=7.5, fontweight='bold', ha='center')
    bx.text(2014.0, 10.6, 'stated in text', color=ACC, fontsize=6.3, ha='left')
    bx.text(2020.5, 5.0, 'counted in\nheaders', color=ACC2, fontsize=6.3, ha='left', linespacing=.95)
    bx.set_xlim(2013.4, 2025.8); bx.set_ylim(0, 18); bx.set_xticks([2014, 2019, 2025])
    bx.set_title('Mean channels', loc='left', fontsize=7, color=INK)
    panel(bx, 'b', x=-.2, y=1.05)

    # c: isobaric share of all datasets
    cx = fig.add_axes([.857, .09, .12, .655])
    cx.fill_between(YEARS, m.isobaric_lo, m.isobaric_hi, color=MUTED, alpha=.22, lw=0)
    cx.plot(YEARS, m.isobaric_pct, color=INK, lw=1.6)
    pk = int(m.isobaric_pct.loc[2016:].idxmax())
    for y, yr, dy in [(pk, pk, 3.0), (2024.2, 2025, -4.4)]:
        cx.text(y, m.isobaric_pct[yr] + dy, f'{m.isobaric_pct[yr]:.0f}%', ha='center', fontsize=7.5, fontweight='bold', color=INK)
    cx.text(2013.2, m.isobaric_pct[2013] - 4.2, f'{m.isobaric_pct[2013]:.0f}%', ha='left', fontsize=7.5, fontweight='bold', color=INK)
    cx.set_ylim(0, 25); yfmt(cx); cx.set_xlim(2012.6, 2025.6); cx.set_xticks([2013, 2019, 2025])
    cx.set_title('Isobaric share', loc='left', fontsize=7, color=INK)
    panel(cx, 'c', x=-.26, y=1.05)
    save(fig, '27-isobaric-plex')
    return dict(itraq0=itraq0, itraq1=itraq1, pro1=pro1, p0=p0.round(1).to_dict(), p1=p1.round(1).to_dict(),
                n_state0=int(n_state.loc[2013:2016].sum()), n_state1=int(n_state.loc[2023:2025].sum()), mean16=y16, mean25=y25)


# ------------------------------------------------------------------------------------------------
# 28. PTM landscape
# ------------------------------------------------------------------------------------------------
MODS = [('phosphorylation', 'Phosphorylation'), ('ubiquitin remnant (GlyGly)', 'Ubiquitin remnant (GlyGly)'),
        ('lysine acetylation', 'Lysine acetylation'), ('other lysine acylations', 'Other lysine acylations'),
        ('methylation', 'Methylation'), ('glycosylation', 'Glycosylation'), ('ADP-ribosylation', 'ADP-ribosylation'),
        ('SUMO remnant', 'SUMO remnant'), ('citrullination', 'Citrullination')]
SRC = [('struct', 'PRIDE PTM field', ACC2, 'o'), ('design', 'Title, description or keywords', ACC3, 's'),
       ('enrich', 'Enrichment named in protocol', ACC, 'D')]


def ptm_field():
    """PRIDE identified-PTM field (names, lower case), same source as the frozen arsenal script."""
    sys.path.insert(0, str(R / 'scripts')); from query import connect
    return connect().execute("""select s.accession, lower(coalesce(array_to_string(list_transform(b.identifiedPTMStrings, x -> x.name), ' | '), '')) ptm
        from study s left join pride_bronze_projects b on b.accession = s.accession""").df()


def fig_ptm(d):
    m = d[d.mq_sites_visible].copy()
    N = len(m)
    m = m.merge(ptm_field(), on='accession', how='left', validate='one_to_one')
    m['ptm_generic_acetyl'] = m.ptm.fillna('').str.contains('acetyl')
    lists = m.site_mods.fillna('').str.split('; ')
    rows = []
    for k, lab in MODS:
        has = lists.apply(lambda l: k in l)
        lo, hi = gci(has, m.cluster)
        rows.append(dict(panel='a', item=lab, source='', n=N, k=int(has.sum()), pct=100 * has.mean(), lo=100 * lo, hi=100 * hi))
    a = pd.DataFrame(rows)
    # b: of datasets with a site table for the modification, share where each metadata source reflects it
    rb = []
    for k, lab in MODS[:3]:
        st = m[m[f'site_{k}']]
        for src, slab, _, _ in SRC:
            v = st[f'{src}_{k}']
            lo, hi = gci(v, st.cluster)
            rb.append(dict(panel='b', item=lab, source=slab, n=len(st), k=int(v.sum()), pct=100 * v.mean(), lo=100 * lo, hi=100 * hi))
    # generic acetyl term in the PTM field (also used for N-terminal acetylation): reported for honesty
    st = m[m['site_lysine acetylation']]
    gen = st.ptm_generic_acetyl
    lo, hi = gci(gen, st.cluster)
    rb.append(dict(panel='b', item='Lysine acetylation', source='PRIDE PTM field, generic acetyl term', n=len(st), k=int(gen.sum()),
                   pct=100 * gen.mean(), lo=100 * lo, hi=100 * hi))
    gen_other = 100 * m.loc[~m['site_lysine acetylation'], 'ptm_generic_acetyl'].mean()
    b = pd.DataFrame(rb)
    # c: share with a phospho site table, by what the metadata says about phospho
    t, e, f = m.design_phosphorylation, m.enrich_phosphorylation, m.struct_phosphorylation
    cats = [('Text + enrichment', t & e), ('PTM field only', ~t & ~e & f),
            ('Text, no enrichment', t & ~e), ('Enrichment, no text', ~t & e),
            ('None of the three', ~t & ~e & ~f)]
    rc = []
    for lab, mask in cats:
        s_ = m[mask]; v = s_.site_phosphorylation
        lo, hi = gci(v, s_.cluster)
        rc.append(dict(panel='c', item=lab.replace('\n', ' '), source='', n=len(s_), k=int(v.sum()), pct=100 * v.mean(), lo=100 * lo, hi=100 * hi))
    c = pd.DataFrame(rc)
    assert sum(mk.sum() for _, mk in cats) == N
    pd.concat([a, b, c]).round(2).to_csv(D / '28-ptm-landscape.csv', index=False)

    bb = b.set_index(['item', 'source']).pct
    HM = 88; fig = plt.figure(figsize=(W, HM * MM))
    headmm(fig, f"PRIDE's PTM field flags {bb[('Phosphorylation', 'PRIDE PTM field')]:.0f}% of deposited phospho site tables "
                f"but only {bb[('Ubiquitin remnant (GlyGly)', 'PRIDE PTM field')]:.0f}% of GlyGly ones",
           f'PRIDE datasets 2013–2025 whose deposited MaxQuant output shows site tables (n = {N:,}; the default Oxidation (M)Sites table makes them visible).\n'
           f'A site table means the modification was searched as variable, not that it was enriched. a, modifications with a site table. b, of datasets with\n'
           f'that site table, share whose metadata names the modification. c, share with a phospho site table, by which metadata sources mention phospho.\n'
           f'Whiskers: 95% CI, group bootstrap over contributor groups. Absence means "not detected". MaxQuant only.', HM)
    top = .695; hgt = .55; bot = top - hgt
    # a
    ax = fig.add_axes([.165, bot, .16, hgt])
    yy = np.arange(len(a))[::-1]
    cols = [ACC2 if i == 0 else '#9fb7cc' for i in range(len(a))]
    ax.barh(yy, a.pct, color=cols, height=.68)
    ax.errorbar(a.pct, yy, xerr=[a.pct - a.lo, a.hi - a.pct], fmt='none', ecolor=INK, elinewidth=.6, capsize=0)
    for y_, r in zip(yy, a.itertuples()):
        ax.text(r.hi + .6, y_, f'{r.pct:.1f}%' if r.pct < 10 else f'{r.pct:.0f}%', va='center', fontsize=6.3, color=INK)
        ax.text(-.6, y_, r.item, va='center', ha='right', fontsize=6.5, color=INK)
    ax.set_yticks([]); ax.set_xlim(0, 22); ax.set_ylim(-.6, len(a) - .4)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%')); ax.set_xticks([0, 10, 20])
    ax.set_title(f'Share of {N:,} datasets', loc='left', fontsize=7, color=INK)
    panel(ax, 'a', x=-.9, y=1.03)
    # b
    bx = fig.add_axes([.455, bot, .2, hgt])
    labs = [lab for _, lab in MODS[:3]]
    for i, lab in enumerate(labs):
        y0 = 2 - i
        bx.axhline(y0, color=HAIR, lw=.5, zorder=0)
        for j, (src, slab, col, mk) in enumerate(SRC):
            r = b[(b.item == lab) & (b.source == slab)].iloc[0]
            yj = y0 + (1 - j) * .22
            bx.plot([r.lo, r.hi], [yj, yj], color=col, lw=1.1, solid_capstyle='butt')
            bx.plot(r.pct, yj, mk, color=col, ms=4.2, mec='white', mew=.5, zorder=3)
            right = r.hi < 82
            if i == 0:  # direct source labels on the first row
                bx.text(r.lo - 3, yj, {'struct': 'PTM field', 'design': 'Text', 'enrich': 'Enrichment'}[src],
                        va='center', ha='right', fontsize=6.1, color=col)
                right = True
            bx.text(r.hi + 3 if right else r.lo - 3, yj, f'{r.pct:.0f}%', va='center', ha='left' if right else 'right',
                    fontsize=6.2, color=col, fontweight='bold')
            if lab == 'Lysine acetylation' and src == 'struct':
                g = b[b.source == 'PRIDE PTM field, generic acetyl term'].iloc[0]
                bx.plot([g.lo, g.hi], [yj, yj], color=col, lw=1.1, solid_capstyle='butt')
                bx.plot(g.pct, yj, 'o', mfc='white', mec=col, ms=4.2, mew=1, zorder=3)
                bx.text(g.lo - 3, yj, f'generic "acetyl"\n{g.pct:.0f}%', va='center', ha='right', fontsize=6, color=col, linespacing=.95)
        n_ = int(b[b.item == lab].n.iloc[0])
        short = {'Ubiquitin remnant (GlyGly)': 'GlyGly', 'Lysine acetylation': 'Lysine\nacetylation'}.get(lab, lab)
        bx.text(-7, y0, short + f'\nn = {n_:,}',
                ha='right', va='center', fontsize=6.5, linespacing=1.05)
    bx.set_yticks([]); bx.set_ylim(-.55, 2.55); bx.set_xlim(-3, 103)
    bx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%')); bx.set_xticks([0, 50, 100])
    bx.set_title('Metadata names it', loc='left', fontsize=7, color=INK)
    panel(bx, 'b', x=-.36, y=1.03)
    bx.text(-.02, -.11, f'Open circle: the field\'s generic "acetylated residue" term,\nalso present in {gen_other:.0f}% of datasets without a lysine-acetyl table',
            transform=bx.transAxes, fontsize=6, color=MUTED, ha='left', va='top', linespacing=1.05)
    # c
    cx = fig.add_axes([.825, bot, .16, hgt])
    yy = np.arange(len(c))[::-1]
    cc = [ACC2, '#6f93b3', '#9fb7cc', '#9fb7cc', '#c9c6be']
    cx.barh(yy, c.pct, color=cc, height=.66)
    cx.errorbar(c.pct, yy, xerr=[c.pct - c.lo, c.hi - c.pct], fmt='none', ecolor=INK, elinewidth=.6)
    for y_, r in zip(yy, c.itertuples()):
        inside = r.pct > 70
        cx.text(4 if inside else r.hi + 3, y_, f'{r.pct:.0f}%' if r.pct >= 10 else f'{r.pct:.1f}%', va='center',
                ha='left', fontsize=6.5, color='white' if inside else INK, fontweight='bold')
        cx.text(-3, y_, f'{r.item}\nn = {r.n:,}', va='center', ha='right', fontsize=6.3, linespacing=1.05)
    cx.set_yticks([]); cx.set_xlim(0, 100); cx.set_ylim(-.6, len(c) - .4)
    cx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%')); cx.set_xticks([0, 50, 100])
    cx.set_title('Phospho site table present', loc='left', fontsize=7, color=INK)
    cx.text(1.0, -.11, 'Text = title, description or keywords;\nenrichment = named in sample protocol', transform=cx.transAxes,
            fontsize=6, color=MUTED, ha='right', va='top', linespacing=1.05)
    panel(cx, 'c', x=-.62, y=1.03)
    save(fig, '28-ptm-landscape')
    return pd.concat([a, b, c])


# ------------------------------------------------------------------------------------------------
# 29. Species
# ------------------------------------------------------------------------------------------------
def norm_species(n):
    """Collapse strains/subspecies to the binomial; keep non-binomial labels (e.g. 'Bacteria') as given."""
    n = re.sub(r'\s*\(.*', '', str(n)).strip().lower()
    if n.startswith('severe acute respiratory syndrome coronavirus 2'):
        return 'SARS-CoV-2'
    t = n.split()
    if len(t) >= 2 and t[1] != 'sp.':
        return ' '.join(t[:2]).capitalize()
    return n.capitalize()


def fig_species():
    sys.path.insert(0, str(R / 'scripts')); from query import connect
    sp = connect().execute('select project_accession accession, species_name from pride_species').df()
    sp['sp'] = sp.species_name.map(norm_species)
    f = F[F.release_year.between(2013, 2025)][['accession', 'release_year', 'pi_group']].copy()
    f['cluster'] = f.pi_group.fillna('acc:' + f.accession)
    sp = sp[sp.accession.isin(f.accession)].drop_duplicates(['accession', 'sp'])
    n_none = int((~f.accession.isin(sp.accession)).sum())
    f['hs'] = f.accession.isin(sp[sp.sp == 'Homo sapiens'].accession)
    f['mm'] = f.accession.isin(sp[sp.sp == 'Mus musculus'].accession)
    f['other'] = f.accession.isin(sp[~sp.sp.isin(['Homo sapiens', 'Mus musculus'])].accession)
    f['neither'] = ~f.hs & ~f.mm
    n_multi = int((sp.groupby('accession').size() > 1).sum())
    # a: yearly shares with CI
    rows = []
    for y in YEARS:
        s_ = f[f.release_year == y]
        for col, lab in [('hs', 'Human'), ('mm', 'Mouse'), ('neither', 'Neither human nor mouse')]:
            lo, hi = gci(s_[col], s_.cluster)
            rows.append(dict(panel='a', release_year=y, series=lab, n=len(s_), k=int(s_[col].sum()), pct=100 * s_[col].mean(),
                             lo=100 * lo, hi=100 * hi))
    a = pd.DataFrame(rows)
    # b: period shares for the next-most-common species
    P0, P1 = (2013, 2016), (2023, 2025)
    f0 = f[f.release_year.between(*P0)]; f1 = f[f.release_year.between(*P1)]
    top = sp[~sp.sp.isin(['Homo sapiens', 'Mus musculus'])].sp.value_counts().head(14).index
    rb = []
    for name in top:
        acc = set(sp[sp.sp == name].accession)
        r = dict(panel='b', series=name)
        for tag, ff in [('p0', f0), ('p1', f1)]:
            v = ff.accession.isin(acc)
            lo, hi = gci(v, ff.cluster)
            r.update({f'{tag}_n': len(ff), f'{tag}_k': int(v.sum()), f'{tag}_pct': 100 * v.mean(), f'{tag}_lo': 100 * lo, f'{tag}_hi': 100 * hi})
        r['ratio'] = r['p1_pct'] / r['p0_pct']
        rb.append(r)
    b = pd.DataFrame(rb).sort_values('ratio')
    # within-group check: groups depositing in both periods
    g0, g1 = f0[f0.pi_group.notna()], f1[f1.pi_group.notna()]
    both = set(g0.pi_group) & set(g1.pi_group)
    w0, w1 = g0[g0.pi_group.isin(both)], g1[g1.pi_group.isin(both)]
    wc = [(100 * w.neither.mean(), *(100 * np.array(gci(w.neither, w.cluster)))) for w in (w0, w1)]
    wrow = pd.DataFrame([dict(panel='within-group', series=f'Neither human nor mouse, {len(both)} groups in both periods',
                              p0_n=len(w0), p0_pct=wc[0][0], p0_lo=wc[0][1], p0_hi=wc[0][2],
                              p1_n=len(w1), p1_pct=wc[1][0], p1_lo=wc[1][1], p1_hi=wc[1][2])])
    pd.concat([a, b, wrow]).round(3).to_csv(D / '29-species.csv', index=False)

    A_ = a.pivot(index='release_year', columns='series', values='pct')
    nei0 = 100 * f0.neither.mean(); nei1 = 100 * f1.neither.mean()
    HM = 100; fig = plt.figure(figsize=(W, HM * MM))
    headmm(fig, f'Proteomics is converging on human and mouse: datasets naming neither fell from {nei0:.0f}% to {nei1:.0f}%',
           f'PRIDE datasets 2013–2025 (n = {len(f):,}; {n_none} list no species). Species from PRIDE species metadata, strains collapsed to the binomial;\n'
           f'a dataset naming several species ({n_multi:,}) counts once for each, so shares do not sum to 100%. Title and b: 2013–16 (n = {len(f0):,}) vs 2023–25\n'
           f'(n = {len(f1):,}). Within the {len(both)} contributor groups depositing in both periods, the neither-share also fell, {wc[0][0]:.0f}% ({wc[0][1]:.0f}–{wc[0][2]:.0f}) → '
           f'{wc[1][0]:.0f}% ({wc[1][1]:.0f}–{wc[1][2]:.0f}).\nShading and bars: 95% CI, group bootstrap. Dataset counts rose for every species shown; only their shares fell.', HM)
    top_, hgt = .745, .645; bot = top_ - hgt
    ax = fig.add_axes([.06, bot, .36, hgt])
    cols = {'Human': ACC2, 'Mouse': ACC, 'Neither human nor mouse': MUTED}
    for lab, col in cols.items():
        s_ = a[a.series == lab].set_index('release_year')
        ax.fill_between(YEARS, s_.lo, s_.hi, color=col, alpha=.18, lw=0)
        ax.plot(YEARS, s_.pct, color=col, lw=1.7)
        ax.text(2013 - .2, s_.pct[2013], f'{s_.pct[2013]:.0f}%', ha='right', va='center', fontsize=7, color=col, fontweight='bold')
        ax.text(2025 + .2, s_.pct[2025], f'{s_.pct[2025]:.0f}%', ha='left', va='center', fontsize=7, color=col, fontweight='bold')
    ax.text(2022, A_.loc[2022, 'Human'] + 4.0, 'Human', color=ACC2, fontsize=7, fontweight='bold', ha='center')
    ax.text(2019, A_.loc[2019, 'Mouse'] - 4.2, 'Mouse', color=ACC, fontsize=7, fontweight='bold', ha='center')
    ax.text(2018.3, 31.5, 'Neither human\nnor mouse', color=MUTED, fontsize=7, fontweight='bold', ha='center', linespacing=.95)
    ax.set_xlim(2011.9, 2026.1); ax.set_ylim(0, 62); yfmt(ax); ax.set_xticks([2013, 2016, 2019, 2022, 2025])
    ax.set_ylabel('Share of datasets naming the species')
    panel(ax, 'a', x=-.1, y=1.02)
    # b: dumbbell, log scale
    bx = fig.add_axes([.6, bot, .3, hgt])
    yy = np.arange(len(b))
    for y_, r in zip(yy, b.itertuples()):
        up = r.ratio > 1.15; down = r.ratio < 1 / 1.15
        col = ACC2 if up else (ACC3 if down else '#9a9790')
        bx.plot([r.p0_lo, r.p0_hi], [y_ + .17, y_ + .17], color='#b9b6ae', lw=.8, solid_capstyle='butt', zorder=1)
        bx.plot([r.p1_lo, r.p1_hi], [y_ - .17, y_ - .17], color=col, lw=.8, solid_capstyle='butt', zorder=1)
        bx.annotate('', xy=(r.p1_pct, y_), xytext=(r.p0_pct, y_),
                    arrowprops=dict(arrowstyle='-|>,head_length=.35,head_width=.18', color=col, lw=1.0, shrinkA=0, shrinkB=2), zorder=2)
        bx.plot(r.p0_pct, y_, 'o', ms=3.2, mfc='white', mec=col, mew=.9, zorder=3)
        bx.plot(r.p1_pct, y_, 'o', ms=3.6, color=col, zorder=3)
        bx.text(.105, y_, 'Bacteria (unspecified)' if r.series == 'Bacteria' else r.series, ha='right', va='center', fontsize=6.4, color=INK, style='italic' if ' ' in r.series else 'normal')
        bx.text(8.6, y_, f'{r.p0_k} → {r.p1_k}', ha='left', va='center', fontsize=6, color=MUTED)
    bx.set_xscale('log'); bx.set_xlim(.11, 8); bx.set_ylim(-.7, len(b) - .3)
    bx.set_xticks([.2, .5, 1, 2, 5]); bx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}%'))
    bx.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ''))
    bx.set_yticks([])
    bx.set_title('Next 14 species, 2013–16 (open) → 2023–25 (filled)', loc='left', fontsize=7, color=INK)
    bx.text(8.6, len(b) - .45, 'datasets', ha='left', va='bottom', fontsize=6, color=MUTED)
    bx.set_xlabel('Share of datasets naming the species (log scale)')
    panel(bx, 'b', x=-.5, y=1.02)
    save(fig, '29-species')
    return a, b, wrow


# ------------------------------------------------------------------------------------------------
# 30. File formats
# ------------------------------------------------------------------------------------------------
FORMATS = [  # (group, label, column builder)
    ('Vendor raw', 'Any vendor raw', lambda f: f.fmt_raw_ext | f.fmt_bruker_d | f.fmt_sciex_wiff),
    ('Vendor raw', 'Thermo .raw', lambda f: f.fmt_raw_ext),
    ('Vendor raw', 'Bruker .d', lambda f: f.fmt_bruker_d),
    ('Vendor raw', 'SCIEX .wiff', lambda f: f.fmt_sciex_wiff),
    ('Open spectra', 'Any open spectra', lambda f: f.fmt_mzml | f.fmt_mzxml | f.fmt_peaklist),
    ('Open spectra', 'mzML', lambda f: f.fmt_mzml),
    ('Open spectra', 'mzXML', lambda f: f.fmt_mzxml),
    ('Open spectra', 'Peak lists (MGF…)', lambda f: f.fmt_peaklist),
    ('Results', 'Standard results', lambda f: f.standards),
    ('Results', 'Tool-native outputs', lambda f: f[[c for c in f.columns if c.startswith('tool_')]].any(axis=1)),
    ('Context', 'FASTA', lambda f: f.fmt_fasta),
    ('Context', 'SDRF', lambda f: f.fmt_sdrf),
    ('Context', 'README', lambda f: f.fmt_readme),
]


def fig_formats():
    from matplotlib.colors import PowerNorm
    f = F[F.release_year.between(2013, 2025)].copy()
    assert (f.standards == f.nongenerated_expanded).all()
    f['cluster'] = f.pi_group.fillna('acc:' + f.accession)
    for _, lab, fn in FORMATS:
        f[lab] = fn(f).fillna(False).astype(bool)
    labs = [lab for _, lab, _ in FORMATS]
    M = 100 * f.groupby('release_year')[labs].mean().reindex(YEARS)
    n_y = f.groupby('release_year').size().reindex(YEARS)
    P0, P1 = (2013, 2016), (2023, 2025)
    f0, f1 = f[f.release_year.between(*P0)], f[f.release_year.between(*P1)]
    per = pd.DataFrame({'p0': 100 * f0[labs].mean(), 'p1': 100 * f1[labs].mean()})
    per['change_pp'] = per.p1 - per.p0
    # b/c: open spectra by PRIDE submission type
    f['complete'] = f.submission_type == 'COMPLETE'
    rows = []
    for y in YEARS:
        for comp in (True, False):
            s_ = f[(f.release_year == y) & (f.complete == comp)]
            for lab in ['Any open spectra', 'mzML', 'Peak lists (MGF…)']:
                lo, hi = gci(s_[lab], s_.cluster)
                rows.append(dict(release_year=y, submission='COMPLETE' if comp else 'other (PARTIAL etc.)', series=lab, n=len(s_),
                                 k=int(s_[lab].sum()), pct=100 * s_[lab].mean(), lo=100 * lo, hi=100 * hi))
    st = pd.DataFrame(rows)
    comp_share = 100 * f.groupby('release_year').complete.mean().reindex(YEARS)
    out = M.add_prefix('pct_').assign(datasets=n_y, complete_pct=comp_share)
    out.round(2).to_csv(D / '30-file-formats.csv')
    st.round(2).to_csv(D / '30-file-formats-by-submission-type.csv', index=False)

    o0, o1 = per.loc['Any open spectra', 'p0'], per.loc['Any open spectra', 'p1']
    # shift-share: 2023-25 within-type rates with the 2013-16 submission-type mix
    r1 = f1.assign(c=f1.submission_type == 'COMPLETE').groupby('c')['Any open spectra'].mean()
    mix0 = (f0.submission_type == 'COMPLETE').mean(); mix1 = (f1.submission_type == 'COMPLETE').mean()
    o1_mix0 = 100 * (mix0 * r1[True] + (1 - mix0) * r1[False])
    per.loc['Any open spectra, 2023-25 rates with 2013-16 COMPLETE share', ['p0', 'p1']] = [o0, o1_mix0]
    z0, z1 = per.loc['mzML', 'p0'], per.loc['mzML', 'p1']
    k0, k1 = per.loc['Peak lists (MGF…)', 'p0'], per.loc['Peak lists (MGF…)', 'p1']
    HM = 138; fig = plt.figure(figsize=(W, HM * MM))
    r0 = f0.assign(c=f0.submission_type == 'COMPLETE').groupby('c')['Any open spectra'].mean()
    o1_rate0 = 100 * (mix1 * r0[True] + (1 - mix1) * r0[False])
    per.loc['Any open spectra, 2013-16 rates with 2023-25 COMPLETE share', ['p0', 'p1']] = [o0, o1_rate0]
    per.round(2).to_csv(D / '30-file-formats-periods.csv')
    mix_lo, mix_hi = sorted([o0 - o1_rate0, o1_mix0 - o1])
    headmm(fig, f'Open spectra fell from {o0:.0f}% to {o1:.0f}% of datasets as COMPLETE submissions halved; mzML still rose',
           f'PRIDE datasets 2013–2025 (n = {len(f):,}) with at least one top-level file of each type (file names; compression stripped; /generated/ paths\n'
           f'and zip contents excluded). Periods: 2013–16 (n = {len(f0):,}) vs 2023–25 (n = {len(f1):,}). Standard results = corrected mzIdentML/mzTab definition.\n'
           f'b–c, open spectra by PRIDE submission type (COMPLETE submissions must include peak lists). COMPLETE fell from {100 * mix0:.0f}% to {100 * mix1:.0f}% of datasets,\n'
           f'accounting for {mix_lo:.0f}–{mix_hi:.0f} of the {o0 - o1:.0f}-point fall (shift-share, both orders). Shading: 95% CI, group bootstrap. Presence ≠ validity or reuse.', HM)
    # a: heatmap
    ax = fig.add_axes([.205, .4, .66, .43])
    nrm = PowerNorm(.5, 0, 100); cmap = plt.get_cmap('Blues')
    groups = [g for g, _, _ in FORMATS]
    ypos = []; y = 0
    for i, g in enumerate(groups):
        if i > 0 and g != groups[i - 1]:
            y += .45
        ypos.append(y); y += 1
    for (g, lab, _), yy in zip(FORMATS, ypos):
        for j, yr in enumerate(YEARS):
            v = M.loc[yr, lab]
            ax.add_patch(Rectangle((j, yy), .94, .9, facecolor=cmap(nrm(v)), edgecolor='none'))
            ax.text(j + .47, yy + .45, f'{v:.0f}' if v >= 9.5 else f'{v:.1f}', ha='center', va='center', fontsize=6,
                    color='white' if nrm(v) > .55 else INK, fontweight='bold' if lab.startswith('Any') else 'normal')
        bold = lab.startswith('Any') or lab.startswith('Standard')
        ax.text(-.2, yy + .45, lab, ha='right', va='center', fontsize=6.5, color=INK, fontweight='bold' if bold else 'normal')
        ch = per.loc[lab, 'change_pp']
        col = ACC2 if ch > 1 else (ACC if ch < -1 else MUTED)
        ax.text(13.3, yy + .45, f'{per.loc[lab, "p0"]:.0f} → {per.loc[lab, "p1"]:.0f}', ha='left', va='center', fontsize=6.2, color=INK)
        ax.text(15.75, yy + .45, f'{ch:+.0f}', ha='right', va='center', fontsize=6.5, color=col, fontweight='bold')
    for g in dict.fromkeys(groups):
        ys = [yy for (gg, _, _), yy in zip(FORMATS, ypos) if gg == g]
        ax.plot([-4.15, -4.15], [min(ys) + .1, max(ys) + .8], color=MUTED, lw=.7, clip_on=False)
        ax.text(-4.3, (min(ys) + max(ys) + .9) / 2, g, rotation=90, ha='right', va='center', fontsize=6.3, color=MUTED, clip_on=False)
    for j, yr in enumerate(YEARS):
        ax.text(j + .47, -.2, f"{yr}" if yr in (2013, 2025) else f"’{str(yr)[2:]}", ha='center', va='bottom', fontsize=6.3, color=INK)
    ax.text(13.3, -.2, '13–16 → 23–25', ha='left', va='bottom', fontsize=6.2, color=INK)
    ax.text(15.75, -.2, 'pp', ha='right', va='bottom', fontsize=6.2, color=INK)
    ax.set_xlim(-.1, 15.8); ax.set_ylim(ypos[-1] + 1, -.3)
    ax.axis('off')
    ax.text(0, ypos[-1] + 1.35, '% of datasets with the file type, by release year; right: period shares and change (percentage points)',
            fontsize=6.3, color=MUTED, ha='left', va='top')
    panel(ax, 'a', x=-.255, y=1.0)
    # b, c
    for sub, title, letter, x0 in [('COMPLETE', 'COMPLETE submissions', 'b', .07), ('other (PARTIAL etc.)', 'All other submissions', 'c', .565)]:
        bx = fig.add_axes([x0, .065, .36, .245])
        ends = {}
        for lab, col in [('Any open spectra', INK), ('Peak lists (MGF…)', ACC3), ('mzML', ACC2)]:
            s_ = st[(st.submission == sub) & (st.series == lab)].set_index('release_year')
            bx.fill_between(YEARS, s_.lo, s_.hi, color=col, alpha=.15, lw=0)
            bx.plot(YEARS, s_.pct, color=col, lw=1.5 if lab != 'Any open spectra' else 1.2)
            ends[lab] = (s_.pct[2025], col)
        pos = repel({k: v[0] for k, v in ends.items()}, 8)
        lo_ = min(pos.values())
        if lo_ < 3:
            pos = {k: v + 3 - lo_ for k, v in pos.items()}
        for lab, (v, col) in ends.items():
            bx.text(2025.35, pos[lab], f'{v:.0f}%', ha='left', va='center', fontsize=6.5, color=col, fontweight='bold')
        ncomp = int((f.complete == (sub == 'COMPLETE')).sum())
        bx.set_xlim(2012.6, 2026.9); bx.set_ylim(0, 100); yfmt(bx); bx.set_xticks([2013, 2016, 2019, 2022, 2025])
        bx.set_title(f'{title} (n = {ncomp:,})', loc='left', fontsize=7, color=INK)
        panel(bx, letter, x=-.1, y=1.04)
        if sub == 'COMPLETE':
            bx.text(2013.2, 97, 'any open spectra', color=INK, fontsize=6.3, va='top')
            bx.text(2020.3, 58, 'peak lists', color=ACC3, fontsize=6.3, va='center')
            bx.text(2019.4, 6, 'mzML', color=ACC2, fontsize=6.3, va='bottom')
    save(fig, '30-file-formats')
    return per, comp_share


if __name__ == '__main__':
    which = sys.argv[1:] or ['27', '28', '29', '30']
    d = strat()
    if '27' in which:
        print(fig_isobaric(d))
    if '28' in which:
        print(fig_ptm(d).round(1).to_string())
    if '29' in which:
        for t in fig_species():
            print(t.round(2).to_string())
    if '30' in which:
        for t in fig_formats():
            print(t.round(1).to_string())
