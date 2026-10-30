"""Gallery figures 34-37: the producers and their outlets (descriptive).

34-contributor-concentration  concentration of datasets and bytes over hashed contributor groups (pi_group = stable
                              hash of one PI e-mail/ORCID); top-share curves (log x), datasets-per-group distribution,
                              and top-1% / top-10% shares and Gini in 2013-17 vs 2021-25 (equal 5-year windows).
35-collaboration              international co-authorship of the papers linked to datasets, PAPER grain (one row per
                              OpenAlex work, 2013-2025): share with authors in >=2 countries by paper year (raw and
                              country-mix standardised), and by the depositing group's country against its output.
36-research-fields            OpenAlex primary field / subfield of the linked paper, DEPOSIT grain, by release year;
                              share of datasets with a linked paper recorded as the denominator.
37-journals                   journals (Europe PMC journal title of the OpenAlex work, joined by PMID) of linked papers:
                              top 15 by linked datasets with yearly share sparklines (paper grain), and concentration.

Uncertainty: 95% percentile intervals resampling contributor groups (pi_group; deposits with no group are their own
cluster; a paper is assigned the most frequent group among its linked deposits), 1,000 reps, seed 20261001.
Claims are descriptive only. Never plots names of people or institutions; countries and journal titles only.
"""
import sys, json
import numpy as np, pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter, NullLocator
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3, panel
from common import SEED
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2] / 'scripts'))
from query import connect

D = G / 'data'; D.mkdir(exist_ok=True)
B = 1000
F = F.copy()
F['grp'] = F.pi_group.fillna('acc:' + F.accession)          # missing group -> own cluster
PCT = FuncFormatter(lambda v, _: f'{v:.0f}%')
_save = save


def save(fig, name):
    """House save, then check the rendered PNG stays within 183 mm x 150 mm."""
    _save(fig, name)
    from PIL import Image
    path = G / 'blocks' / f'{name}.png' if JOURNAL else G / f'{name}.png'
    w, h = Image.open(path).size
    dpi = 600 if JOURNAL else 300
    wmm, hmm = w / dpi * 25.4, h / dpi * 25.4
    print(f'{name}: {wmm:.1f} x {hmm:.1f} mm')
    assert wmm <= 183.5 and hmm <= 150, f'{name} exceeds 183 x 150 mm'


def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return (2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum())


def top_share(x, q):
    x = np.sort(np.asarray(x, float))[::-1]; k = max(1, int(round(len(x) * q)))
    return x[:k].sum() / x.sum()


def ci(a):
    return np.nanpercentile(a, [2.5, 97.5], axis=0)


def boot_share(grp, num, den, by, reps=B, seed=SEED, return_draws=False):
    """Group bootstrap of sum(num)/sum(den) within each level of `by`. Returns DataFrame est, lo, hi."""
    rng = np.random.default_rng(seed)
    g_codes, g_uni = pd.factorize(grp); b_codes, b_uni = pd.factorize(by, sort=True)
    G_, Y = len(g_uni), len(b_uni)
    N = np.zeros((G_, Y)); Dn = np.zeros((G_, Y))
    np.add.at(N, (g_codes, b_codes), np.asarray(num, float)); np.add.at(Dn, (g_codes, b_codes), np.asarray(den, float))
    est = N.sum(0) / Dn.sum(0)
    bs = np.empty((reps, Y))
    for r in range(reps):
        w = np.bincount(rng.integers(0, G_, G_), minlength=G_)
        bs[r] = (w @ N) / (w @ Dn)
    lo, hi = ci(bs)
    out = pd.DataFrame({'est': est, 'lo': lo, 'hi': hi, 'num': N.sum(0), 'den': Dn.sum(0)}, index=b_uni)
    return (out, pd.DataFrame(bs, columns=b_uni)) if return_draws else out


def label_end(ax, x, y, s, c, **kw):
    ax.text(x, y, s, color=c, fontsize=6.5, va='center', ha='left', **kw)


# ----------------------------------------------------------------------------------------------------------------
# 34. Contributor concentration
# ----------------------------------------------------------------------------------------------------------------
def fig_concentration():
    s = pd.read_csv(A / 'H1-H3-scale-instruments/H1-dataset-scale.csv.gz', usecols=['accession', 'bytes_all'])
    d = F[['accession', 'pi_group', 'release_year']].merge(s, on='accession', how='left')
    assert d.bytes_all.notna().all()
    n_nogroup = int(d.pi_group.isna().sum())
    d = d.dropna(subset=['pi_group'])
    per = {'2013–17': d[d.release_year <= 2017], '2021–25': d[d.release_year >= 2021]}
    rng = np.random.default_rng(SEED)
    rows = []
    curves = {}
    for lab, x in [('All years', d)] + list(per.items()):
        g = x.groupby('pi_group').agg(n=('accession', 'size'), b=('bytes_all', 'sum'))
        for col, what in [('n', 'datasets'), ('b', 'bytes')]:
            v = g[col].to_numpy()
            est = dict(gini=gini(v), top1=top_share(v, .01), top10=top_share(v, .10))
            bs = {k: [] for k in est}
            if lab != 'All years' or True:
                for _ in range(B):
                    vv = v[rng.integers(0, len(v), len(v))]
                    bs['gini'].append(gini(vv)); bs['top1'].append(top_share(vv, .01)); bs['top10'].append(top_share(vv, .10))
            for k in est:
                lo, hi = ci(np.array(bs[k]))
                rows.append(dict(period=lab, measure=what, stat=k, est=est[k], lo=lo, hi=hi, groups=len(v), datasets=len(x),
                                 total=v.sum()))
            vs = np.sort(v)[::-1]
            curves[(lab, what)] = (np.arange(1, len(vs) + 1) / len(vs) * 100, np.cumsum(vs) / vs.sum() * 100)
    S = pd.DataFrame(rows)
    S.round(5).to_csv(D / '34-contributor-concentration.csv', index=False)
    gall = d.groupby('pi_group').size()
    hist = gall.value_counts().sort_index()
    pd.DataFrame({'datasets_per_group': hist.index, 'groups': hist.values}).to_csv(D / '34-contributor-concentration-histogram.csv', index=False)

    def st(p, m, k):
        r = S[(S.period == p) & (S.measure == m) & (S.stat == k)].iloc[0]; return r
    n_groups = len(gall); n_ds = int(gall.sum())
    one = (gall == 1).sum(); one_ds = one / n_ds
    t1n, t1b = st('All years', 'datasets', 'top1'), st('All years', 'bytes', 'top1')
    k1 = int(round(n_groups * .01))

    # output classes for panel b
    gb = d.groupby('pi_group').agg(n=('accession', 'size'), b=('bytes_all', 'sum'))
    bins = [0, 1, 2, 4, 9, 19, 49, 1e9]; blab = ['1', '2', '3–4', '5–9', '10–19', '20–49', '≥50']
    gb['cls'] = pd.cut(gb.n, bins, labels=blab)
    cls = gb.groupby('cls', observed=False).agg(groups=('n', 'size'), datasets=('n', 'sum'), bytes=('b', 'sum'))
    shares = cls / cls.sum() * 100
    cls.join(shares.add_suffix('_pct')).round(4).to_csv(D / '34-contributor-concentration-classes.csv', index_label='datasets_per_group')
    # sensitivity: 2021-25 subsampled to the 2013-17 number of groups
    rs = np.random.default_rng(SEED)
    e_n = d[d.release_year <= 2017].groupby('pi_group').size().to_numpy(); l_n = d[d.release_year >= 2021].groupby('pi_group').size().to_numpy()
    sub = np.array([top_share(rs.choice(l_n, len(e_n), replace=False), .01) for _ in range(B)])
    pd.DataFrame({'stat': ['top1_datasets_2021_25_subsampled_to_2013_17_groups'], 'mean': [sub.mean()],
                  'lo': [np.percentile(sub, 2.5)], 'hi': [np.percentile(sub, 97.5)], 'groups': [len(e_n)]}).round(4) \
        .to_csv(D / '34-contributor-concentration-sensitivity.csv', index=False)

    ge, gl = st('2013–17', 'datasets', 'gini'), st('2021–25', 'datasets', 'gini')
    fig = plt.figure(figsize=(W, 88 * MM))
    headline(fig, f'The top 1% of contributor groups released {t1n.est*100:.0f}% of PRIDE datasets and {t1b.est*100:.0f}% of its bytes',
             f'{n_ds:,} datasets released 2013–2025 from {n_groups:,} contributor groups (hashed PI e-mail/ORCID; top 1% = {k1} groups; {n_nogroup} without one PI).\n'
             f'Bytes = all files listed per dataset. Gini (datasets) {ge.est:.2f} in 2013–17 and {gl.est:.2f} in 2021–25. Lines in c: 95% bootstrap over groups.')
    top = .755; h = .59; bot = top - h

    # a: top-share curves on log x
    ax = fig.add_axes([.07, bot, .295, h])
    sty = {('2013–17', 'datasets'): (ACC2, .4, 1.2), ('2021–25', 'datasets'): (ACC2, 1, 1.8),
           ('2013–17', 'bytes'): (ACC, .4, 1.2), ('2021–25', 'bytes'): (ACC, 1, 1.8)}
    for key, (c, a, lw) in sty.items():
        xx, yy = curves[key]; ax.plot(xx, yy, color=c, alpha=a, lw=lw, solid_capstyle='butt')
    ax.set_xscale('log'); ax.set_xlim(.07, 100); ax.set_ylim(0, 100); ax.xaxis.set_major_locator(FixedLocator([.1, 1, 10, 100]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}%')); ax.yaxis.set_major_formatter(PCT)
    for q in [1, 10]:
        ax.axvline(q, color=HAIR, lw=.6, zorder=0)
    ax.set_xlabel('Largest x% of contributor groups'); ax.set_ylabel('Share they supplied')
    ax.text(.1, 46, 'Bytes', color=ACC, fontsize=7, fontweight='bold')
    ax.text(.1, 39, 'dark: 2021–25\nlight: 2013–17', color=MUTED, fontsize=6.2, va='top', linespacing=1.1)
    ax.text(3.2, 14, 'Datasets', color=ACC2, fontsize=7, fontweight='bold')
    panel(ax, 'a', x=-.17)
    ax.set_title('Concentration curves', fontsize=7, loc='left', color=INK)

    # b: by output class: share of groups vs datasets vs bytes
    bx = fig.add_axes([.44, bot, .30, h])
    xs = np.arange(len(blab)); wd = .27
    for j, (col, c, lab) in enumerate([('groups', '#a9a6a0', 'Groups'), ('datasets', ACC2, 'Datasets'), ('bytes', ACC, 'Bytes')]):
        bx.bar(xs + (j - 1) * wd, shares[col], width=wd, color=c, lw=0, label=lab)
    for i in [0, len(blab) - 1]:
        for j, col in enumerate(['groups', 'datasets', 'bytes']):
            v = shares[col].iloc[i]
            bx.text(xs[i] + (j - 1) * wd, v + 1.2, f'{v:.0f}' if v >= 1 else f'{v:.1f}', ha='center', va='bottom', fontsize=6, color=INK)
    bx.set_xticks(xs); bx.set_xticklabels(blab); bx.set_xlim(-.6, len(blab) - .4)
    bx.set_ylim(0, 75); bx.yaxis.set_major_formatter(PCT)
    bx.set_xlabel('Datasets released by the group, 2013–2025'); bx.set_ylabel('Share of the total')
    bx.legend(loc='upper right', handlelength=.9, handletextpad=.4, borderaxespad=.3)
    bx.text(.42, 64, f'{one:,} groups deposited once:\n{one/n_groups*100:.0f}% of groups, {one_ds*100:.0f}% of datasets,\n{shares["bytes"].iloc[0]:.0f}% of bytes', fontsize=6.3, color=INK, va='top', linespacing=1.15)
    ge50 = gb[gb.n >= 50]
    bx.text(xs[-1] + .45, 27, f'≥50: {len(ge50)} groups\n(max {gb.n.max()} datasets)', fontsize=6.2, color=MUTED, ha='right', va='bottom', linespacing=1.1)
    panel(bx, 'b', x=-.15)
    bx.set_title('Who supplies what', fontsize=7, loc='left', color=INK)

    # c: dumbbells of top-1% / top-10% shares, early vs late
    cx = fig.add_axes([.815, bot, .165, h])
    items = [('datasets', 'top1', 'Top 1%'), ('datasets', 'top10', 'Top 10%'),
             ('bytes', 'top1', 'Top 1%'), ('bytes', 'top10', 'Top 10%')]
    ypos = [4.2, 3.2, 1.4, .4]
    for (m, k, txt), y in zip(items, ypos):
        c = ACC2 if m == 'datasets' else ACC
        e, l = st('2013–17', m, k), st('2021–25', m, k)
        for r, a, dy in [(e, .4, -.16), (l, 1, .16)]:
            cx.plot([r.lo * 100, r.hi * 100], [y + dy, y + dy], color=c, alpha=a, lw=1.1, zorder=2, solid_capstyle='butt')
            cx.plot(r.est * 100, y + dy, 'o', color=c, alpha=1 if a == 1 else .45, ms=3.8, zorder=3, mec='white', mew=.4)
        if l.hi < .6:
            cx.text(max(e.hi, l.hi) * 100 + 4, y, f'{e.est*100:.0f} → {l.est*100:.0f}%', fontsize=6.3, color=c, fontweight='bold', va='center')
        else:
            cx.text(min(e.lo, l.lo) * 100 - 4, y, f'{e.est*100:.0f} → {l.est*100:.0f}%', fontsize=6.3, color=c, fontweight='bold', va='center', ha='right')
        cx.text(-3, y, txt, fontsize=6.3, color=INK, va='center', ha='right')
    cx.text(50, 4.95, 'Datasets', color=ACC2, fontsize=7, fontweight='bold', ha='center')
    cx.text(50, 2.15, 'Bytes', color=ACC, fontsize=7, fontweight='bold', ha='center')
    cx.set_ylim(-.3, 5.5); cx.set_yticks([]); cx.set_xlim(0, 100); cx.xaxis.set_major_formatter(PCT)
    cx.set_xticks([0, 50, 100])
    cx.set_xlabel('Share supplied')
    cx.set_title('2013–17 (light) → 2021–25 (dark)', fontsize=7, loc='left', color=INK, x=-.32)
    panel(cx, 'c', x=-.33)
    save(fig, '34-contributor-concentration')
    return S


# ----------------------------------------------------------------------------------------------------------------
# 35. International collaboration in linked papers (paper grain)
# ----------------------------------------------------------------------------------------------------------------
def paper_table():
    """One row per OpenAlex work linked (primary link, F.work_id) to a 2013-2025 deposit."""
    c = connect()
    ai = c.sql("select work_id, country_code from lit.openalex_author_institutions where country_code is not null").df()
    au = c.sql("select work_id, countries_json from lit.openalex_authorships").df()
    w = c.sql("""select w.work_id, w.publication_year, w.work_type, e.journal
                 from lit.openalex_works w left join lit.europepmc_works e on w.pmid = e.pmid""").df()
    au['cc'] = au.countries_json.apply(lambda x: json.loads(x) if isinstance(x, str) else [])
    ex = au.explode('cc').dropna(subset=['cc'])[['work_id', 'cc']].rename(columns={'cc': 'country_code'})
    cc = pd.concat([ai, ex]).drop_duplicates()
    nc = cc.groupby('work_id').country_code.nunique().rename('n_countries')
    dep = F.dropna(subset=['work_id'])
    mode = lambda s: s.mode().iloc[0] if s.notna().any() else np.nan
    wk = dep.groupby('work_id').agg(n_deposits=('accession', 'size'), country=('country', mode), grp=('grp', mode),
                                    first_release=('release_year', 'min'))
    P = wk.join(w.set_index('work_id')).join(nc)
    assert w.work_id.is_unique and len(P) == dep.work_id.nunique()
    return P


def fig_collaboration():
    P = paper_table()
    n_all = len(P)
    P = P[P.publication_year.between(2013, 2025)]
    n_yr = len(P)
    P = P[P.n_countries.notna()].copy()
    n_cc = len(P)
    P['intl'] = (P.n_countries >= 2).astype(float); P['tri'] = (P.n_countries >= 3).astype(float)
    a2 = boot_share(P.grp, P.intl, np.ones(len(P)), P.publication_year)
    a3 = boot_share(P.grp, P.tri, np.ones(len(P)), P.publication_year, seed=SEED + 1)
    mc = P.groupby('publication_year').n_countries.mean()
    # country-mix standardised (pooled mix; unresolved country is its own stratum)
    topc = P.country.value_counts().head(20).index
    P['stratum'] = np.where(P.country.isin(topc), P.country, np.where(P.country.isna(), 'Unresolved', 'Other'))
    wts = P.stratum.value_counts(normalize=True)
    r = P.groupby(['publication_year', 'stratum']).intl.mean().unstack()
    std = (r * wts).sum(axis=1) / (r.notna() * wts).sum(axis=1)
    A_ = pd.DataFrame({'papers': a2.den, 'intl_share': a2.est, 'intl_lo': a2.lo, 'intl_hi': a2.hi,
                       'ge3_share': a3.est, 'ge3_lo': a3.lo, 'ge3_hi': a3.hi, 'mean_countries': mc, 'intl_share_country_std': std,
                       'china_share': P.groupby('publication_year').country.apply(lambda s: (s == 'China').mean())})
    A_.index.name = 'paper_year'
    A_.round(4).to_csv(D / '35-collaboration-by-year.csv')
    # by depositing country
    Q = P.dropna(subset=['country'])
    cn = Q.country.value_counts(); keep = cn[cn >= 150].index
    Qk = Q[Q.country.isin(keep)]
    bc = boot_share(Qk.grp, Qk.intl, np.ones(len(Qk)), Qk.country, seed=SEED + 2)
    bc['papers'] = bc.den
    early = Qk[Qk.publication_year <= 2017].groupby('country').intl.agg(['mean', 'size'])
    late = Qk[Qk.publication_year >= 2021].groupby('country').intl.agg(['mean', 'size'])
    bc = bc.join(early.add_prefix('p2013_17_')).join(late.add_prefix('p2021_25_'))
    bc.index.name = 'country'
    bc.round(4).to_csv(D / '35-collaboration-by-country.csv')
    from scipy.stats import spearmanr
    rho = spearmanr(np.log(bc.papers), bc.est).correlation

    EUROPE = {'Germany', 'United Kingdom', 'France', 'Netherlands', 'Spain', 'Switzerland', 'Denmark', 'Italy', 'Sweden', 'Austria',
              'Belgium', 'Czechia', 'Norway', 'Poland', 'Finland', 'Ireland', 'Portugal', 'Hungary', 'Greece', 'Slovenia', 'Estonia'}
    Q['europe'] = Q.country.isin(EUROPE)
    reg = boot_share(Q.grp, Q.intl, np.ones(len(Q)), Q.europe.map({True: 'Europe', False: 'Elsewhere'}), seed=SEED + 3)
    reg.index.name = 'region'; reg.round(4).to_csv(D / '35-collaboration-by-region.csv')
    bc['europe'] = bc.index.isin(EUROPE)
    fig = plt.figure(figsize=(W, 100 * MM))
    ov = P.intl.mean()
    us, cn_ = bc.est['United States'], bc.est['China']
    eu, el = reg.loc['Europe'], reg.loc['Elsewhere']
    headline(fig, f'European groups’ papers are {eu.est*100:.0f}% international; US groups’ {us*100:.0f}%, Chinese groups’ {cn_*100:.0f}%',
             f'Paper grain: {n_cc:,} OpenAlex works published 2013–2025 and linked as the primary paper of a PRIDE dataset (of {n_yr:,}; {n_yr-n_cc} lack author countries).\n'
             f'International = authors’ institutions in ≥2 countries. Country in b = the depositing group’s country. Bands and bars: 95% bootstrap over groups.')
    top = .76; h = .64; bot = top - h
    ax = fig.add_axes([.065, bot, .40, h])
    yrs = a2.index.to_numpy()
    for t, c, lab in [(a2, ACC2, '≥2 countries'), (a3, ACC2, '≥3 countries')]:
        al = 1 if lab.startswith('≥2') else .45
        ax.fill_between(yrs, t.lo * 100, t.hi * 100, color=c, alpha=.15 * al + .03, lw=0)
        ax.plot(yrs, t.est * 100, color=c, lw=1.7 if al == 1 else 1.2, alpha=al)
        ax.text(2025.3, t.est.iloc[-1] * 100, f'{lab}\n{t.est.iloc[-1]*100:.0f}%', color=c, alpha=max(al, .7), fontsize=6.3, va='center',
                fontweight='bold' if al == 1 else 'normal', linespacing=1.05)
    ax.plot(std.index, std * 100, color=INK, lw=.8, ls=(0, (2.5, 1.5)))
    ax.annotate('≥2 countries if the country mix\nhad been the same every year', xy=(2015, std.loc[2015] * 100), xytext=(2016.2, 33),
                fontsize=6.2, color=INK, va='top', linespacing=1.05, arrowprops=dict(arrowstyle='-', color=INK, lw=.5))
    ax.set_xlim(2012.6, 2025.2); ax.set_ylim(0, 60); ax.yaxis.set_major_formatter(PCT)
    ax.set_xticks(range(2013, 2026, 3)); ax.set_xlabel('Paper publication year'); ax.set_ylabel('Share of linked papers')
    ax.set_title(f'International papers: {a2.est.iloc[0]*100:.0f}% (2013) → {a2.est.max()*100:.0f}% ({a2.est.idxmax()}) → {a2.est.iloc[-1]*100:.0f}% (2025)',
                 fontsize=7, loc='left', color=INK)
    ax.text(.03, .04, f'Mean author countries per paper: {mc.iloc[0]:.2f} (2013), {mc.iloc[-1]:.2f} (2025)', transform=ax.transAxes,
            fontsize=6.2, color=MUTED)
    panel(ax, 'a', x=-.12)

    bx = fig.add_axes([.665, bot, .245, h])
    o = bc.sort_values('est')
    ys = np.arange(len(o))
    for y, (cty, rr) in zip(ys, o.iterrows()):
        col = ACC2 if rr.europe else ACC
        bx.plot([rr.lo * 100, rr.hi * 100], [y, y], color=col, lw=1, alpha=.35, solid_capstyle='butt')
        bx.plot(rr.est * 100, y, 'o', color=col, ms=3.6, mec='white', mew=.4, zorder=3)
        bx.text(101, y, f'{int(rr.papers):,}', fontsize=6, color=MUTED, va='center', ha='left')
    bx.set_yticks(ys); bx.set_yticklabels(o.index, fontsize=6.3)
    for tl, is_eu in zip(bx.get_yticklabels(), o.europe):
        tl.set_color(INK if is_eu else ACC)
    bx.tick_params(axis='y', length=0)
    bx.set_ylim(-.7, len(o) - .3); bx.set_xlim(0, 100); bx.xaxis.set_major_formatter(PCT); bx.set_xticks([0, 25, 50, 75, 100])
    bx.axvline(ov * 100, color=HAIR, lw=.7, zorder=0)
    bx.text(1.01, 1.005, 'papers', transform=bx.transAxes, fontsize=6, color=MUTED, va='bottom', ha='left')
    bx.set_xlabel('Share international, 2013–2025 (countries ≥150 papers)')
    t1 = bx.text(-.3, 1.035, f'Europe {eu.est*100:.0f}%', transform=bx.transAxes, fontsize=7, color=ACC2, fontweight='bold', va='bottom')
    bx.annotate(f'  elsewhere {el.est*100:.0f}%', xycoords=t1, xy=(1, 0), va='bottom', fontsize=7, color=ACC, fontweight='bold')
    bx.text(ov * 100 + 1, 1.2, f'all\n{ov*100:.0f}%', fontsize=6, color=MUTED, va='center')
    panel(bx, 'b', x=-.31)
    save(fig, '35-collaboration')
    return A_, bc


# ----------------------------------------------------------------------------------------------------------------
# 36. Research fields of linked papers (deposit grain)
# ----------------------------------------------------------------------------------------------------------------
FIELD_SHORT = {'Biochemistry, Genetics and Molecular Biology': 'Biochemistry, genetics & molecular biology',
               'Medicine': 'Medicine', 'Chemistry': 'Chemistry', 'Agricultural and Biological Sciences': 'Agricultural & biological sciences',
               'Immunology and Microbiology': 'Immunology & microbiology'}


def fig_fields():
    n_all = len(F)
    cov = F.groupby('release_year').agg(datasets=('accession', 'size'), linked=('primary_subfield', lambda x: x.notna().sum()))
    cov['share'] = cov.linked / cov.datasets
    d = F.dropna(subset=['primary_subfield']).copy()
    n_lk = len(d)
    one = np.ones(len(d))
    fields = list(FIELD_SHORT)
    FS = {}
    for i, f in enumerate(fields):
        FS[f] = boot_share(d.grp, (d.primary_field == f).astype(float), one, d.release_year, seed=SEED + 10 + i)
    chem_spec = (d[d.primary_field == 'Chemistry'].primary_subfield == 'Spectroscopy').mean()
    # subfield shares early (2013-15) vs late (2023-25)
    d['per'] = np.where(d.release_year <= 2015, 'early', np.where(d.release_year >= 2023, 'late', 'mid'))
    e = d[d.per != 'mid']
    topsf = d.primary_subfield.value_counts().head(14).index
    rows = []
    for i, sf in enumerate(topsf):
        b, bd = boot_share(e.grp, (e.primary_subfield == sf).astype(float), np.ones(len(e)), e.per, seed=SEED + 30 + i, return_draws=True)
        rr_lo, rr_hi = np.percentile(bd['late'] / bd['early'], [2.5, 97.5])
        rows.append(dict(subfield=sf, ratio_lo=rr_lo, ratio_hi=rr_hi, early=b.est['early'], early_lo=b.lo['early'], early_hi=b.hi['early'],
                         late=b.est['late'], late_lo=b.lo['late'], late_hi=b.hi['late'],
                         n_early=b.num['early'], n_late=b.num['late'], den_early=b.den['early'], den_late=b.den['late']))
    SF = pd.DataFrame(rows); SF['change'] = SF.late - SF.early; SF['ratio'] = SF.late / SF.early
    # publication-grain sensitivity (one row per linked paper, first release year)
    pg = d.sort_values('release_year').drop_duplicates('work_id')
    pgs = pd.crosstab(pg.release_year, pg.primary_subfield == 'Spectroscopy', normalize='index')[True]
    pgm = pd.crosstab(pg.release_year, pg.primary_field == 'Medicine', normalize='index')[True]
    out = pd.DataFrame({f'{FIELD_SHORT[f]}_{k}': FS[f][k] for f in fields for k in ['est', 'lo', 'hi']})
    out['spectroscopy_share'] = pd.crosstab(d.release_year, d.primary_subfield == 'Spectroscopy', normalize='index')[True]
    out['spectroscopy_datasets'] = pd.crosstab(d.release_year, d.primary_subfield == 'Spectroscopy')[True]
    out['other_linked_datasets'] = pd.crosstab(d.release_year, d.primary_subfield == 'Spectroscopy')[False]
    out['pubgrain_spectroscopy_share'] = pgs; out['pubgrain_medicine_share'] = pgm
    out = out.join(cov); out.index.name = 'release_year'
    out.round(4).to_csv(D / '36-research-fields-by-year.csv')
    SF.round(4).to_csv(D / '36-research-fields-subfields.csv', index=False)

    sp = out.spectroscopy_share
    med = FS['Medicine'].est
    fig = plt.figure(figsize=(W, 100 * MM))
    headline(fig, f'MS-methods papers fell from {sp.iloc[0]*100:.0f}% to {sp.iloc[-1]*100:.0f}% of linked datasets; medicine rose from '
                  f'{med.iloc[0]*100:.0f}% to {med.iloc[-1]*100:.0f}%',
             f'OpenAlex primary field/subfield of the paper linked to each PRIDE dataset, by release year: {n_lk:,} of {n_all:,} datasets (2013–2025) have a linked paper.\n'
             f'“Methods” = subfield Spectroscopy ({chem_spec*100:.0f}% of Chemistry). Dataset grain; bands and bars: 95% bootstrap over contributor groups.')
    top = .75; h = .50; bot = top - h
    ax = fig.add_axes([.06, bot, .33, h])
    col = {'Biochemistry, Genetics and Molecular Biology': INK, 'Medicine': ACC2, 'Chemistry': ACC,
           'Agricultural and Biological Sciences': '#5a9e6f', 'Immunology and Microbiology': ACC3}
    yrs = np.array(YEARS)
    ends = {}
    for f in fields:
        t = FS[f]; c = col[f]
        ax.fill_between(yrs, t.lo * 100, t.hi * 100, color=c, alpha=.13, lw=0)
        ax.plot(yrs, t.est * 100, color=c, lw=1.7 if f in ('Medicine', 'Chemistry') else 1.1)
        ends[f] = t.est.iloc[-1] * 100
    order = sorted(fields, key=lambda f: ends[f]); pos = {}; prev = -9
    for f in order:
        pos[f] = max(ends[f], prev + 3.6, 2.2); prev = pos[f]
    short = {'Biochemistry, Genetics and Molecular Biology': 'Biochem. & mol. biology', 'Medicine': 'Medicine',
             'Chemistry': 'Chemistry (MS methods)', 'Agricultural and Biological Sciences': 'Agri. & biological sci.',
             'Immunology and Microbiology': 'Immunol. & microbiol.'}
    for f in fields:
        ax.text(2025.4, pos[f], f'{short[f]} {ends[f]:.0f}%', color=col[f], fontsize=6.2, va='center',
                fontweight='bold' if f in ('Medicine', 'Chemistry') else 'normal')
    ax.set_xlim(2012.7, 2025.3); ax.set_ylim(0, 55); ax.yaxis.set_major_formatter(PCT)
    ax.set_xticks(range(2013, 2026, 2)); ax.tick_params(axis='x', labelbottom=False)
    ax.set_ylabel('Share of linked datasets')
    ax.set_title('Primary field of the linked paper', fontsize=7, loc='left', color=INK)
    panel(ax, 'a', x=-.11)
    # coverage strip (denominator)
    cx = fig.add_axes([.06, bot - .17, .33, .11])
    cx.bar(yrs, cov.share * 100, width=.72, color=HAIR, lw=0)
    for y in [2013, 2021, 2025]:
        cx.text(y, cov.share[y] * 100 - 6, f'{cov.share[y]*100:.0f}', ha='center', va='top', fontsize=6, color=INK)
    cx.set_ylim(0, 100); cx.set_yticks([0, 100]); cx.yaxis.set_major_formatter(PCT)
    cx.set_xlim(2012.7, 2025.3); cx.set_xticks(range(2013, 2026, 2))
    cx.text(0, 1.04, 'Datasets with a linked paper (the denominator in a and b)', transform=cx.transAxes, fontsize=6.2, color=MUTED, va='bottom')

    # b: relative change of subfield shares (ratio, log), with the two shares as text
    bx = fig.add_axes([.715, .08, .15, .67])
    o = SF.sort_values('ratio').reset_index(drop=True)
    for i, r in o.iterrows():
        c = ACC if r.ratio_hi < 1 else (ACC2 if r.ratio_lo > 1 else MUTED)
        bx.plot([r.ratio_lo, r.ratio_hi], [i, i], color=c, lw=1, alpha=.45, solid_capstyle='butt')
        bx.plot(r.ratio, i, 'o', ms=3.8, color=c, mec='white', mew=.4, zorder=3)
        bx.text(1.04, i, f'{r.early*100:.1f}% → {r.late*100:.1f}%', transform=bx.get_yaxis_transform(), fontsize=6.1,
                color=c if c != MUTED else INK, va='center', ha='left')
    bx.axvline(1, color=INK, lw=.6, ls=(0, (2, 2)))
    bx.set_yticks(range(len(o))); bx.set_yticklabels(o.subfield.replace({'Pulmonary and Respiratory Medicine': 'Pulmonary & respiratory'}), fontsize=6.2)
    bx.tick_params(axis='y', length=0)
    bx.set_xscale('log'); bx.set_xlim(.2, 4.5)
    bx.xaxis.set_major_locator(FixedLocator([.25, .5, 1, 2, 4])); bx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: {.25: '×¼', .5: '×½'}.get(v, f'×{v:g}')))
    bx.xaxis.set_minor_locator(NullLocator())
    bx.set_ylim(-.7, len(o) - .3); bx.invert_yaxis()
    nE, nL = int(o.den_early.iloc[0]), int(o.den_late.iloc[0])
    bx.set_xlabel(f'Share 2023–25 ÷ share 2013–15 (log)\nn = {nE:,} → {nL:,} linked datasets;\ngrey: 95% interval includes ×1', linespacing=1.3)
    bx.text(-.95, 1.03, 'Top 14 subfields: relative change', transform=bx.transAxes, fontsize=7, color=INK, va='bottom')
    bx.text(1.04, 1.03, '2013–15 → 2023–25', transform=bx.transAxes, fontsize=6.1, color=MUTED, va='bottom')
    panel(bx, 'b', x=-.97)
    save(fig, '36-research-fields')
    return out, SF


# ----------------------------------------------------------------------------------------------------------------
# 37. Journals of linked papers
# ----------------------------------------------------------------------------------------------------------------
JNAME = {'nature communications': 'Nature Communications', 'journal of proteome research': 'Journal of Proteome Research',
         'molecular & cellular proteomics : mcp': 'Molecular & Cellular Proteomics', 'scientific reports': 'Scientific Reports',
         'international journal of molecular sciences': 'Int. J. Molecular Sciences', 'proteomics': 'Proteomics',
         'journal of proteomics': 'Journal of Proteomics', 'the journal of biological chemistry': 'J. Biological Chemistry',
         'cell reports': 'Cell Reports', 'nucleic acids research': 'Nucleic Acids Research', 'plos one': 'PLOS ONE',
         'proceedings of the national academy of sciences of the united states of america': 'PNAS',
         'analytical chemistry': 'Analytical Chemistry', 'nature': 'Nature', 'elife': 'eLife',
         'frontiers in microbiology': 'Frontiers in Microbiology', 'frontiers in plant science': 'Frontiers in Plant Science',
         'iscience': 'iScience', 'molecular cell': 'Molecular Cell', 'communications biology': 'Communications Biology'}


def fig_journals():
    P = paper_table()
    n_link = len(P)
    P = P[P.publication_year.between(2013, 2025)].copy()
    n_yr = len(P)
    P['jn'] = P.journal.str.lower().str.strip()
    P = P[P.jn.notna()].copy()
    n_j = len(P)
    P['prot'] = P.jn.str.contains('proteom').astype(float)
    # deposit grain ranking: linked datasets per journal (2013-2025 releases with a 2013-2025 paper)
    dep = F.dropna(subset=['work_id']).merge(P[['jn']], left_on='work_id', right_index=True)
    dsn = dep.groupby('jn').size().sort_values(ascending=False)
    top15 = dsn.head(15).index
    yrs = np.arange(2013, 2026)
    tot = P.groupby('publication_year').size().reindex(yrs)
    ct = pd.crosstab(P.publication_year, P.jn).reindex(yrs).fillna(0)
    sh = ct[top15].div(tot, axis=0)
    # concentration: share of papers in each year's own top-10 journals (group bootstrap), proteomics-titled share
    rng = np.random.default_rng(SEED)
    g_codes, g_uni = pd.factorize(P.grp)
    yr_idx = P.publication_year.to_numpy() - 2013
    j_codes, j_uni = pd.factorize(P.jn)
    def top10_by_year(w):
        out = np.empty(len(yrs))
        ww = w[g_codes]
        for k in range(len(yrs)):
            m = yr_idx == k
            c = np.bincount(j_codes[m], weights=ww[m], minlength=len(j_uni))
            out[k] = np.sort(c)[::-1][:10].sum() / c.sum()
        return out
    est10 = top10_by_year(np.ones(len(g_uni)))
    bs = np.array([top10_by_year(np.bincount(rng.integers(0, len(g_uni), len(g_uni)), minlength=len(g_uni)).astype(float)) for _ in range(B)])
    lo10, hi10 = ci(bs)
    pr = boot_share(P.grp, P.prot, np.ones(len(P)), P.publication_year, seed=SEED + 5)
    neff = P.groupby('publication_year').jn.apply(lambda x: 1 / ((x.value_counts(normalize=True)) ** 2).sum())
    nj = P.groupby('publication_year').jn.nunique()
    C = pd.DataFrame({'papers_with_journal': tot, 'top10_share': est10, 'top10_lo': lo10, 'top10_hi': hi10,
                      'proteomics_titled_share': pr.est, 'proteomics_lo': pr.lo, 'proteomics_hi': pr.hi,
                      'distinct_journals': nj, 'effective_journals': neff}, index=yrs)
    C.index.name = 'paper_year'
    C.round(4).to_csv(D / '37-journals-concentration.csv')
    T = sh.T.copy(); T.columns = [f'share_{y}' for y in yrs]
    T.insert(0, 'linked_datasets', dsn[top15]); T.insert(1, 'linked_papers', ct[top15].sum())
    T.insert(0, 'journal', [JNAME.get(j, j.title()) for j in top15]); T['proteomics_titled'] = [('proteom' in j) for j in top15]
    T.round(5).to_csv(D / '37-journals-top15.csv', index=False)
    top10_all = P.jn.value_counts().head(10).sum() / n_j

    fig = plt.figure(figsize=(W, 118 * MM))
    headline(fig, f'Proteomics journals carried {pr.est.iloc[0]*100:.0f}% of dataset-linked papers in 2013 and {pr.est.iloc[-1]*100:.0f}% in 2025',
             f'Journal (Europe PMC, joined by PMID) of the paper linked to PRIDE datasets: {n_j:,} papers published 2013–2025 ({n_yr - n_j} without a\n'
             f'journal, mostly preprints, excluded). Ranked by linked datasets; sparklines = share of each year’s papers (paper grain). Proteomics\n'
             f'journals = title contains “proteom” ({P[P.prot == 1].jn.nunique()} titles). Bands: 95% bootstrap over contributor groups.')
    # a: sparkline table
    x0, x1 = .03, .60
    top = .735; rowh = .66 / 15
    fig.text(x0 + .26, top + .012, 'Linked\ndatasets', fontsize=6.2, color=MUTED, ha='right', va='bottom', linespacing=1.0)
    fig.text(x0 + .285, top + .012, 'Share of each year’s linked papers, 2013 → 2025', fontsize=6.2, color=MUTED, ha='left', va='bottom')
    ymax = sh.max().max()
    for i, j in enumerate(top15):
        yb = top - (i + 1) * rowh
        prot = 'proteom' in j
        c = ACC if prot else ACC2
        fig.text(x0, yb + rowh * .45, JNAME.get(j, j.title()), fontsize=6.5, color=ACC if prot else INK, va='center', fontweight='bold' if prot else 'normal')
        fig.text(x0 + .26, yb + rowh * .45, f'{dsn[j]:,}', fontsize=6.3, color=INK, va='center', ha='right')
        sx = fig.add_axes([x0 + .285, yb + rowh * .1, .21, rowh * .8])
        v = sh[j].to_numpy() * 100
        sx.fill_between(yrs, 0, v, color=c, alpha=.18, lw=0)
        sx.plot(yrs, v, color=c, lw=1.1)
        sx.plot(yrs[[0, -1]], v[[0, -1]], 'o', color=c, ms=2.2)
        sx.set_xlim(2012.5, 2025.5); sx.set_ylim(0, ymax * 100 * 1.05)
        sx.axis('off')
        sx.axhline(0, color=HAIR, lw=.5)
        fig.text(x0 + .505, yb + rowh * .45, f'{v[0]:.1f}% → {v[-1]:.1f}%', fontsize=6.1, color=c, va='center')
        if i < 14:
            fig.add_artist(Line2D([x0, x0 + .6], [yb, yb], color=HAIR, lw=.4))
    fig_letter(fig, x0 - .005, top + .045, 'a', va='bottom', ha='right')
    fig.text(x0 + .005, top + .045, 'Top 15 journals by linked datasets', fontsize=7, color=INK, va='bottom')
    fig.text(x0 + .285, top - 15 * rowh - .015, f'Common vertical scale: 0–{ymax*100:.0f}% of a year’s papers', fontsize=6, color=MUTED, va='top')

    # b: concentration
    bx = fig.add_axes([.72, .08, .26, .655])
    bx.fill_between(yrs, lo10 * 100, hi10 * 100, color=INK, alpha=.12, lw=0)
    bx.plot(yrs, est10 * 100, color=INK, lw=1.6)
    bx.fill_between(yrs, pr.lo * 100, pr.hi * 100, color=ACC, alpha=.15, lw=0)
    bx.plot(yrs, pr.est * 100, color=ACC, lw=1.6)
    bx.text(2025.3, est10[-1] * 100 + 2.5, f'Top 10 journals\nof the year {est10[-1]*100:.0f}%', color=INK, fontsize=6.3, va='bottom', ha='right', fontweight='bold')
    bx.text(2013.3, est10[0] * 100 + 2, f'{est10[0]*100:.0f}%', color=INK, fontsize=6.3, va='bottom')
    bx.text(2025.3, pr.est.iloc[-1] * 100 - 3, f'Proteomics\njournals {pr.est.iloc[-1]*100:.0f}%', color=ACC, fontsize=6.3, va='top', ha='right', fontweight='bold')
    bx.text(2013.3, pr.est.iloc[0] * 100 + 2, f'{pr.est.iloc[0]*100:.0f}%', color=ACC, fontsize=6.3, va='bottom')
    bx.set_xlim(2012.6, 2025.4); bx.set_ylim(0, 80); bx.yaxis.set_major_formatter(PCT); bx.set_xticks(range(2013, 2026, 3))
    bx.set_xlabel('Paper publication year'); bx.set_ylabel('Share of linked papers')
    fig.text(.72, top + .045, f'Concentration (journals per year: {nj.iloc[0]} → {nj.iloc[-1]})',
             fontsize=7, color=INK, va='bottom')
    fig_letter(fig, .69, top + .045, 'b', va='bottom', ha='right')
    save(fig, '37-journals')
    return C, T


if __name__ == '__main__':
    which = sys.argv[1:] or ['34', '35', '36', '37']
    if '34' in which:
        S = fig_concentration(); print(S.round(3).to_string())
    if '35' in which:
        a, b = fig_collaboration(); print(a.round(3).to_string()); print(b.round(3).to_string())
    if '36' in which:
        a, b = fig_fields(); print(a.round(3).to_string()); print(b.round(3).to_string())
    if '37' in which:
        a, b = fig_journals(); print(a.round(3).to_string()); print(b.iloc[:, :4].to_string())
