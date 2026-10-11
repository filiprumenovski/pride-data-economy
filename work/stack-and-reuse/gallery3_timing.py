"""Gallery 3, figures 38-41: timing and regional arrival (descriptive only).

38 38-instrument-arrival   first release year of each newer instrument model in deposits from each continent,
                           observed vs a volume-matched null (ledger H3.6).
39 39-data-age-by-country  age of the data at submission (latest raw-file-name date -> submission) by country,
                           all vs like-for-like, and by instrument model (ledger H2.2, H2.5-H2.7).
40 40-weekly-rhythm        the clock behind the recorded submission timestamp (UTC hour by depositor country,
                           2014-19 vs 2023-25), the local working week since the switch, and the lab clock
                           read from raw-file names (ledger H2.1, H2.14; new: time of day).
41 41-paper-and-release    release date minus linked-paper date, bulk release days vs other days, and the
                           composition of each release year (ledger H2.9-H2.10; Attacks/release-timing).

Inputs (frozen): Arsenal/H1-H3-scale-instruments/H3-region-*.csv, Arsenal/H2-time/Dataset-time.csv.gz and
Validity-acquisition-weekday.csv, Arsenal/H11-clocks/Dataset-clocks.csv.gz, Deposit-features.parquet,
pride_files.submission_date (read-only through scripts/query.py; time of day, UTC as stored).
Uncertainty: 95% percentile interval over 2,000 resamples of whole pi_group clusters (datasets without a group
hash are their own cluster), seed 20261001. Plotted tables -> Gallery/data/<NN>-<slug>.csv.
Not claimed: deposits are not ownership or adoption records; first appearance is not purchase; timestamps
describe when a record was created, not why.
Run from work/stack-and-reuse/: python gallery3_timing.py [38 39 40 41]
"""
import sys, json
from zoneinfo import ZoneInfo
import numpy as np, pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3
from matplotlib.ticker import FuncFormatter
from common import R, O, SEED

D = G / 'data'; D.mkdir(exist_ok=True)
REPS = 2000
H3 = A / 'H1-H3-scale-instruments'
H2 = A / 'H2-time'
CONTS = ['Europe', 'North America', 'Asia', 'Oceania', 'South America', 'Africa']
WD = ['Monday', 'Tuesday', 'Wednesday', 'Thursday', 'Friday', 'Saturday', 'Sunday']


# ------------------------------------------------------------------ helpers ----------------------------------
def boot_weights(groups, reps=REPS, seed=SEED):
    """Per-replicate dataset weights from resampling whole clusters (multinomial cluster counts)."""
    codes, uniq = pd.factorize(pd.Series(groups).astype(str))
    rng = np.random.default_rng(seed); k = len(uniq)
    cnt = np.stack([np.bincount(rng.integers(0, k, k), minlength=k) for _ in range(reps)]).astype(np.int32)
    return cnt[:, codes]


def gb_median(v, g, q=.5):
    """Median and 95% group-bootstrap interval."""
    v = np.asarray(v, float); o = np.argsort(v); v = v[o]
    Wt = boot_weights(np.asarray(g)[o]); cw = np.cumsum(Wt, axis=1)
    b = v[np.argmax(cw >= q * cw[:, -1:], axis=1)]
    return float(np.median(v)), float(np.quantile(b, .025)), float(np.quantile(b, .975))


def gb_share(flag, g):
    f = np.asarray(flag, float); Wt = boot_weights(g)
    b = (Wt @ f) / Wt.sum(axis=1)
    return float(f.mean()), float(np.quantile(b, .025)), float(np.quantile(b, .975))


def continents():
    geo = json.loads((R / 'work/stack-and-reuse/ne_50m_admin_0_countries.geojson').read_text())
    cont, iso = {}, {}
    for f in geo['features']:
        p = f['properties']
        for k in ('ADMIN', 'NAME', 'NAME_LONG'):
            cont[p[k]] = p['CONTINENT']; iso[p[k]] = p['ISO_A2_EH']
    return cont, iso


def short(m):
    return (m.replace('Orbitrap Fusion Lumos', 'Fusion Lumos').replace('Orbitrap Fusion ETD', 'Fusion ETD')
             .replace('Orbitrap Exploris', 'Exploris').replace('Orbitrap Eclipse', 'Eclipse').replace('Orbitrap Astral', 'Astral'))


# ================================================================== 38 instrument arrival ===================
def fig38():
    r = pd.read_csv(H3 / 'H3-region-first-appearance.csv')
    sm = pd.read_csv(H3 / 'H3-region-summary.csv').set_index('continent')
    life = pd.read_csv(H3 / 'H3-model-life-cycles.csv').set_index('model')
    r['vendor'] = r.model.map(life.vendor); r['datasets_total'] = r.model.map(life.datasets)
    r['exp_lag'] = r.expected_median_first_year - r.global_first_year
    r['censor_lag'] = 2025 - r.global_first_year
    r['status'] = np.select([r.first_year.isna() & (r.p_late < .05), r.first_year.isna(),
                             r.p_late < .05, r.p_early < .05],
                            ['not seen, expected', 'not seen', 'later than null', 'earlier than null'], 'as expected')
    order = (r.drop_duplicates('model').sort_values(['global_first_year', 'datasets_total'], ascending=[True, False]).model.tolist())
    r.to_csv(D / '38-instrument-arrival.csv', index=False)
    # summary numbers (computed, not typed)
    med_lag = r.groupby('continent').lag_vs_global.median()
    med_exp = r.groupby('continent').exp_lag.median()
    n_late = int((r.p_late < .05).sum()); n_pairs = len(r)
    seen = r.groupby('continent').first_year.apply(lambda s: s.notna().sum())
    bb = r.dropna(subset=['lag_vs_global', 'exp_lag'])
    both = bb.groupby('continent').size(); late = (bb.lag_vs_global > bb.exp_lag).groupby(bb.continent).sum()
    assert med_lag['Africa'] == sm.loc['Africa', 'median_lag_vs_global'] and med_lag['Asia'] == 1

    fig = plt.figure(figsize=(W, 112 * MM))
    headline(fig, f'Africa\'s {med_lag["Africa"]:.0f}-year lag in new instrument models is about what its '
                  f'{sm.loc["Africa", "continent_datasets_2013_25"]} datasets predict',
             f'First release year of {r.model.nunique()} instrument models (≥100 datasets; first ≥5 in 2014 or later) in deposits from each continent\n'
             f'(PI country; {int(sm.continent_datasets_2013_25.sum()):,} of 36,249 datasets), in years after its first appearance anywhere; lag = median over models seen.\nRing: median first year '
             'expected if labels were independent of continent within each year (volume-matched null). A deposit is not a purchase record.')
    L, Rr, B, T = .15, .995, .085, .725
    gap = .009; pw = (Rr - L - gap * 5) / 6
    ny = len(order); ypos = {m: ny - 1 - i for i, m in enumerate(order)}
    XMAX = 12.3
    col = {'later than null': ACC, 'earlier than null': ACC2, 'as expected': INK}
    axes = []
    for j, ct in enumerate(CONTS):
        ax = fig.add_axes([L + j * (pw + gap), B, pw, T - B]); axes.append(ax)
        x = r[r.continent == ct]
        for _, q in x.iterrows():
            y = ypos[q.model]
            ax.add_patch(plt.Rectangle((q.censor_lag + .5, y - .5), XMAX, 1, color='#efeeea', lw=0, zorder=0))
            if not np.isnan(q.exp_lag):
                ax.plot(q.exp_lag, y, 'o', ms=4.4, mfc='white', mec=MUTED, mew=.8, zorder=2)
            if np.isnan(q.first_year):
                c = ACC if q.status == 'not seen, expected' else MUTED
                ax.text(q.censor_lag + .55, y, '›', color=c, fontsize=7.5, va='center', ha='left', fontweight='bold', zorder=3)
                continue
            c = col[q.status]
            if not np.isnan(q.exp_lag) and q.lag_vs_global != q.exp_lag:
                ax.plot([q.exp_lag, q.lag_vs_global], [y, y], color=c if c != INK else MUTED, lw=.9, zorder=1)
            ax.plot(q.lag_vs_global, y, 'o', ms=3.4, color=c, mec='white', mew=.4, zorder=3)
        ax.set_xlim(-.7, XMAX); ax.set_ylim(-.7, ny - .3)
        ax.set_xticks([0, 5, 10]); ax.set_yticks([])
        ax.tick_params(axis='x', length=2)
        for v in [0, 5, 10]:
            ax.axvline(v, color=HAIR, lw=.4, zorder=0)
        nd = int(sm.loc[ct, 'continent_datasets_2013_25'])
        ax.text(0, 1.105, ct, transform=ax.transAxes, fontsize=7, fontweight='bold', color=INK, va='bottom')
        ax.text(0, 1.045, f'{nd:,} datasets\nlag {med_lag[ct]:.0f} y, expected {med_exp[ct]:.0f} y',
                transform=ax.transAxes, fontsize=6, color=MUTED, va='bottom', linespacing=1.05)
        ax.text(0, 1.008, f'later than ring: {late[ct]}/{both[ct]}', transform=ax.transAxes, fontsize=6, va='bottom',
                color=ACC if late[ct] / both[ct] >= 2 / 3 else MUTED, fontweight='bold' if late[ct] / both[ct] >= 2 / 3 else 'normal')
        if j == 0:
            ax.set_yticks(range(ny)); ax.set_yticklabels([short(m) for m in order[::-1]], fontsize=6.5)
            ax.tick_params(axis='y', length=0, pad=19)
            for m in order:
                g0 = int(r.loc[r.model == m, 'global_first_year'].iloc[0])
                ax.text(-.035, ypos[m], f'{g0}', transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6, color=MUTED)
            ax.text(-.035, ny - .45, 'first\nanywhere', transform=ax.get_yaxis_transform(), ha='right', va='bottom', fontsize=6, color=MUTED, linespacing=1.0)
    fig.text(L + 2.5 * (pw + gap), .03, 'Years after the model first appears in any PRIDE deposit (release year)',
             ha='center', fontsize=7, color=INK)
    # legend
    hx = [Line2D([], [], marker='o', ls='', color=INK, ms=3.4, label='observed first year'),
          Line2D([], [], marker='o', ls='', mfc='white', mec=MUTED, ms=4.4, label='expected (volume-matched)'),
          Line2D([], [], marker='o', ls='-', color=ACC, ms=3.4, lw=.9, label=f'later than expected, p<.05† ({n_late} of {n_pairs} pairs)'),
          Line2D([], [], marker='o', ls='-', color=ACC2, ms=3.4, lw=.9, label='earlier, p<.05†'),
          Line2D([], [], marker='$›$', ls='', color=MUTED, ms=5, label='not seen by 2025'),
          plt.Rectangle((0, 0), 1, 1, color='#efeeea', label='after 2025')]
    fig.legend(handles=hx, loc='lower left', bbox_to_anchor=(.005, -.045), ncol=6, fontsize=6, handletextpad=.3, columnspacing=1.1)
    footnote(fig, .005, -.075, '† Hypergeometric null treating datasets as independent; it ignores contributor clustering, so p-values are too small '
             '(≈4.5 of 90 pairs expected below .05 by chance).', fontsize=6, color=MUTED)
    save(fig, '38-instrument-arrival')
    print('38', med_lag.to_dict(), med_exp.to_dict(), n_late, n_pairs, late.to_dict(), both.to_dict())


# ================================================================== 39 data age by country ==================
def age_rows(sub, by, levels=None, minn=1):
    rows = []
    for k, x in sub.groupby(by):
        if (levels is not None and k not in levels) or len(x) < minn:
            continue
        m, lo, hi = gb_median(x.age_latest, x.grp)
        rows.append({by: k, 'datasets': len(x), 'groups': x.grp.nunique(), 'median_days': m, 'ci_lo': lo, 'ci_hi': hi,
                     'p25': x.age_latest.quantile(.25), 'p75': x.age_latest.quantile(.75),
                     'group_weighted_median': float(np.median(x.groupby('grp').age_latest.median()))})
    return pd.DataFrame(rows)


def fig39():
    t = pd.read_csv(H2 / 'Dataset-time.csv.gz')
    f = F[['accession', 'country', 'vendor_single', 'dia_expanded', 'instruments', 'human']]
    d = t.merge(f, on='accession', validate='one_to_one')
    d['grp'] = d.pi_group.fillna('acc:' + d.accession)
    a = d[d.analysed].copy()
    n_raw = int((d.n_raw > 0).sum()); n_an = len(a)
    assert n_an == 11868 and a.age_latest.median() == 492
    a['single'] = np.where(a.instruments.fillna('').str.contains(r'\|'), None, a.instruments.fillna('').str.strip())
    ll = a[(~a.dia_expanded) & (a.vendor_single == 'Thermo Fisher') & a.sub_year.between(2016, 2025)]
    cn = ll.country.value_counts(); keep = cn[cn >= 100].index.tolist()
    C = age_rows(ll, 'country', keep).rename(columns=lambda c: c if c == 'country' else 'll_' + c)
    Ca = age_rows(a, 'country', keep).rename(columns=lambda c: c if c == 'country' else 'all_' + c)
    cov = d[d.n_raw > 0].groupby('country').agg(raw_datasets=('analysed', 'size'), coverage=('analysed', 'mean'))
    C = C.merge(Ca, on='country').merge(cov, left_on='country', right_index=True).sort_values('ll_median_days')
    llm, lll, llh = gb_median(ll.age_latest, ll.grp)
    C['vs_pooled'] = np.where(C.ll_ci_lo > llm, 'older', np.where(C.ll_ci_hi < llm, 'younger', 'overlaps'))
    # instrument models (single-model datasets submitted 2020-2025, >= 60) and DIA within model
    im = a[a.single.notna() & (a.single != '') & a.sub_year.between(2020, 2025)]
    big = im.single.value_counts(); big = big[big >= 60].index
    M = age_rows(im, 'single', set(big)).sort_values('median_days')
    life = pd.read_csv(H3 / 'H3-model-life-cycles.csv').set_index('model')
    M['vendor'] = M.single.map(life.vendor); M['first_year_any'] = M.single.map(life.first_year_any)
    im['dia'] = np.where(im.dia_expanded, 'DIA', 'not DIA')
    dm = []
    for mdl, x in im[im.single.isin(big)].groupby('single'):
        if (x.dia == 'DIA').sum() >= 25 and (x.dia == 'not DIA').sum() >= 25:
            r_ = age_rows(x, 'dia'); r_['single'] = mdl; dm.append(r_)
    a['dia'] = np.where(a.dia_expanded, 'DIA', 'not DIA')
    r_ = age_rows(a, 'dia'); r_['single'] = 'All analysed datasets'; dm.append(r_)
    DM = pd.concat(dm)
    if FROZEN_DATA:
        C, M, DM = [pd.read_csv(D / f'39-data-age-by-country{suffix}.csv')
                    for suffix in ['', '-models', '-dia']]
        C['coverage'] = C.all_datasets / C.raw_datasets
        # Original block's pooled interval; the earlier country tables omitted it.
        pooled = pd.read_csv(D / '39-data-age-by-country-pooled.csv').iloc[0]
        llm, lll, llh = pooled[['median_days', 'ci_lo', 'ci_hi']]
    C.round(3).to_csv(D / '39-data-age-by-country.csv', index=False)
    M.round(3).to_csv(D / '39-data-age-by-country-models.csv', index=False)
    DM.round(3).to_csv(D / '39-data-age-by-country-dia.csv', index=False)
    # sensitivity: human only, like-for-like (printed, not drawn)
    hs = ll[ll.human == True].groupby('country').age_latest.agg(['median', 'size'])
    print('39 human-only like-for-like:', hs[hs['size'] >= 80].sort_values('median').to_dict())
    old = C[C.vs_pooled == 'older']; yng = C[C.vs_pooled == 'younger']
    names = ', '.join(old.country.str.replace('United Kingdom', 'UK').tolist()[:-1]) + ' and ' + old.country.str.replace('United Kingdom', 'UK').iloc[-1]
    fig = plt.figure(figsize=(W, 118 * MM))
    headline(fig, f'Like for like, data from the {names} are {old.ll_median_days.min() / 365.25:.1f}–{old.ll_median_days.max() / 365.25:.1f} '
                  f'years old at deposit, against {llm / 365.25:.1f} pooled',
             f'Days from the latest acquisition date in raw-file names to PRIDE submission; only {n_an:,} of {n_raw:,} datasets with raw files '
             f'({n_an / n_raw:.0%}) have dated names.\n(a) Thermo-only, not DIA, submitted 2016–2025, countries with ≥100 such datasets ({len(ll):,}; '
             f'the pooled median covers all of them). No country is clearly\nbelow it. Dots: median; bars: 95% interval over 2,000 resamples of '
             'contributor groups. A file-name date is a lab label (possibly a naming, not a run, date).')
    ax = fig.add_axes([.12, .09, .375, .67])
    n = len(C); y = np.arange(n)[::-1]
    ax.axvspan(lll, llh, color=HAIR, alpha=.45, lw=0)
    ax.axvline(llm, color=MUTED, lw=.7)
    ax.text(llm, n - .35, f'pooled {llm:.0f} d', ha='center', va='bottom', fontsize=6, color=MUTED,
            bbox=dict(boxstyle='square,pad=.1', fc='white', ec='none'))
    ax.axvline(365.25, color=HAIR, lw=.6, ls=(0, (3, 2))); ax.axvline(730.5, color=HAIR, lw=.6, ls=(0, (3, 2)))
    ax.text(365.25, n - .35, '1 year', ha='center', va='bottom', fontsize=6, color=MUTED)
    ax.text(730.5, n - .35, '2 years', ha='center', va='bottom', fontsize=6, color=MUTED)
    labs = []
    for yy, (_, q) in zip(y, C.iterrows()):
        hot = q.vs_pooled != 'overlaps'
        cc = INK if hot else '#8f8d88'
        ax.plot([q.ll_ci_lo, q.ll_ci_hi], [yy, yy], color=cc, lw=1.4 if hot else 1.0, solid_capstyle='butt')
        ax.plot(q.ll_median_days, yy, 'o', ms=4.4 if hot else 3.8, color=cc, mec='white', mew=.5, zorder=3)
        ax.plot(q.all_median_days, yy, 'D', ms=2.6, mfc='white', mec=MUTED, mew=.7, zorder=2)
        for xx, txt, col_ in [(1.135, f'{q.ll_datasets:,}', INK), (1.245, f'{q.ll_groups}', INK),
                              (1.36, f'{q.coverage:.0%}', ACC if q.coverage < .25 else INK)]:
            ax.text(xx, yy, txt, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=6, color=col_)
    for xx, lab in [(1.135, 'datasets'), (1.245, 'groups'), (1.36, 'dated*')]:
        ax.text(xx, n - .35, lab, transform=ax.get_yaxis_transform(), ha='right', va='bottom', fontsize=6, color=MUTED)
    ax.set_yticks(y); ax.set_yticklabels(C.country.str.replace('United Kingdom', 'UK').str.replace('United States', 'USA'))
    for tl, (_, q) in zip(ax.get_yticklabels(), C.iterrows()):
        if q.vs_pooled != 'overlaps':
            tl.set_fontweight('bold')
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-.7, n + .3); ax.set_xlim(250, 900)
    ax.set_xlabel('Median age of data at submission (days)')
    ax.set_title('By PI country', fontsize=7, loc='left', color=INK, fontweight='bold')
    hx = [Line2D([], [], marker='o', ls='-', color=INK, ms=4, lw=1.3, label='interval above pooled'),
          Line2D([], [], marker='o', ls='-', color='#8f8d88', ms=3.6, lw=1, label='overlaps pooled'),
          Line2D([], [], marker='D', ls='', mfc='white', mec=MUTED, ms=2.8, label='all datasets of the country')]
    ax.legend(handles=hx, loc='upper right', bbox_to_anchor=(1, .955), fontsize=6, handletextpad=.3, borderaxespad=.3, labelspacing=.3)
    ax.text(0, -.115, '* share of the country\'s raw-file datasets with dated names; red <25%', transform=ax.transAxes,
            fontsize=6, color=MUTED, va='top')
    panel(ax, 'a', x=-.2, y=1.0)
    # b: instrument models
    X0, XW = .755, .235
    bx = fig.add_axes([X0, .42, XW, .34])
    nm = len(M); ym = np.arange(nm)[::-1]
    for yy, (_, q) in zip(ym, M.iterrows()):
        cc = VEND.get(q.vendor, MUTED)
        bx.plot([q.ci_lo, q.ci_hi], [yy, yy], color=cc, lw=1.2, alpha=.85, solid_capstyle='butt')
        bx.plot(q.median_days, yy, 'o', ms=3.4, color=cc, mec='white', mew=.4, zorder=3)
    bx.set_yticks(ym); bx.set_yticklabels([short(m).replace('LTQ Orbitrap', 'LTQ Orb.') for m in M.single], fontsize=6)
    bx.tick_params(axis='y', length=0, pad=2)
    bx.set_ylim(-.7, nm - .3); bx.set_xlim(0, 1300); bx.set_xticks([0, 365.25, 730.5, 1095.75])
    bx.set_xticklabels(['0', '1 y', '2 y', '3 y'])
    for v_ in [365.25, 730.5]:
        bx.axvline(v_, color=HAIR, lw=.5, zorder=0)
    bx.set_title('By instrument model (single-model\ndatasets submitted 2020–2025)', fontsize=7, loc='left', color=INK, fontweight='bold', linespacing=1)
    for v, cc in VEND.items():
        if v in M.vendor.values:
            bx.plot([], [], 'o', color=cc, ms=3.4, label=v.replace(' Fisher', ''))
    bx.legend(loc='upper right', fontsize=6, handletextpad=.1, borderaxespad=.2, labelspacing=.2)
    bx.text(.03, .03, 'newer models cannot\nhave old data', transform=bx.transAxes, ha='left', va='bottom', fontsize=6, color=MUTED, style='italic')
    panel(bx, 'b', x=-.48, y=1.04)
    # c: DIA vs not, within model
    cx = fig.add_axes([X0, .09, XW, .225])
    order = ['All analysed datasets'] + [m for m in M.single if m in set(DM.single)]
    yc = np.arange(len(order))[::-1]
    for yy, mdl in zip(yc, order):
        z = DM[DM.single == mdl].set_index('dia')
        cx.plot([z.median_days['DIA'], z.median_days['not DIA']], [yy, yy], color=HAIR, lw=1.5, zorder=1)
        cx.plot(z.median_days['DIA'], yy, 'o', ms=3.6, color=INK, mec='white', mew=.4, zorder=3)
        cx.plot(z.median_days['not DIA'], yy, 'o', ms=3.4, mfc='white', mec=INK, mew=.8, zorder=3)
    cx.set_yticks(yc); cx.set_yticklabels(['All, 2013–25'] + [short(m) for m in order[1:]], fontsize=6)
    cx.tick_params(axis='y', length=0, pad=2)
    cx.set_ylim(-.6, len(order) - .4); cx.set_xlim(0, 1300); cx.set_xticks([0, 365.25, 730.5, 1095.75]); cx.set_xticklabels(['0', '1 y', '2 y', '3 y'])
    for v_ in [365.25, 730.5]:
        cx.axvline(v_, color=HAIR, lw=.5, zorder=0)
    cx.set_title('DIA vs not, within model (≥25 each)', fontsize=7, loc='left', color=INK, fontweight='bold')
    cx.plot([], [], 'o', color=INK, ms=3.6, label='DIA'); cx.plot([], [], 'o', mfc='white', mec=INK, ms=3.4, label='not DIA')
    cx.legend(loc='upper right', fontsize=6, handletextpad=.1, borderaxespad=.2, labelspacing=.2)
    cx.set_xlabel('Median age of data at submission')
    panel(cx, 'c', x=-.48, y=1.04)
    save(fig, '39-data-age-by-country')
    print('39', C[['country', 'll_median_days', 'll_ci_lo', 'll_ci_hi', 'all_median_days', 'coverage', 'vs_pooled', 'll_group_weighted_median']].round(2).to_string())
    print(llm, lll, llh); print(M[['single', 'median_days', 'ci_lo', 'ci_hi', 'datasets']].to_string()); print(DM.to_string())


# ================================================================== 40 weekly rhythm ========================
# Rows for panel a, ordered east -> west; local 09:00-17:00 window in UTC hours (standard time; envelope for
# countries spanning several zones). Used only to draw reference boxes; the shares are computed from UTC stamps.
ROWS40 = [('Australia', 'AU', 23, 9, 'AWST–AEST'), ('China', 'CN', 1, 9, 'UTC+8'), ('India', 'IN', 3.5, 11.5, 'UTC+5:30'),
          ('Israel', 'IL', 7, 15, 'UTC+2'), ('Germany', 'DE', 8, 16, 'UTC+1'), ('United Kingdom', 'GB', 9, 17, 'UTC+0'),
          ('Brazil', 'BR', 12, 20, 'UTC−3'), ('United States', 'US', 14, 25, 'ET–PT')]
MULTI_TZ = {'US', 'CA', 'AU', 'BR', 'RU', 'MX', 'ID', 'KZ', 'MN', 'CD', 'CL', 'AR'}


def stamps():
    """One submission timestamp per dataset (all files of a dataset share it; checked)."""
    sys.path.insert(0, str(R / 'scripts')); from query import connect
    c = connect()
    t = c.execute("""select project_accession accession, min(submission_date) ts, max(submission_date) ts_max,
        bool_and(strftime(submission_date, '%H:%M:%S') = '00:00:00') midnight
        from pride_files where project_accession in (select accession from study) group by 1""").df()
    assert (t.ts == t.ts_max).all() and len(t) == 36249
    t['ts'] = pd.to_datetime(t.ts, utc=True)
    t['midnight'] = t.ts.eq(t.ts.dt.normalize())
    f = F[['accession', 'country', 'pi_group', 'submission_date', 'lag_days']]
    t = t.merge(f, on='accession', validate='one_to_one')
    assert (t.ts.dt.tz_localize(None).dt.normalize() == t.submission_date).all()
    t['grp'] = t.pi_group.fillna('acc:' + t.accession)
    cont, iso = continents()
    t['continent'] = t.country.map(cont); t['iso'] = t.country.map(iso)
    zt = {}
    for line in open('/usr/share/zoneinfo/zone.tab'):
        if not line.startswith('#'):
            c_, _, tz = line.rstrip('\n').split('\t')[:3]; zt.setdefault(c_, []).append(tz)
    t['tz'] = t.iso.map(lambda c_: zt[c_][0] if c_ in zt and c_ not in MULTI_TZ else None)
    t['year'] = t.ts.dt.year; t['h'] = t.ts.dt.hour
    return t


def fig40():
    t = stamps()
    n_all = len(t); x = t[~t.midnight].copy(); n_time = len(x)
    x['office'] = x.h.between(8, 16)
    early = x[x.year.between(2014, 2019)]; late = x[x.year.between(2023, 2025)]
    sh_early = early.office.mean(); sh_2025 = x[x.year == 2025].office.mean()
    # a: UTC hour by country, two periods
    H = []
    for per, z in [('2014–2019', early), ('2023–2025', late)]:
        for name, *_ in ROWS40:
            q = z[z.country == name]
            v = q.h.value_counts(normalize=True).reindex(range(24)).fillna(0)
            for h_, sh in v.items():
                H.append({'period': per, 'country': name, 'utc_hour': h_, 'share': sh, 'datasets': len(q)})
    H = pd.DataFrame(H)
    # b: share stamped 08:00-16:59 UTC by year and region
    reg = {'Europe': 'Europe', 'North America': 'Americas', 'South America': 'Americas', 'Asia': 'Asia–Oceania', 'Oceania': 'Asia–Oceania'}
    x['region'] = x.continent.map(reg)
    Bt = []
    for (rg, yr), z in x[x.region.notna() & x.year.between(2014, 2025)].groupby(['region', 'year']):
        sh, lo, hi = gb_share(z.office, z.grp)
        Bt.append({'region': rg, 'year': yr, 'datasets': len(z), 'share_08_17_utc': sh, 'lo': lo, 'hi': hi})
    Bt = pd.DataFrame(Bt)
    if FROZEN_DATA:
        Bt = pd.read_csv(D / '40-weekly-rhythm-office-share.csv')
    # c: weekday shares under three clocks
    sing = x[x.tz.notna()].copy()
    sing['lwd'] = [ts.tz_convert(ZoneInfo(z)).day_name() for ts, z in zip(sing.ts, sing.tz)]
    e_ = x[x.year.between(2014, 2019)]; l_ = sing[sing.year.between(2023, 2025)]
    acq = pd.read_csv(H2 / 'Validity-acquisition-weekday.csv').set_index('subset')
    Cw = pd.DataFrame({'submission stamp 2014–2019 (UTC)': e_.ts.dt.day_name().value_counts(normalize=True).reindex(WD),
                       'submission stamp 2023–2025 (local)': l_.lwd.value_counts(normalize=True).reindex(WD),
                       'acquisition date in raw-file names': acq.loc['all', WD].astype(float)})
    il_e = e_[e_.country == 'Israel']; il_l = l_[l_.country == 'Israel']
    Iw = pd.DataFrame({'2014–2019 stamp (UTC)': il_e.ts.dt.day_name().value_counts(normalize=True).reindex(WD).fillna(0),
                       '2023–2025 stamp (local)': il_l.lwd.value_counts(normalize=True).reindex(WD).fillna(0),
                       'acquisition (file names)': acq.loc['Israel', WD].astype(float)})
    wk = lambda v: v['Saturday'] + v['Sunday']
    # December, truncation-safe (submission years 2014-2023; later years miss still-private datasets)
    mz = t[t.year.between(2014, 2023)]; dec = (mz.ts.dt.month == 12)
    dec_sh, dec_lo, dec_hi = gb_share(dec, mz.grp)
    H.to_csv(D / '40-weekly-rhythm.csv', index=False)
    Bt.round(4).to_csv(D / '40-weekly-rhythm-office-share.csv', index=False)
    pd.concat({'all': Cw, 'Israel': Iw}, axis=1).round(4).to_csv(D / '40-weekly-rhythm-weekday.csv')

    fig = plt.figure(figsize=(W, 140 * MM))
    headline(fig, f'Until 2019, {sh_early:.0%} of PRIDE submissions were time-stamped 08–17 UTC, wherever they came from',
             f'Submission timestamp of {n_time:,} of {n_all:,} datasets (the rest carry no time of day), as stored (UTC). From about 2020 the stamps '
             f'follow the\ndepositors\' own working day ({sh_2025:.0%} within 08–17 UTC in 2025). What the stamp records is not documented; weekday '
             'patterns in pre-2020\nsubmission dates describe this single clock, not depositors.')
    cmap = LinearSegmentedColormap.from_list('h', ['#ffffff', '#c9d9e8', '#6f9cc3', '#2e6f9e', '#163a57'])
    vmax = .2
    for j, per in enumerate(['2014–2019', '2023–2025']):
        ax = fig.add_axes([.105 + j * .45, .53, .39, .285])
        M_ = H[H.period == per].pivot(index='country', columns='utc_hour', values='share').reindex([r_[0] for r_ in ROWS40])
        ax.imshow(M_.to_numpy(), aspect='auto', cmap=cmap, vmin=0, vmax=vmax, extent=(0, 24, len(ROWS40), 0), interpolation='nearest')
        for i, (name, iso_, a0, a1, lab) in enumerate(ROWS40):
            segs = [(a0, a1)] if a1 <= 24 and a0 < a1 else [(a0, 24), (0, (a1 % 24) if a1 > 24 or a0 > a1 else a1)]
            for s0, s1 in segs:
                ax.add_patch(plt.Rectangle((s0, i + .08), s1 - s0, .84, fill=False, ec=ACC, lw=.7, zorder=3))
        ax.set_yticks(np.arange(len(ROWS40)) + .5)
        if j == 0:
            nn = H[H.period == per].drop_duplicates('country').set_index('country').datasets
            ax.set_yticklabels([f'{r_[0].replace("United Kingdom", "UK").replace("United States", "USA")}' for r_ in ROWS40], fontsize=6.5)
        else:
            ax.set_yticklabels([])
        for i, (name, *_ , lab) in enumerate(ROWS40):
            nn = int(H[(H.period == per) & (H.country == name)].datasets.iloc[0])
            ax.text(24.25, i + .5, f'{nn:,}', fontsize=6, va='center', ha='left', color=MUTED)
        ax.tick_params(axis='y', length=0)
        ax.set_xticks([0, 6, 12, 18, 24]); ax.set_xlabel('Hour of the submission stamp (UTC)')
        ax.axvspan(8, 17, ymin=1.0, ymax=1.015, color=MUTED, lw=0, clip_on=False)
        ax.text(24.25, -.15, 'n', fontsize=6, va='bottom', ha='left', color=MUTED)
        ax.set_title(f'{per}', fontsize=7, fontweight='bold', loc='left', color=INK, pad=6)
        ax.text(12.5, -.3, '08–17 UTC', fontsize=6, ha='center', va='bottom', color=MUTED)
        panel(ax, 'a' if j == 0 else 'b', x=-.2 if j == 0 else -.04, y=1.03)
    # colour key + box key
    kx = fig.add_axes([.19, .425, .12, .012])
    kx.imshow(np.linspace(0, vmax, 100)[None, :], aspect='auto', cmap=cmap, extent=(0, vmax, 0, 1))
    kx.set_yticks([]); kx.set_xticks([0, .1, .2]); kx.set_xticklabels(['0', '10', '≥20%'], fontsize=6)
    kx.tick_params(length=1.5, pad=1)
    kx.text(-.01, .5, 'share of the row', transform=kx.transAxes, fontsize=6, ha='right', va='center', color=MUTED)
    fig.patches.append(plt.Rectangle((.38, .425), .018, .012, transform=fig.transFigure, fill=False, ec=ACC, lw=.7))
    fig.text(.402, .431, 'depositor-country local 09:00–17:00 (standard time; envelope of zones for Australia and USA)', fontsize=6, va='center', color=MUTED)
    # c: switch by region
    cx = fig.add_axes([.105, .075, .25, .255])
    rc = {'Europe': ACC2, 'Americas': ACC, 'Asia–Oceania': ACC3}
    for rg, cc in rc.items():
        z = Bt[Bt.region == rg]
        cx.fill_between(z.year, z.lo * 100, z.hi * 100, color=cc, alpha=.18, lw=0)
        cx.plot(z.year, z.share_08_17_utc * 100, color=cc, lw=1.5)
        cx.text(2025.3, z.share_08_17_utc.iloc[-1] * 100, rg, color=cc, fontsize=6, va='center', fontweight='bold')
    cx.set_ylim(0, 102); yfmt(cx); cx.set_xlim(2013.6, 2025.2); cx.set_xticks([2014, 2017, 2020, 2023, 2025])
    cx.set_ylabel('Stamped 08:00–16:59 UTC')
    cx.set_title('Share in the 08–17 UTC window, by region', fontsize=7, loc='left', color=INK)
    panel(cx, 'c', x=-.2)
    # d: weekday under three clocks
    dx = fig.add_axes([.555, .075, .235, .255])
    ser = [('submission stamp 2014–2019 (UTC)', MUTED, 'o', 'stamp 2014–19'), ('submission stamp 2023–2025 (local)', ACC2, 'o', 'stamp 2023–25, local'),
           ('acquisition date in raw-file names', INK, 's', 'raw-file dates')]
    xs = np.arange(7)
    for k, (col_, cc, mk, lab) in enumerate(ser):
        v = Cw[col_] * 100
        dx.plot(xs, v, color=cc, lw=1.1, marker=mk, ms=3, mec='white', mew=.3, zorder=3)
        dx.text(-.2, 10.5 - 3.3 * k, f'{lab}: {wk(Cw[col_]):.0%}', fontsize=6, color=cc, va='center', ha='left',
                fontweight='bold' if k == 1 else 'normal')
    dx.set_xticks(xs); dx.set_xticklabels([w_[:3] for w_ in WD]); dx.set_ylim(0, 26); yfmt(dx)
    dx.set_xlim(-.4, 6.4)
    dx.set_ylabel('Share of datasets')
    dx.set_title('Weekday under three clocks', fontsize=7, loc='left', color=INK)
    dx.text(-.2, 13.8, 'Saturday + Sunday share', fontsize=6, color=MUTED, va='center')
    dx.text(0, -.17, f'Local: principal zone of single-zone countries ({len(l_):,} datasets).\nRaw-file dates: each dated day once per dataset.',
            transform=dx.transAxes, fontsize=6, color=MUTED, va='top')
    panel(dx, 'd', x=-.2)
    # e: Israel (Sun-Thu working week)
    ex = fig.add_axes([.88, .075, .105, .255])
    for k, (col_, cc) in enumerate(zip(Iw.columns, [MUTED, ACC2, INK])):
        ex.barh(np.arange(7) + (k - 1) * .27, Iw[col_] * 100, height=.26, color=cc)
    ex.set_yticks(np.arange(7)); ex.set_yticklabels([w_[:3] for w_ in WD]); ex.invert_yaxis()
    ex.set_xticks([0, 20, 40]); ex.set_xticklabels(['0', '20', '40%'])
    ex.set_title('Israel (Sun–Thu week)', fontsize=7, loc='left', color=INK)
    ex.text(0, -.17, f'colours as d;\nn = {len(il_e)}, {len(il_l)}, {int(acq.loc["Israel", "datasets"])}', transform=ex.transAxes, fontsize=6, color=MUTED, va='top')
    panel(ex, 'e', x=-.32)
    save(fig, '40-weekly-rhythm')
    print('40', n_all, n_time, sh_early, sh_2025, Cw.apply(wk).to_dict(), Iw.to_dict(), (dec_sh, dec_lo, dec_hi), len(l_))


# ================================================================== 41 paper and release ====================
def sl(v):
    v = np.asarray(v, float); return np.sign(v) * np.log10(1 + np.abs(v))


def fig41():
    c = pd.read_csv(A / 'H11-clocks/Dataset-clocks.csv.gz', usecols=['accession', 'release_date', 'release_year', 'grp', 'paper_to_rel'])
    c['release_date'] = pd.to_datetime(c.release_date)
    f = F[['accession', 'lag_days']]
    c = c.merge(f, on='accession', validate='one_to_one')
    cnt = pd.to_datetime(F.release_date).value_counts()          # bulk = day with >= 50 releases in the full cohort
    c['bulk'] = c.release_date.map(cnt) >= 50
    b1 = pd.read_csv(O / 'Attacks/release-timing/B1-bulk-days.csv', parse_dates=['release_date'])
    assert set(c.release_date[c.bulk]) == set(b1.release_date)
    pl = c[c.lag_days >= 0]; n_pl = len(pl)
    p = pl[pl.paper_to_rel.notna()].copy(); n_p = len(p)
    assert n_p == 27001
    p['cat'] = np.select([p.paper_to_rel < -30, p.paper_to_rel <= 30, ~p.bulk], ['>30 d before', 'within ±30 d', '>30 d after, other day'],
                         '>30 d after, bulk day')
    # a: distributions (bulk days vs other days, same release years 2018-2025)
    q = p[p.release_year.between(2018, 2025)]
    S_ = {}
    for k, z in [('bulk', q[q.bulk]), ('other', q[~q.bulk]), ('other_all_years', p[~p.bulk])]:
        m, lo, hi = gb_median(z.paper_to_rel, z.grp); a_, alo, ahi = gb_share(z.paper_to_rel > 30, z.grp)
        w_, wlo, whi = gb_share(z.paper_to_rel.abs() <= 30, z.grp)
        S_[k] = dict(n=len(z), groups=z.grp.nunique(), median=m, lo=lo, hi=hi, after30=a_, after_lo=alo, after_hi=ahi, within30=w_, within_lo=wlo, within_hi=whi)
    S_ = pd.DataFrame(S_).T
    assert S_.loc['bulk', 'median'] == 95 and S_.loc['other_all_years', 'median'] == 11
    # integer-day bins, about 0.1 wide on the log axis but never narrower than one day; bar height = share per 0.1 log unit
    pos = np.unique(np.round(10 ** np.arange(0, 4.31, .1)).astype(int)) + .5
    edges_d = np.concatenate([-pos[::-1], [-.5], pos])
    edges = sl(edges_d)
    hist = []
    for k, z in [('bulk', q[q.bulk]), ('other', q[~q.bulk])]:
        h, _ = np.histogram(z.paper_to_rel.clip(edges_d[0] + .1, edges_d[-1] - .1), bins=edges_d)
        hist.append(pd.DataFrame({'set': k, 'bin_lo_days': edges_d[:-1], 'bin_hi_days': edges_d[1:], 'datasets': h, 'share': h / h.sum(),
                                  'share_per_0.1_log': h / h.sum() / (np.diff(edges) / .1)}))
    hist = pd.concat(hist)
    per = {}
    for lab, (y0, y1) in {'2013–2020': (2013, 2020), '2021–2025': (2021, 2025)}.items():
        z = p[p.release_year.between(y0, y1)]
        per[lab] = {'after': (z.paper_to_rel > 30).mean(), 'after_other': ((z.paper_to_rel > 30) & ~z.bulk).mean(),
                    'after_bulk': ((z.paper_to_rel > 30) & z.bulk).mean(), 'n': len(z)}
        per[lab]['after_ci'] = gb_share(z.paper_to_rel > 30, z.grp)[1:]
    pe, pl_ = per['2013–2020'], per['2021–2025']
    # b: composition by release year
    Y = pd.crosstab(p.release_year, p.cat, normalize='index')[['>30 d before', 'within ±30 d', '>30 d after, other day', '>30 d after, bulk day']]
    nY = p.release_year.value_counts().sort_index()
    allY = pl.release_year.value_counts().sort_index()
    cov = (nY / allY)
    wy = []
    for yr, z in p.groupby('release_year'):
        sh, lo, hi = gb_share(z.paper_to_rel.abs() <= 30, z.grp); wy.append((yr, sh, lo, hi))
    wy = pd.DataFrame(wy, columns=['release_year', 'within30', 'lo', 'hi']).set_index('release_year')
    out = Y.join(wy).assign(datasets_with_paper=nY, all_datasets=allY, share_with_dated_paper=cov)
    out.round(4).to_csv(D / '41-paper-and-release.csv')
    hist.round(5).to_csv(D / '41-paper-and-release-histogram.csv', index=False)
    S_.round(4).to_csv(D / '41-paper-and-release-summary.csv')

    fig = plt.figure(figsize=(W, 104 * MM))
    sb, so = S_.loc['bulk'], S_.loc['other']
    headline(fig, f'Datasets public >30 days after their paper rose from {pe["after"]:.0%} to {pl_["after"]:.0%} after 2020, '
                  f'all on bulk days',
             f'Public release minus publication date of the linked paper (OpenAlex): {n_p:,} of {n_pl:,} datasets have a dated paper. Bulk day: a calendar '
             f'day with ≥50 PRIDE\nreleases (62 days, all since Oct 2018). Release years 2013–20 vs 2021–25: >30 d after on other days '
             f'{pe["after_other"]:.0%} vs {pl_["after_other"]:.0%}; on bulk days {pe["after_bulk"]:.0%} vs {pl_["after_bulk"]:.0%}.\n'
             'Intervals: 95% over 2,000 resamples of contributor groups. Descriptive only: who set a release date is not recorded.')
    ax = fig.add_axes([.07, .12, .5, .62])
    ax.axvspan(sl(-30), sl(30), color=HAIR, alpha=.5, lw=0)
    for k, cc, sgn in [('bulk', ACC, 1), ('other', ACC2, -1)]:
        h = hist[hist.set == k]['share_per_0.1_log'].to_numpy() * 100
        ax.bar(edges[:-1], sgn * h, width=np.diff(edges), align='edge', color=cc, lw=.2, ec='white', alpha=.95)
    ax.axhline(0, color=INK, lw=.6)
    for k, cc, sgn in [('bulk', ACC, 1), ('other', ACC2, -1)]:
        r_ = S_.loc[k]; xm = sl(r_['median'])
        ax.plot([xm, xm], [0, sgn * 9.6], color=cc, lw=.8, ls=(0, (2, 1.5)))
        ax.text(xm + .05, sgn * 9.4, f'median {r_["median"]:+.0f} d ({r_["lo"]:.0f}–{r_["hi"]:.0f})', color=cc, fontsize=6.5,
                va='top' if sgn > 0 else 'bottom', fontweight='bold')
        lab = 'bulk days' if k == 'bulk' else 'other days'
        ax.text(-3.0, sgn * 8.2, f'{lab}, 2018–2025\n{int(r_["n"]):,} datasets, {int(r_["groups"]):,} groups', color=cc, fontsize=6.5,
                va='center', fontweight='bold', linespacing=1.05)
        ax.text(4.3, sgn * 6.2, f'{r_["after30"]:.0%} more than\n30 d after paper\n({r_["after_lo"]:.0%}–{r_["after_hi"]:.0%})', color=cc,
                fontsize=6.5, va='center', ha='right', linespacing=1.05)
    ax.text(0, 10.2, '±30 d', ha='center', va='bottom', fontsize=6, color=MUTED)
    for k, cc, sgn in [('bulk', ACC, 1), ('other', ACC2, -1)]:
        z = q[q.bulk] if k == 'bulk' else q[~q.bulk]
        d01 = z.paper_to_rel.between(0, 1).mean()
        ax.annotate(f'day 0–1: {d01:.1%}', xy=(sl(.5), sgn * (1.6 if k == 'other' else .3)), xytext=(sl(-12), sgn * (3.6 if k == 'other' else 2.2)),
                    fontsize=6, color=cc, ha='center', va='center', arrowprops=dict(arrowstyle='-', color=cc, lw=.5))
    ticks = [-1000, -100, -30, 0, 30, 100, 365, 1000, 3650]
    ax.set_xticks(sl(ticks)); ax.set_xticklabels(['−1000', '−100', '−30', '0', '30', '100', '1 y', '1000', '10 y'])
    ax.set_xlim(-3.1, 4.4); ax.set_ylim(-11, 11)
    ax.set_yticks([-10, -5, 0, 5, 10]); ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{abs(v):.0f}%'))
    ax.set_ylabel('Share of datasets per 0.1 log₁₀ unit')
    ax.set_xlabel('Days from paper publication to public release (symmetric log scale)')
    ax.text(-3.0, -10.6, '← public before the paper', fontsize=6, color=MUTED, va='bottom')
    ax.text(4.3, -10.6, 'public after the paper →', fontsize=6, color=MUTED, va='bottom', ha='right')
    panel(ax, 'a', x=-.09)
    # b: composition by release year
    bx = fig.add_axes([.665, .12, .25, .62])
    cols = {'>30 d before': '#d6d4cd', 'within ±30 d': '#8d8a83', '>30 d after, other day': ACC2, '>30 d after, bulk day': ACC}
    bot = np.zeros(len(Y))
    for k, cc in cols.items():
        v = Y[k].to_numpy() * 100
        bx.bar(Y.index, v, bottom=bot, color=cc, width=.78, lw=0)
        bot += v
    for yr in Y.index:
        bx.text(yr, 101, f'{cov[yr]:.0%}', ha='center', va='bottom', fontsize=6, color=MUTED, rotation=90)
    bx.text(2012.3, 101, 'with\npaper', ha='right', va='bottom', fontsize=6, color=MUTED, linespacing=1)
    bx.set_ylim(0, 100); yfmt(bx); bx.set_xlim(2012.4, 2025.6)
    bx.set_xticks([2013, 2016, 2019, 2022, 2025])
    bx.set_ylabel('Datasets with a dated paper, by release year')
    lab_y = {}
    last = Y.loc[2025] * 100; cum = last.cumsum() - last / 2
    for k, cc in cols.items():
        lab_y[k] = cum[k]
    lab_y = repel(lab_y, 9)
    for k, cc in cols.items():
        txt = {'>30 d before': 'public >30 d\nbefore paper', 'within ±30 d': 'within ±30 d', '>30 d after, other day': '>30 d after,\nother day',
               '>30 d after, bulk day': '>30 d after,\nbulk day'}[k]
        bx.text(2025.75, lab_y[k], f'{txt} {last[k]:.0f}%', fontsize=6, color=cc if k != '>30 d before' else MUTED, va='center',
                fontweight='bold', linespacing=1)
    panel(bx, 'b', x=-.2)
    save(fig, '41-paper-and-release')
    print('41', per, S_.round(3).to_string()); print(out.round(3).to_string())


if __name__ == '__main__':
    todo = sys.argv[1:] or ['38', '39', '40', '41']
    for k in todo:
        globals()['fig' + k]()
