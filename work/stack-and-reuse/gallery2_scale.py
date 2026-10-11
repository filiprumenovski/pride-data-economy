"""Gallery 2 — scale and instruments: figures 14 (fleet history), 15 (run ridgelines), 16 (depth vs breadth).

Descriptive only. Inputs are frozen Arsenal outputs (H1/H3) and the runs attack (Attacks/runs). A deposited
instrument label is not an ownership record and a model's first PRIDE year is not its purchase or launch year.
Each figure's plotted table is written to outputs/Stack-and-reuse/Gallery/data/<NN>-<slug>.csv.
Run from work/stack-and-reuse:  python gallery2_scale.py
"""
import numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.stats import gaussian_kde
from matplotlib.ticker import FuncFormatter, NullFormatter, FixedLocator
import matplotlib.colors as mc
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3, VEND
from common import SEED

H = A / 'H1-H3-scale-instruments'
RUNS = O / 'Attacks/runs'
D = G / 'data'; D.mkdir(exist_ok=True)
VEND_ALL = {**VEND, 'Applied Biosystems': '#8a7f73'}


def shade(hex_, t):
    """t in [0,1]: 0 = light tint of the vendor colour, 1 = deep shade."""
    r, g, b = mc.to_rgb(hex_)
    light = np.array([r, g, b]) * .45 + .55
    dark = np.array([r, g, b]) * .62
    return tuple(light + (dark - light) * t)


def wmedian(x, w):
    o = np.argsort(x); x, w = x[o], w[o]; c = np.cumsum(w)
    return x[np.searchsorted(c, c[-1] / 2)]


def gboot_median(df, col, reps=1000, seed=SEED):
    """95% percentile interval for a median, resampling whole pi_group clusters (accession when missing)."""
    g = df.pi_group.fillna(df.accession).to_numpy()
    u, inv = np.unique(g, return_inverse=True)
    x = df[col].to_numpy(float); rng = np.random.default_rng(seed); out = []
    for _ in range(reps):
        cnt = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u))
        w = cnt[inv]; k = w > 0
        out.append(wmedian(x[k], w[k]))
    return np.quantile(out, .025), np.quantile(out, .975)


# ----------------------------------------------------------------------------------------------------------------
# 14. Fleet history: every instrument model with >= 50 public datasets, as a ridgeline ordered by first appearance
# ----------------------------------------------------------------------------------------------------------------
def fig_fleet():
    m = pd.read_csv(H / 'H3-model-by-year.csv')
    div = pd.read_csv(H / 'H3-fleet-diversity-by-year.csv').set_index('release_year')
    tot = m.groupby(['model', 'vendor']).agg(datasets=('datasets', 'sum')).reset_index()
    keep = tot[tot.datasets >= 50].copy()
    first3 = m[m.datasets >= 3].groupby('model').release_year.min()      # same convention as gallery.m_first
    sh = m.pivot_table(index='model', columns='release_year', values='credit_share_labelled', aggfunc='sum').reindex(columns=YEARS).fillna(0) * 100
    n = m.pivot_table(index='model', columns='release_year', values='datasets', aggfunc='sum').reindex(columns=YEARS).fillna(0)
    keep['first'] = keep.model.map(first3)
    keep['peak_share'] = keep.model.map(sh.max(axis=1)); keep['peak_year'] = keep.model.map(sh.idxmax(axis=1))
    keep['n2025'] = keep.model.map(n[2025]).astype(int); keep['share2025'] = keep.model.map(sh[2025])
    keep['years_present_after_first'] = [int((n.loc[k, [y for y in YEARS if y >= f]] > 0).sum()) for k, f in zip(keep.model, keep['first'])]
    keep['years_possible'] = 2026 - keep['first']
    keep = keep.sort_values(['first', 'peak_year', 'datasets'], ascending=[True, True, False]).reset_index(drop=True)
    K = len(keep); alive = int((keep.n2025 > 0).sum()); never_gap = int((keep.years_present_after_first == keep.years_possible).sum())
    # vintage: share of 2025 model-labelled credit going to models (all models, not only these 37) first seen by 2015
    f_all = m[m.datasets >= 3].groupby('model').release_year.min()
    c25 = m[m.release_year == 2025].set_index('model').credit
    old25 = 100 * c25[c25.index.map(lambda x: f_all.get(x, 9999) <= 2015)].sum() / m[m.release_year == 2025].model_labelled_datasets.iloc[0]
    e13, e25 = div.effective_models_inv_simpson[2013], div.effective_models_inv_simpson[2025]
    nlab = int(m.groupby('release_year').model_labelled_datasets.first().loc[YEARS].sum())

    # colours: vendor family, shade within vendor by first year (light = earlier, deep = later)
    col = {}
    for v, g in keep.groupby('vendor'):
        fy = g['first']; lo, hi = 2013, 2024
        for j, (k, f) in enumerate(zip(g.model, fy)):
            col[k] = shade(VEND_ALL[v], np.clip((f - lo) / (hi - lo) * .8 + .12 + (.07 if j % 2 else -.03), 0, 1))

    fig = plt.figure(figsize=(W, 150 * MM))
    headline(fig, f'Instrument models are added but none drop out: all {K} common models still appear in 2025 deposits',
             f'Rows: every instrument model named in ≥50 public PRIDE datasets released 2013–2025 ({nlab:,} model-labelled datasets; fractional credit when a dataset\n'
             f'names several models), ordered by first year with ≥3 datasets. Ridge height = √(share of that year’s model-labelled datasets), common scale. Thin\n'
             'line = years with ≥1 dataset. Rows starting in 2013 are left-censored. A deposited label is not an ownership record; first year is not purchase.')
    if not JOURNAL: fig.texts[0].set_y(.992); fig.texts[1].set_y(.962); fig.texts[1].set_fontsize(6.5)
    x0, x1 = .035, .72; ytop = .875
    ax = fig.add_axes([x0, .045, x1 - x0, ytop - .045])
    sp = 1.0; hmax = 1.6 * sp; smax = 40.0                     # 40% share -> 1.6 rows tall
    hs = lambda v: np.sqrt(np.clip(v, 0, None) / smax) * hmax
    xs = np.linspace(2013, 2025, 400)
    disp = {'6520A Quadrupole Time-of-Flight LC/MS': '6520A Q-TOF LC/MS', '4800 Proteomics Analyzer': '4800 Proteomics Analyzer'}
    for i, r in keep.iterrows():
        k = r.model; base = (K - 1 - i) * sp
        y = sh.loc[k].to_numpy()
        pres = [yy for yy in YEARS if n.loc[k, yy] > 0]
        f = hs(PchipInterpolator(YEARS, y)(xs)); f[xs < min(pres)] = 0      # no ridge before the first observed dataset
        z = 200 - i * 2                                       # upper rows in front: thin ridges stay visible over big ones
        ax.fill_between(xs, base, base + f, where=f > 1e-4, color=col[k], lw=0, zorder=z)
        ax.plot(xs, np.where(f > 1e-4, base + f, np.nan), color=mc.to_hex(shade(VEND_ALL[r.vendor], 1.0)), lw=.45, zorder=z + .5)
        ax.plot([min(pres), 2025], [base, base], color=col[k], lw=.8, zorder=z + .6, solid_capstyle='butt')
        ax.plot(2025, base, 'o', ms=2, color=col[k], zorder=z + .7, mec='none')
        lx = min(pres) - .2 if min(pres) > 2013 else 2012.8
        ax.text(lx, base + .12, disp.get(k, k), ha='right', va='bottom', fontsize=6, color=mc.to_hex(shade(VEND_ALL[r.vendor], .95)), zorder=300)
    ytopd = (K - 1) * sp + 2.9
    ax.set_xlim(2008.9, 2025.45); ax.set_ylim(-.6, ytopd)
    ax.set_yticks([])
    ax.set_xticks(range(2013, 2026, 2))
    # vendor key (direct, coloured words)
    xk = 2025.0
    for v in ['Applied Biosystems', 'Agilent', 'Waters', 'SCIEX', 'Bruker', 'Thermo Fisher']:
        t = ax.text(xk, ytopd - .25, v, ha='right', va='top', fontsize=6.3, fontweight='bold', color=mc.to_hex(shade(VEND_ALL[v], .8)))
        fig.canvas.draw(); bb = t.get_window_extent().transformed(ax.transData.inverted()); xk = bb.x0 - .35
    # scale key for ridge height
    n13 = int((keep['first'] == 2013).sum())
    # right columns: peak share bar, 2025 datasets
    axr = fig.add_axes([x1 + .035, .045, .11, ytop - .045], sharey=ax)
    axn = fig.add_axes([x1 + .16, .045, .1, ytop - .045], sharey=ax)
    for a in (axr, axn):
        a.tick_params(axis='y', left=False, labelleft=False)
    for i, r in keep.iterrows():
        base = (K - 1 - i) * sp + .3
        c = col[r.model]
        axr.barh(base, r.peak_share, height=.6, color=c, lw=0)
        axr.text(r.peak_share + .8, base, f'{r.peak_share:.1f}%' if r.peak_share < 10 else f'{r.peak_share:.0f}%', va='center', fontsize=6, color=INK)
        axn.barh(base, r.n2025, height=.6, color=c, lw=0)
        axn.text(r.n2025 + 25, base, f'{r.n2025:,}', va='center', fontsize=6, color=INK)
    axr.set_xlim(0, 52); axr.set_xticks([0, 20, 40]); axr.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%'))
    axn.set_xlim(0, 1350); axn.set_xticks([0, 500, 1000])
    axr.text(2, ytopd - .25, 'Peak share', fontsize=6.5, color=INK, va='top')
    axn.text(40, ytopd - .25, 'Datasets in 2025', fontsize=6.5, color=INK, va='top')
    # inset: effective number of models, placed in the empty lower-left triangle (data coordinates of ax)
    ai = ax.inset_axes([2010.2, 1.7, 5.9, 5.4], transform=ax.transData)
    e = div.loc[YEARS]
    ai.fill_between(YEARS, e.eff_ci_lo, e.eff_ci_hi, color=ACC2, alpha=.18, lw=0)
    ai.plot(YEARS, e.effective_models_inv_simpson, color=ACC2, lw=1.4)
    ai.text(2013, e.effective_models_inv_simpson[2013] - 1.2, f'{e13:.1f}', ha='left', va='top', fontsize=8, fontweight='bold', color=ACC2)
    ai.text(2024.5, e.effective_models_inv_simpson[2025] + .3, f'{e25:.1f}', ha='right', va='bottom', fontsize=8, fontweight='bold', color=ACC2)
    ai.set_ylim(0, 16.5); ai.set_xlim(2012.7, 2025.3); ai.set_xticks([2013, 2019, 2025]); ai.set_yticks([0, 5, 10, 15])
    ai.tick_params(labelsize=6); ai.set_zorder(400)
    ai.set_title('Effective number of models\n(inverse Simpson; 95% group CI)', fontsize=6.3, loc='left', color=INK, pad=2)
    ax.text(2009.0, 10.9, f'In 2025, {old25:.0f}% of model-labelled credit\nwent to models already in PRIDE by 2015', ha='left', va='top',
            fontsize=6.3, color=MUTED, linespacing=1.1)
    panel(ax, 'a', x=-.01); panel(axr, 'b', x=-.03); panel(axn, 'c', x=-.03); panel(ai, 'd', x=-.13, y=1.17)
    out = keep.assign(color=[mc.to_hex(col[k]) for k in keep.model])
    out = out.merge(sh.add_prefix('share_').reset_index(), on='model').merge(n.add_prefix('n_').reset_index(), on='model')
    out.to_csv(D / '14-fleet-history.csv', index=False)
    save(fig, '14-fleet-history')
    return dict(models=K, alive_2025=alive, never_gap=never_gap, eff=(e13, e25), old25=old25, n13=n13)


# ----------------------------------------------------------------------------------------------------------------
# 15. Runs per dataset: the whole distribution moved left while the >=100-run tail grew in absolute numbers
# ----------------------------------------------------------------------------------------------------------------
def kde_log(v, grid, bw=.085):
    """Gaussian KDE of log10(values) evaluated on grid (log10 units); integrates to 1."""
    x = np.log10(v); h, e = np.histogram(x, bins=np.arange(-.3, 4.3, .005))
    c = (e[:-1] + e[1:]) / 2
    k = np.exp(-.5 * ((grid[:, None] - c[None, :]) / bw) ** 2) @ h
    return k / (k.sum() * (grid[1] - grid[0]))


def fig_runs():
    d = pd.read_csv(H / 'H1-dataset-scale.csv.gz')
    pop = d[d.release_year.between(2013, 2025)]; nall = len(pop)
    d = pop[pop.runs_lb >= 1].copy()                                       # claim population: >=1 identified run or raw bundle
    s = pd.read_csv(H / 'H1-scale-by-year.csv').set_index('release_year')
    ci = pd.read_csv(H / 'H1-year-medians-with-CI.csv').set_index('release_year')
    c = pd.read_csv(RUNS / 'Contrasts-2013-16-vs-2023-25.csv').set_index('label')
    allc = c.loc['all datasets with runs (claim population)']
    grid = np.linspace(-.25, 4.05, 700)
    rows = []
    for y in YEARS:
        x = d[d.release_year == y].runs_lb
        assert len(x) == s.datasets_with_raw[y] and x.median() == s.median_runs_lb[y]
        rows.append(dict(release_year=y, datasets=len(x), median=x.median(), ci_lo=ci.ci_lo[y], ci_hi=ci.ci_hi[y],
                         ge100=int((x >= 100).sum()), share_ge100=100 * (x >= 100).mean(), ge100_ci_lo=100 * s.ge100_ci_lo[y],
                         ge100_ci_hi=100 * s.ge100_ci_hi[y], p10=x.quantile(.1), p90=x.quantile(.9)))
    t = pd.DataFrame(rows).set_index('release_year')
    dens = {y: kde_log(d[d.release_year == y].runs_lb.to_numpy(float), grid) for y in YEARS}
    peak = max(v.max() for v in dens.values())

    fig = plt.figure(figsize=(W, 118 * MM))
    headline(fig, f'Typical datasets hold half as many raw runs as in 2013–16, yet the number with ≥100 runs rose from {t.ge100[2013]} in 2013 to {t.ge100[2025]} in 2025',
             f'PRIDE datasets released 2013–2025 with ≥1 identified raw run or raw bundle (n = {len(d):,} of {nall:,}). Runs = vendor raw acquisitions, a lower bound when\n'
             f'runs are bundled; not samples (TMT puts up to 18 samples in one run). Ridges: distribution per release year (equal area, log axis); white tick = median.\n'
             f'Pooled medians {allc.median_early:.0f} (2013–16, n = {allc.n_early:,}) and {allc.median_late:.0f} (2023–25, n = {allc.n_late:,}); intervals are 95% group bootstraps over hashed PI groups.')
    if not JOURNAL: fig.texts[0].set_y(.99); fig.texts[1].set_y(.948); fig.texts[1].set_fontsize(6.5)
    top = .82; bot = .1
    ax = fig.add_axes([.075, bot, .44, top - bot])
    sp = 1.0; hmax = 2.0 * sp; ny = len(YEARS)
    body, tail = '#9fb6c9', ACC
    for i, y in enumerate(YEARS):
        base = (ny - 1 - i) * sp; f = dens[y] / peak * hmax; xg = 10 ** grid; z = 10 + i * 3
        ax.fill_between(xg, base, base + f, color='white', lw=0, zorder=z)
        ax.fill_between(xg, base, base + f, where=xg <= 100, color=body, lw=0, zorder=z + .1)
        ax.fill_between(xg, base, base + f, where=xg >= 100, color=tail, lw=0, zorder=z + .1, alpha=.9)
        ax.plot(xg, base + f, color='#3d5a73', lw=.55, zorder=z + .2)
        ax.plot([xg[0], xg[-1]], [base, base], color='#3d5a73', lw=.4, zorder=z + .2)
        md = t['median'][y]; fm = np.interp(np.log10(md), grid, f)
        ax.plot([md, md], [base, base + fm], color='white', lw=1.2, zorder=z + .3)
        ax.text(.5, base + .12, str(y), ha='right', va='bottom', fontsize=6.5, color=INK, zorder=z + .4)
        ax.text(.5, base - .02, f'n = {t.datasets[y]:,}', ha='right', va='top', fontsize=6, color=MUTED, zorder=z + .4)
    yb = [(ny - 1 - i) * sp for i in range(ny)]
    ytopd = (ny - 1) * sp + hmax + .5
    # pooled period medians as reference lines through all ridges
    for v, lab in [(allc.median_early, '2013–16'), (allc.median_late, '2023–25')]:
        ax.plot([v, v], [-.35, ytopd - .55], color=INK, lw=.5, ls=(0, (3, 2)), zorder=200, alpha=.75)
        ax.text(v, ytopd - .5, f'{v:.0f}', ha='center', va='bottom', fontsize=7.5, fontweight='bold', color=INK)
    ax.text(np.sqrt(allc.median_early * allc.median_late), ytopd + .35, 'pooled median\n2023–25 ← 2013–16', ha='center', va='bottom', fontsize=6, color=MUTED, linespacing=1)
    ax.set_xscale('log'); ax.set_xlim(.11, 10000); ax.set_ylim(-.75, ytopd + 1.15)
    ax.xaxis.set_major_locator(FixedLocator([1, 3, 10, 30, 100, 300, 1000, 3000]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}')); ax.xaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_minor_locator(FixedLocator([m_ * 10 ** e_ for e_ in range(0, 4) for m_ in range(2, 10)]))
    ax.set_yticks([])
    ax.set_xlabel('Raw runs per dataset (log scale)')
    ax.axvline(100, color=tail, lw=.6, ls=(0, (2, 2)), zorder=1)
    ax.text(115, ytopd - .5, '≥100 runs', color=tail, fontsize=6.5, fontweight='bold', va='bottom')
    # count column: absolute number of >=100-run datasets, with share
    axc = fig.add_axes([.535, bot, .12, top - bot], sharey=ax)
    axc.tick_params(axis='y', left=False, labelleft=False)
    for i, y in enumerate(YEARS):
        b = yb[i] + .3
        axc.barh(b, t.ge100[y], height=.55, color=tail, lw=0, alpha=.9)
        axc.text(t.ge100[y] + 12, b, f'{t.ge100[y]}', va='center', fontsize=6.3, color=INK, fontweight='bold' if y in (2013, 2025) else 'normal')
        axc.text(t.ge100[y] + 12, b - .42, f'{t.share_ge100[y]:.1f}%', va='center', fontsize=6, color=MUTED)
    axc.set_xlim(0, 760); axc.set_xticks([0, 250, 500])
    axc.text(18, ytopd - .5, '≥100-run datasets\n(count; % of year)', fontsize=6.5, va='bottom', color=INK)
    # dumbbells: pooled 2013-16 vs 2023-25 medians in subsets (frozen runs attack)
    sel = [('all datasets with runs (claim population)', 'All datasets'),
           ('like-for-like: Thermo, DIA & TMT not detected, human, publication linked', 'Like-for-like¹'),
           ('grain: publication component', 'Summed per publication'),
           ('fraction families collapsed (sample proxy)', 'Fractions collapsed²'),
           ('vendor4=Thermo Fisher', 'Thermo Fisher'),
           ('dia_s=DIA detected', 'DIA detected'),
           ('tmt_s=TMT detected', 'TMT detected'),
           ('vendor4=SCIEX', 'SCIEX'),
           ('vendor4=Bruker', 'Bruker³')]
    dd = c.loc[[k for k, _ in sel]].copy(); dd['name'] = [n_ for _, n_ in sel]
    dd['falls'] = dd.median_ratio_hi < 1
    axd = fig.add_axes([.8, bot, .17, top - bot])
    yy = np.arange(len(dd))[::-1]
    for yv, (_, r) in zip(yy, dd.iterrows()):
        colr = ACC2 if r.falls else MUTED
        axd.plot([r.median_early, r.median_late], [yv, yv], color=colr, lw=1.2, alpha=.5, zorder=2)
        axd.plot([r.median_early_lo, r.median_early_hi], [yv + .17, yv + .17], color='#b9b9b9', lw=.9, zorder=2)
        axd.plot([r.median_late_lo, r.median_late_hi], [yv - .17, yv - .17], color=colr, lw=.9, zorder=2, alpha=.7)
        axd.scatter(r.median_early, yv, s=16, facecolor='white', edgecolor=colr, lw=.9, zorder=3)
        axd.scatter(r.median_late, yv, s=16, color=colr, zorder=3, lw=0)
        axd.text(4.3, yv, r['name'], ha='right', va='center', fontsize=6.2, color=INK if r.falls else MUTED)
        axd.text(85, yv, f'×{r.median_ratio:.2f}', ha='left', va='center', fontsize=6, color=colr)
    axd.set_xscale('log'); axd.set_xlim(5, 70); axd.set_ylim(-.7, len(dd) + .75)
    axd.xaxis.set_major_locator(FixedLocator([5, 10, 20, 40])); axd.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    axd.xaxis.set_minor_formatter(NullFormatter()); axd.xaxis.set_minor_locator(FixedLocator([]))
    axd.set_yticks([])
    axd.set_xlabel('Median runs per dataset')
    axd.text(.04, len(dd) + .6, 'Median runs, 2013–16 (○)\nvs 2023–25 (●); 95% CIs', transform=axd.get_yaxis_transform(), fontsize=6.5, ha='left', va='top', color=INK)
    footnote(fig, .685, bot - .085, '¹ Thermo, no DIA/TMT, human, publication-linked.\n² Fraction-numbered runs counted once (sample proxy).\n'
             f'³ Bruker 2013–16 n = {int(dd.loc["vendor4=Bruker", "n_early"])}. Grey: CI of the ratio includes 1.', fontsize=6, color=MUTED, va='top', ha='left')
    panel(ax, 'a', x=-.02); panel(axc, 'b', x=-.06); panel(axd, 'c', x=-.62)
    t.reset_index().to_csv(D / '15-run-ridgelines.csv', index=False)
    dd.reset_index().to_csv(D / '15-run-ridgelines-subsets.csv', index=False)
    pd.DataFrame({'log10_runs': grid, **{str(y): dens[y] for y in YEARS}}).to_csv(D / '15-run-ridgelines-density.csv', index=False)
    save(fig, '15-run-ridgelines')
    return t[['datasets', 'median', 'ge100', 'share_ge100']].loc[[2013, 2025]].to_dict(), dict(allc[['median_early', 'median_late', 'n_early', 'n_late']])


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['14', '15', '16']
    if '14' in which:
        print('14', fig_fleet())
    if '15' in which:
        print('15', fig_runs())
