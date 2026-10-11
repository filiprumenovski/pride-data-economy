"""Gallery figure 32: sequence databases (FASTA) at deposit (descriptive).

32-database-age  a  Age of the main deposited FASTA at submission (submission date minus a high-confidence date in the
                    file name; largest dated FASTA per dataset), by submission period: density outlines + medians with
                    95% group-bootstrap intervals.
                 b  Like-for-like naming check, 2019-2025 submissions: median age by named software family, separately for
                    YYYY-MM-DD names (FragPipe's auto-download style) and all other date formats.
                 c  Share of datasets depositing any FASTA file (top level or inside a listed zip), by release year.
                 d  Source named in FASTA file names, % of datasets with FASTA, by release period (flags overlap).

Sources: Arsenal/H4-H7-databases-storage (H4-dataset-database-date.csv.gz, H4-fasta-rate-by-year.csv,
H4-sources-by-period-dataset-pct.csv), Deposit-features + software flags (common.features()).
A name date can be a download, a UniProt release or a tool build date; dates > 31 days after submission were already
excluded upstream (H4.6). "No source named" means the file names do not name one. Descriptive only.
Uncertainty: 95% percentile intervals resampling hashed pi_group clusters (no group = own cluster), 2,000 resamples,
seed 20261001. Plotted tables -> Gallery/data/32-database-age*.csv.
"""
import warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd
import matplotlib.patheffects as pe
from matplotlib.ticker import FuncFormatter
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3
from common import O, SEED, features
from gallery3_versions import headmm, group_boot_median, darken, TOOL_COL as _TC
TOOL_COL = {**_TC, 'Spectronaut': '#e8a33d', 'Mascot': '#7f8c96'}

D = G / 'data'; D.mkdir(exist_ok=True)
H4 = A / 'H4-H7-databases-storage'
HALO = [pe.withStroke(linewidth=2, foreground='white')]


def fig_database_age():
    d = features()[['accession', 'family', 'release_year', 'submission_date', 'pi_group']]
    b = pd.read_csv(H4 / 'H4-dataset-database-date.csv.gz').merge(d, on='accession', how='inner', validate='one_to_one')
    b['grp'] = b.pi_group.fillna('acc:' + b.accession)
    h = b[b.date_conf == 'high'].copy()
    h['sy'] = h.submission_date.dt.year
    n_all, g_all = len(h), h.grp.nunique()
    med_all = h.years.median(); ci_all = group_boot_median(h.years, h.grp)
    gt2 = 100 * (h.years > 2).mean(); gt1 = 100 * (h.years > 1).mean(); gt5 = 100 * (h.years > 5).mean()
    fasta_ds = pd.read_csv(H4 / 'H4-sources-by-period-dataset-pct.csv').iloc[-1, 1:].astype(float).sum()

    # a: by submission period
    PER = [('2013–18', 2012, 2018), ('2019–21', 2019, 2021), ('2022–25', 2022, 2025)]
    PCOL = ['#a9bccd', '#5f87aa', '#1d4f75']
    grid = np.linspace(0, 8, 321)
    arows, dens = [], {}
    for (lab, y0, y1) in PER:
        x = h[h.sy.between(y0, y1)]
        lo, hi = group_boot_median(x.years, x.grp)
        arows.append({'period': lab, 'datasets': len(x), 'groups': x.grp.nunique(), 'median_years': x.years.median(), 'ci_lo': lo, 'ci_hi': hi,
                      'pct_gt2y': 100 * (x.years > 2).mean(), 'pct_gt5y': 100 * (x.years > 5).mean()})
        v = x.years.clip(lower=0).to_numpy()
        bwk = .9 * min(v.std(), (np.quantile(v, .75) - np.quantile(v, .25)) / 1.34) * len(v) ** -.2
        # reflected Gaussian KDE (ages are >= 0)
        k = (np.exp(-.5 * ((grid[:, None] - v[None]) / bwk) ** 2) + np.exp(-.5 * ((grid[:, None] + v[None]) / bwk) ** 2)).sum(1)
        dens[lab] = k / (len(v) * bwk * np.sqrt(2 * np.pi))
    AR = pd.DataFrame(arows)

    # b: naming check by software, 2019-25
    fams = ['MaxQuant', 'Proteome Discoverer', 'Mascot', 'Spectronaut', 'DIA-NN', 'FragPipe']
    hb = h[h.sy.between(2019, 2025)].copy(); hb['ymd'] = hb.date_format.eq('YYYY-MM-DD')
    brows = []
    for f in fams:
        for fmt, sel in (('YYYY-MM-DD', True), ('other formats', False)):
            x = hb[(hb.family == f) & (hb.ymd == sel)]
            lo, hi = group_boot_median(x.years, x.grp) if len(x) >= 10 else (np.nan, np.nan)
            brows.append({'family': f, 'name_format': fmt, 'datasets': len(x), 'groups': x.grp.nunique(), 'median_years': x.years.median(),
                          'ci_lo': lo, 'ci_hi': hi})
        x = hb[hb.family == f]
        lo, hi = group_boot_median(x.years, x.grp)
        brows.append({'family': f, 'name_format': 'all', 'datasets': len(x), 'groups': x.grp.nunique(), 'median_years': x.years.median(),
                      'ci_lo': lo, 'ci_hi': hi, 'pct_ymd_names': 100 * x.ymd.mean()})
    BR = pd.DataFrame(brows)

    # c: FASTA deposit rate
    r = pd.read_csv(H4 / 'H4-fasta-rate-by-year.csv')
    rc = r[r.measure == 'fasta_anywhere'].set_index('release_year')
    # d: sources
    s = pd.read_csv(H4 / 'H4-sources-by-period-dataset-pct.csv').set_index('name flag (pct of datasets with FASTA)')
    nden = s.loc['datasets with FASTA (denominator, n)']
    SRC = [('only sources not named', 'No FASTA names a source', MUTED),
           ('UniProt (any mention)', 'UniProt (any)', ACC2),
           ('UniProt proteome ID (UP0...)', 'UniProt proteome ID', '#5a9e6f'),
           ('Swiss-Prot / reviewed (explicit)', 'Swiss-Prot', '#8d6cab'),
           ('Contaminants (cRAP / contaminant list)', 'Contaminant list', ACC3),
           ('Genome-, transcript- or sample-derived (custom, variant, ORF, translation, metagenome, assembly...)', 'Genome/sample-derived', ACC),
           ('Ensembl / GENCODE', 'Ensembl/GENCODE', '#3d3d3d')]

    pd.concat([AR.assign(panel='a'), BR.assign(panel='b')], ignore_index=True).round(3).to_csv(D / '32-database-age.csv', index=False)
    pd.DataFrame({'age_years': grid, **{k: v for k, v in dens.items()}}).round(5).to_csv(D / '32-database-age-density.csv', index=False)
    rc.reset_index().to_csv(D / '32-database-age-fasta-rate.csv', index=False)
    s.loc[[k for k, _, _ in SRC] + ['datasets with FASTA (denominator, n)']].to_csv(D / '32-database-age-sources.csv')

    H = 140
    fig = plt.figure(figsize=(W, H * MM))
    headmm(fig, f'Deposited sequence databases are typically {med_all:.1f} years old, and {gt2:.0f}% are older than two years',
           f'PRIDE datasets released 2013–2025. a, b: age of the main FASTA at submission, from a high-confidence date in its file name ({n_all:,} datasets, '
           f'{g_all:,}\ngroups; median {med_all:.2f} y, 95% CI {ci_all[0]:.2f}–{ci_all[1]:.2f}; >1 y {gt1:.0f}%, >5 y {gt5:.0f}%). '
           'A name date may be a download, release or tool build date, so ages are approximate.', H)

    # ---- a
    ax = fig.add_axes([.06, .555, .5, .3])
    for (lab, *_), c in zip(PER, PCOL):
        ax.plot(grid, dens[lab], color=c, lw=1.5 if lab != '2022–25' else 1.8)
    ymax = max(v.max() for v in dens.values())
    for i, ((lab, *_), c) in enumerate(zip(PER, PCOL)):
        rr = AR.iloc[i]; yy = ymax * (1.2 - .1 * i)
        ax.plot([rr.ci_lo, rr.ci_hi], [yy, yy], color=c, lw=2.2, solid_capstyle='butt')
        ax.plot(rr.median_years, yy, '|', ms=6, mew=1.3, color=darken(c, .8))
        ax.text(2.3, yy, f'{lab}: median {rr.median_years:.2f} y ({rr.ci_lo:.2f}–{rr.ci_hi:.2f}), >2 y {rr.pct_gt2y:.0f}%, n = {int(rr.datasets):,}',
                fontsize=6.3, va='center', color=darken(c, .85 if i == 2 else (.7 if i == 1 else .62)))
    for xv, lab in ((1, '1 y'), (2, '2 y'), (5, '5 y')):
        ax.axvline(xv, color=HAIR, lw=.6, zorder=0)
    ax.set_xlim(0, 8); ax.set_ylim(0, ymax * 1.32); ax.set_yticks([])
    ax.set_xlabel('Database age at submission (years)')
    ax.set_ylabel('Density', fontsize=6.5)
    ax.text(2.3, ymax * .91, 'bars: median and 95% group-bootstrap CI', fontsize=6, va='center', color=MUTED)
    ax.text(0, 1.04, 'Age distribution is stable across submission periods', transform=ax.transAxes, fontsize=7.5, fontweight='bold', color=INK)
    panel(ax, 'a', x=-.035, y=1.04)

    # ---- b
    bx = fig.add_axes([.70, .555, .245, .3])
    for j, f in enumerate(fams):
        yy = len(fams) - 1 - j
        q = BR[BR.family == f].set_index('name_format')
        c = TOOL_COL[f]
        a_, o_, al = q.loc['YYYY-MM-DD'], q.loc['other formats'], q.loc['all']
        for row, dy in ((a_, .16), (o_, -.16)):
            bx.plot([row.ci_lo, row.ci_hi], [yy + dy] * 2, color=c, lw=.8, zorder=2, clip_on=True)
            filled = row.name == 'other formats'
            bx.plot(row.median_years, yy + dy, 'o', ms=4.4, mfc=c if filled else 'white', mec=c, mew=1, zorder=3)
        bx.plot(al.median_years, yy, marker='|', ms=11, mew=1.1, color=INK, zorder=4)
        nm = f + (' †' if a_.groups < 5 else '')
        bx.text(-.06, yy, nm, ha='right', va='center', fontsize=6.3, color=darken(c, .75), transform=bx.get_yaxis_transform())
        bx.text(1.03, yy, f'{int(a_.datasets)} / {int(o_.datasets)}', ha='left', va='center', fontsize=6, color=MUTED, transform=bx.get_yaxis_transform())
    mas = BR[(BR.family == 'Mascot') & (BR.name_format == 'YYYY-MM-DD')].iloc[0]
    bx.text(1.03, len(fams) - .25, 'n', ha='left', va='center', fontsize=6, color=MUTED, transform=bx.get_yaxis_transform())
    bx.set_yticks([]); bx.set_ylim(-.6, len(fams) - .2); bx.set_xlim(0, 2.8)
    bx.set_xlabel('Median age at submission (years)'); bx.grid(axis='x', color=HAIR, lw=.4); bx.set_axisbelow(True)
    bx.text(-.5, 1.04, 'YYYY-MM-DD names read younger in every tool', transform=bx.transAxes, fontsize=7.5, fontweight='bold', color=INK)
    for xx, kw, tx in ((-.5, dict(mfc='white', mec=INK, mew=.9, ms=4), 'YYYY-MM-DD names'), (.1, dict(color=INK, ms=4), 'other dates'),
                       (.53, dict(marker='|', ms=7, mew=1.2, color=INK), 'all names')):
        bx.plot(xx, -.25, **({'marker': 'o'} | kw), transform=bx.transAxes, clip_on=False, ls='none')
        bx.text(xx + .04, -.25, tx, fontsize=6, va='center', transform=bx.transAxes, color=INK)
    bx.text(-.5, -.3, f'2019–25 submissions; lines: 95% CI; † YYYY-MM-DD names from {int(mas.groups)} groups', fontsize=6,
            va='top', transform=bx.transAxes, color=MUTED, linespacing=1.2)
    panel(bx, 'b', x=-.52, y=1.04)

    # ---- c
    cx = fig.add_axes([.06, .085, .31, .32])
    yrs = YEARS
    cx.fill_between(yrs, rc.ci_lo[yrs], rc.ci_hi[yrs], color=ACC2, alpha=.16, lw=0)
    cx.plot(yrs, rc.pct[yrs], color=ACC2, lw=1.7, marker='o', ms=2.4)
    for yv in (2013, 2025):
        cx.text(yv + (.35 if yv == 2013 else -.1), rc.pct[yv] + (-.3 if yv == 2013 else 2.5), f'{rc.pct[yv]:.1f}%', fontsize=7.5,
                fontweight='bold', color=ACC2, ha='left' if yv == 2013 else 'center', va='top' if yv == 2013 else 'baseline', path_effects=HALO)
    cx.set_ylim(0, 30); yfmt(cx); cx.set_xlim(2012.5, 2025.5); cx.set_xticks([2013, 2016, 2019, 2022, 2025])
    cx.text(0, 1.04, 'Datasets depositing any FASTA', transform=cx.transAxes, fontsize=7.5, fontweight='bold', color=INK)
    cx.text(.97, .06, f'Release year; top level or inside listed zips.\nn = {int(rc.datasets[2013]):,} (2013) to {int(rc.datasets[2025]):,} (2025); 95% CI',
            transform=cx.transAxes, fontsize=6, color=MUTED, va='bottom', ha='right', linespacing=1.2)
    panel(cx, 'c', x=-.1, y=1.04)

    # ---- d
    dx = fig.add_axes([.465, .085, .3, .32])
    P = ['2013-16', '2017-19', '2020-22', '2023-25']; xp = np.arange(4)
    ends = {}
    for key, lab, c in SRC:
        vals = s.loc[key, P].astype(float).to_numpy()
        dx.plot(xp, vals, color=c, lw=1.8 if key in ('only sources not named', 'UniProt proteome ID (UP0...)') else 1.1, marker='o', ms=2.4)
        ends[lab] = vals[-1]
    lab_y = repel(ends, 3.4)
    for key, lab, c in SRC:
        v0, v1 = s.loc[key, P[0]], s.loc[key, P[-1]]
        dx.annotate(f'{lab}  {v0:.0f}→{v1:.0f}%', xy=(3, ends[lab]), xytext=(3.2, lab_y[lab]), fontsize=6.2, va='center', color=darken(c, .8) if c != MUTED else '#55534f',
                    fontweight='bold' if key in ('only sources not named', 'UniProt proteome ID (UP0...)') else 'normal', annotation_clip=False,
                    arrowprops=dict(arrowstyle='-', color=HAIR, lw=.4, shrinkA=0, shrinkB=1))
    dx.set_xticks(xp); dx.set_xticklabels([f'{p.replace("-", "–")}\nn = {int(nden[p]):,}' for p in P], fontsize=6.2, linespacing=1.2)
    dx.set_ylim(0, 60); yfmt(dx); dx.set_xlim(-.2, 3.05)
    dx.text(0, 1.04, 'Source named in FASTA file names', transform=dx.transAxes, fontsize=7.5, fontweight='bold', color=INK)
    dx.set_ylabel('% of datasets with FASTA (flags overlap)', fontsize=6.5)
    panel(dx, 'd', x=-.1, y=1.04)
    save(fig, '32-database-age')
    return dict(n=n_all, groups=g_all, med=med_all, ci=ci_all, gt1=gt1, gt2=gt2, gt5=gt5, AR=AR, BR=BR, rc=rc, fasta_ds=fasta_ds)


if __name__ == '__main__':
    r = fig_database_age()
    print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k not in ('AR', 'BR', 'rc')})
    print(r['AR'].round(2).to_string()); print(r['BR'].round(2).to_string())
