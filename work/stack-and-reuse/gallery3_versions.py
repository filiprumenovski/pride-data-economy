"""Gallery figures 31 and 33: software release lines and software-release age in public deposits (descriptive).

31-version-streams      Which release line of MaxQuant, Proteome Discoverer, DIA-NN and FragPipe is recorded in
                        deposits submitted each year. Stacked shares by submission year (= deposit year);
                        each band starts at the line's first appearance in PRIDE (second contributor group to
                        record it, robust to single re-processed deposits). Years with < 15 versioned datasets
                        are not drawn.
33-tool-age-at-release  For each tool with native version files, the age of the exact software release
                        recorded in a deposit when the dataset goes public: release date minus the first PRIDE
                        submission recording the same release ("first seen"). Software can only be older than
                        its first PRIDE appearance, so every age is a lower bound. Quantile dot plots (50 dots
                        = 2% each), medians with 95% group-bootstrap intervals, and each tool's ceiling.

Sources: Arsenal/H9-versions/Dataset-tool-version.csv.gz (native files: MaxQuant parameters.txt/mqpar.xml,
DIA-NN logs, fragpipe.workflow, Skyline .sky, Mascot .dat; one row per dataset x tool) and
Arsenal/H9-versions/Standard-file-software.csv.gz (Proteome Discoverer named as the writer of a deposited
mzIdentML/mzTab file). Release dates come from Deposit-features (cohort: 36,249 datasets released 2013-2025).
A recorded version is what the deposited files say, not an installation or licence record; datasets that
deposit these files are a selected subset (denominators in Gallery/data). No trend in software age is claimed
(ledger verdict H9.2: the trend depends on the anchor). Uncertainty: 95% percentile intervals resampling hashed
pi_group clusters (deposits without a group are their own cluster), 2,000 resamples, seed 20261001.
Plotted tables -> Gallery/data/<NN>-<slug>.csv.
"""
import warnings; warnings.filterwarnings('ignore', category=UserWarning)
import numpy as np, pandas as pd
import matplotlib.patheffects as pe
from matplotlib.ticker import FuncFormatter
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3
from common import O, SEED, cluster_boot

D = G / 'data'; D.mkdir(exist_ok=True)
H9 = A / 'H9-versions'
REPS = 2000
HALO = [pe.withStroke(linewidth=2, foreground='white')]
TOOL_COL = {'MaxQuant': '#2e6f9e', 'Proteome Discoverer': '#8d6cab', 'DIA-NN': '#d1495b', 'FragPipe': '#5a9e6f',
            'Mascot': '#7f8c96', 'Skyline': '#e8a33d'}


def headmm(fig, title, sub, h_mm):
    """headline with offsets fixed in millimetres (as in gallery2_*); warns if a line overflows 181 mm."""
    if record_head(title, sub):
        return
    fig.text(.012, 1 - 1.6 / h_mm, title, fontsize=10.5, fontweight='bold', va='top', color=INK)
    fig.text(.012, 1 - 6.4 / h_mm, sub, fontsize=7, va='top', color=MUTED, linespacing=1.3)
    r = fig.canvas.get_renderer()
    for tx in fig.texts[-2:]:
        wmm = tx.get_window_extent(r).width / fig.dpi * 25.4
        if wmm > 181:
            print(f'WARNING headline text {wmm:.0f} mm wide: {tx.get_text()[:60]}')


def group_boot_median(v, g, reps=REPS, seed=SEED):
    """95% percentile CI of the median, resampling whole clusters g."""
    v = np.asarray(v, float); g = pd.Series(g).astype(str).to_numpy()
    o = np.argsort(g, kind='stable'); v, g = v[o], g[o]
    _, starts, sizes = np.unique(g, return_index=True, return_counts=True)
    rng = np.random.default_rng(seed); k = len(starts); out = np.empty(reps)
    for r in range(reps):
        s = rng.integers(0, k, k); ln = sizes[s]
        idx = np.repeat(starts[s], ln) + (np.arange(ln.sum()) - np.repeat(np.cumsum(ln) - ln, ln))
        out[r] = np.median(v[idx])
    return np.quantile(out, .025), np.quantile(out, .975)


def decyear(ts):
    ts = pd.to_datetime(ts)
    return ts.dt.year + (ts.dt.dayofyear - .5) / 365.25


def load_versions():
    """One row per dataset x tool with a recorded version: native tool files plus PD from standard files."""
    v = pd.read_csv(H9 / 'Dataset-tool-version.csv.gz')
    v = v[~((v.software == 'Skyline') & (v.release.astype(str) == '64'))]          # '(64-bit)' misparsed as a version
    v = v[['accession', 'software', 'version', 'release', 'major_minor', 'submission_date', 'pi_group', 'first_seen']].copy()
    v['major_minor'] = v.major_minor.astype(float)
    s = pd.read_csv(H9 / 'Standard-file-software.csv.gz')
    p = s[s['product'] == 'Proteome Discoverer'].copy()
    p['release'] = p.version.str.extract(r'(\d+\.\d+\.\d+)')[0].fillna(p.version.str.extract(r'(\d+\.\d+)')[0])
    p['major_minor'] = p.version.str.extract(r'(\d+\.\d+)')[0].astype(float)
    p = p.dropna(subset=['major_minor'])
    # several PD entries per dataset (one per mzIdentML node): keep the newest recorded version (conservative
    # for any "old line persists" statement)
    p = p.sort_values('major_minor').drop_duplicates('accession', keep='last')
    p['software'] = 'Proteome Discoverer'
    p = p[['accession', 'software', 'version', 'release', 'major_minor', 'submission_date', 'pi_group']]
    v = pd.concat([v, p], ignore_index=True)
    v['grp'] = v.pi_group.fillna('acc:' + v.accession)
    v['sub_year'] = v.submission_date.str[:4].astype(int)
    return v.merge(F[['accession', 'release_date', 'release_year']], on='accession', how='inner', validate='many_to_one')


LINES = {
    'MaxQuant': [('≤1.4', lambda m: m < 1.5), ('1.5', lambda m: m == 1.5), ('1.6', lambda m: m == 1.6), ('2.0', lambda m: m == 2.0),
                 ('2.1', lambda m: m == 2.1), ('2.2–2.3', lambda m: (m >= 2.2) & (m < 2.4)), ('2.4', lambda m: m == 2.4),
                 ('2.5–2.6', lambda m: m >= 2.5)],
    'Proteome Discoverer': [('1.x', lambda m: m < 2), ('2.0–2.2', lambda m: (m >= 2) & (m < 2.3)), ('2.3', lambda m: m == 2.3),
                            ('2.4', lambda m: m == 2.4), ('2.5', lambda m: m == 2.5), ('3.0', lambda m: m == 3.0),
                            ('3.1–3.2', lambda m: m >= 3.1)],
    'DIA-NN': [('1.7', lambda m: m < 1.8), ('1.8', lambda m: m == 1.8), ('1.9', lambda m: m == 1.9), ('2.x', lambda m: m >= 2)],
    'FragPipe': [('19', lambda m: m < 20), ('20', lambda m: m == 20), ('21', lambda m: (m >= 21) & (m < 22)), ('22', lambda m: m >= 22)],
}


def robust_first_seen(x):
    """Date the line is first recorded by a second distinct contributor group (single early outliers -- e.g. a
    2022 submission re-processed with MaxQuant 2.6 before its 2025 release -- would otherwise backdate it)."""
    f = x.sort_values('submission_date').drop_duplicates('grp').submission_date
    return f.iloc[1] if len(f) > 1 else f.iloc[0]


# --------------------------------------------------------------------------------------------
# 31. Release-line streams
def fig_version_streams():
    v = load_versions()
    tools = ['MaxQuant', 'Proteome Discoverer', 'DIA-NN', 'FragPipe']
    flag = {'MaxQuant': 'maxquant_expanded', 'Proteome Discoverer': 'pd_expanded', 'DIA-NN': 'diann_expanded', 'FragPipe': 'fragpipe_expanded'}
    X0, X1, NMIN = 2015, 2025, 15
    rows, labels25, stats = [], {}, {}
    for t in tools:
        x = v[v.software == t].copy()
        x['line'] = None
        for name, rule in LINES[t]:
            x.loc[rule(x.major_minor) & x.line.isna(), 'line'] = name
        assert x.line.notna().all(), t
        fs = x.groupby('line').apply(robust_first_seen)
        tab = x[x.sub_year.between(X0, X1)].pivot_table(index='sub_year', columns='line', values='accession', aggfunc='count').fillna(0)
        n = tab.sum(axis=1)
        sy = F.submission_date.dt.year
        named = F[F[flag[t]].fillna(False).astype(bool)].groupby(sy).size()
        for yr in tab.index:
            for name, _ in LINES[t]:
                c = tab.loc[yr, name] if name in tab else 0
                rows.append({'tool': t, 'sub_year': yr, 'line': name, 'datasets': int(c), 'versioned_datasets': int(n[yr]),
                             'datasets_naming_tool': int(named.get(yr, 0)), 'share_pct': round(100 * c / n[yr], 2),
                             'drawn': bool(n[yr] >= NMIN), 'line_first_seen_2nd_group': str(pd.Timestamp(fs[name]).date())})
        x25 = x[x.sub_year == 2025]
        stats[t] = {'n25': len(x25), 'named25': int(named.get(2025, 0)), 'fs': fs, 'groups25': x25.grp.nunique(), 'x': x}
    out = pd.DataFrame(rows)

    def share_ci(t, rule, yr=2025):
        x = v[(v.software == t) & (v.sub_year == yr)]
        k = rule(x.major_minor).astype(float).to_numpy()
        lo, hi = cluster_boot(k, x.grp.to_numpy(), reps=REPS, seed=SEED)
        return 100 * k.mean(), 100 * lo, 100 * hi, len(x), int(k.sum()), x[rule(x.major_minor)].grp.nunique()
    mq16 = share_ci('MaxQuant', lambda m: m == 1.6)
    mq1 = share_ci('MaxQuant', lambda m: m < 2)
    dn18 = share_ci('DIA-NN', lambda m: m == 1.8)
    pd2 = share_ci('Proteome Discoverer', lambda m: m < 3)
    fp22 = share_ci('FragPipe', lambda m: m >= 22)
    mq16_rel = v[(v.software == 'MaxQuant') & (v.release_year == 2025)]
    mq16_rel = (100 * (mq16_rel.major_minor == 1.6).mean(), len(mq16_rel))
    fs16 = pd.Timestamp(stats['MaxQuant']['fs']['1.6'])
    fs16_first = pd.Timestamp(v[(v.software == 'MaxQuant') & (v.major_minor == 1.6)].submission_date.min())
    yrs_on = 2025 - fs16_first.year
    pd.DataFrame([dict(zip(['metric', 'share_pct', 'ci_lo', 'ci_hi', 'versioned_datasets_2025', 'datasets_on_line', 'groups_on_line'],
                           [k, *vals])) for k, vals in [('MaxQuant 1.6', mq16), ('MaxQuant 1.x', mq1), ('DIA-NN 1.8', dn18),
                                                         ('Proteome Discoverer 2.x', pd2), ('FragPipe 22', fp22)]]).round(2) \
        .to_csv(D / '31-version-streams-2025-shares.csv', index=False)
    out.to_csv(D / '31-version-streams.csv', index=False)

    # e: MaxQuant 1.x share by data age at submission (latest raw-file acquisition date in file names, Arsenal/H11)
    clk = pd.read_csv(A / 'H11-clocks/Dataset-clocks.csv.gz', usecols=['accession', 'acq_to_sub'])
    z = v[(v.software == 'MaxQuant') & v.sub_year.between(2023, 2025)].merge(clk, on='accession', how='left', validate='one_to_one')
    z['old'] = z.major_minor < 2
    z['dage'] = pd.cut(z.acq_to_sub / 365.25, [-np.inf, 1, 2, 3, np.inf], labels=['< 1 year', '1–2 years', '2–3 years', '> 3 years'])
    z['dage'] = z.dage.cat.add_categories('no dated raw files').fillna('no dated raw files').astype(str)
    erows = []
    for kk in ['< 1 year', '1–2 years', '2–3 years', '> 3 years', 'no dated raw files', 'all']:
        q = z if kk == 'all' else z[z.dage == kk]
        lo, hi = cluster_boot(q.old.astype(float).to_numpy(), q.grp.to_numpy(), reps=REPS, seed=SEED)
        erows.append({'data_age_at_submission': kk, 'datasets': len(q), 'groups': q.grp.nunique(), 'on_1x': int(q.old.sum()),
                      'share_pct': 100 * q.old.mean(), 'ci_lo': 100 * lo, 'ci_hi': 100 * hi})
    E = pd.DataFrame(erows); E.round(2).to_csv(D / '31-version-streams-maxquant-by-data-age.csv', index=False)

    H = 128
    fig = plt.figure(figsize=(W, H * MM))
    words = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}
    headmm(fig, f'Old release lines linger: {mq16[0]:.0f}% of 2025 MaxQuant deposits ran 1.6, first seen {words[yrs_on]} years earlier',
           'Share of PRIDE datasets recording each software release line, by submission year (datasets released 2013–2025). Versions are read from deposited\n'
           'parameter files, logs and workflows (MaxQuant, DIA-NN, FragPipe) or mzIdentML writer fields (Proteome Discoverer). Black triangles: a second\n'
           f'contributor group first records the line in PRIDE. Grey numbers: datasets with a version. 2025 coverage: MaxQuant {stats["MaxQuant"]["n25"]} of '
           f'{stats["MaxQuant"]["named25"]} datasets\nnaming the tool, Proteome Discoverer {stats["Proteome Discoverer"]["n25"]} of {stats["Proteome Discoverer"]["named25"]}, '
           f'DIA-NN {stats["DIA-NN"]["n25"]} of {stats["DIA-NN"]["named25"]}, FragPipe {stats["FragPipe"]["n25"]} of {stats["FragPipe"]["named25"]}.', H)
    U = .0315                                                    # figure width per year: identical time scale in every panel
    Y0 = {'MaxQuant': 2015, 'Proteome Discoverer': 2016, 'DIA-NN': 2022, 'FragPipe': 2023}
    left = {'MaxQuant': .06, 'Proteome Discoverer': .545, 'DIA-NN': .06, 'FragPipe': .285}
    bottom = {'MaxQuant': .485, 'Proteome Discoverer': .485, 'DIA-NN': .075, 'FragPipe': .075}
    AH = .28
    letters = dict(zip(tools, 'abcd'))
    for t in tools:
        y0 = Y0[t]
        ax = fig.add_axes([left[t], bottom[t], U * (2025 - y0), AH])
        o = out[(out.tool == t) & out.drawn & (out.sub_year >= y0)]
        tab = o.pivot(index='sub_year', columns='line', values='share_pct')[[l for l, _ in LINES[t]]]
        yrs = tab.index.to_numpy(); assert yrs[0] == y0, (t, yrs)
        fs = stats[t]['fs']
        fsd = {l: float(decyear(pd.Series([fs[l]])).iloc[0]) for l in tab.columns}
        xm = yrs + .5
        grid = np.array(sorted(set(xm) | {f for f in fsd.values() if xm[0] < f < xm[-1]}))
        Y = np.zeros((len(tab.columns), len(grid)))
        for i, l in enumerate(tab.columns):
            f = fsd[l]
            y = np.interp(grid, xm, tab[l].to_numpy())
            if f > xm[0]:
                nxt = xm[xm >= f]
                if len(nxt):
                    y = np.where(grid < f, 0, np.where(grid <= nxt[0], np.interp(grid, [f, nxt[0]], [0, tab[l].loc[int(nxt[0] - .5)]]), y))
            Y[i] = y
        Y = 100 * Y / Y.sum(0, keepdims=True)
        k = len(tab.columns); base = TOOL_COL[t]
        cols = [mix(base, i, k) for i in range(k)]
        cum = np.vstack([np.zeros(len(grid)), np.cumsum(Y, 0)])
        for i in range(k):
            ax.fill_between(grid, cum[i], cum[i + 1], color=cols[i], lw=0)
            ax.plot(grid, cum[i + 1], color='white', lw=.5)
        for i, l in enumerate(tab.columns):
            f = fsd[l]
            if xm[0] < f < xm[-1]:
                ax.plot(f, np.interp(f, grid, cum[i + 1]) if np.interp(f, grid, cum[i]) > 99 else np.interp(f, grid, cum[i]),
                        marker='v', ms=3.4, color=INK, mec='white', mew=.4, zorder=5, clip_on=False)
        last = Y[:, -1]; mid = (cum[:-1, -1] + cum[1:, -1]) / 2
        lab = repel({l: mid[i] for i, l in enumerate(tab.columns) if last[i] >= .5}, 7.6)
        lc = darken(base, .62)
        for i, l in enumerate(tab.columns):
            if l not in lab:
                continue
            bold = 'bold' if last[i] == last.max() else 'normal'
            ax.annotate(f'{l}  {last[i]:.0f}%', xy=(grid[-1], mid[i]), xytext=(grid[-1] + .32, lab[l]), textcoords='data',
                        fontsize=6.3, va='center', ha='left', color=lc, fontweight=bold, annotation_clip=False,
                        arrowprops=dict(arrowstyle='-', color='#9a9890', lw=.45, shrinkA=0, shrinkB=1))
        for yr in yrs:
            nn = int(o[o.sub_year == yr].versioned_datasets.iloc[0])
            ax.text(yr + .5, 102.5, f'{nn / 1000:.1f}k' if nn >= 1000 else f'{nn}', ha='center', va='bottom', fontsize=6, color=MUTED)
        ax.set_xlim(xm[0], xm[-1]); ax.set_ylim(0, 100)
        ax.set_xticks(xm, minor=True)
        tk = [y for y in yrs if (y - 2015) % 2 == 0] if len(yrs) > 5 else list(yrs)
        ax.set_xticks([y + .5 for y in tk]); ax.set_xticklabels([f'{y}' if len(yrs) > 5 else f"'{y % 100:02d}" for y in tk])
        ax.tick_params(axis='x', which='minor', length=1.5)
        ax.set_yticks([0, 25, 50, 75, 100]); yfmt(ax)
        if t != 'MaxQuant' and t != 'DIA-NN':
            ax.set_yticklabels([])
        st = stats[t]
        ax.text(0, 1.115, t, transform=ax.transAxes, fontsize=7.5, fontweight='bold', color=darken(base, .75), va='bottom', ha='left')
        panel(ax, letters[t], x=-.035 / (U * (2025 - y0)) * .32 if False else -.012 / (U * (2025 - y0)), y=1.115)

    # e: MaxQuant 1.x by data age
    ex = fig.add_axes([.62, .075, .285, AH])
    rows_e = ['< 1 year', '1–2 years', '2–3 years', '> 3 years', 'no dated raw files']
    c = TOOL_COL['MaxQuant']
    for j, kk in enumerate(rows_e):
        r = E[E.data_age_at_submission == kk].iloc[0]; yy = len(rows_e) - 1 - j + (0 if j < 4 else -.25)
        cc = c if j < 4 else MUTED
        ex.plot([r.ci_lo, r.ci_hi], [yy] * 2, color=cc, lw=1.6, solid_capstyle='butt')
        ex.plot(r.share_pct, yy, 'o', ms=4.6, color=cc, mec='white', mew=.6, zorder=3)
        ex.text(r.ci_hi + 2.5, yy, f'{r.share_pct:.0f}%', fontsize=6.5, fontweight='bold', color=darken(cc, .75), va='center')
        ex.text(-3, yy + .1, kk if j < 4 else 'no dated\nraw files', ha='right', va='center', fontsize=6.3, color=INK if j < 4 else MUTED,
                linespacing=1.0)
        ex.text(-3, yy - (.28 if j < 4 else .42), f'n = {int(r.datasets):,}', ha='right', va='center', fontsize=6, color=MUTED)
    ex.set_yticks([]); ex.set_xlim(0, 103); ex.set_ylim(-.85, len(rows_e) - .25)
    ex.set_xticks([0, 25, 50, 75, 100]); ex.xaxis.set_major_formatter(FuncFormatter(lambda v_, _: f'{v_:.0f}%'))
    ex.grid(axis='x', color=HAIR, lw=.4); ex.set_axisbelow(True)
    ex.axhline(.375, color=HAIR, lw=.6)
    ex.text(-.27, 1.115, 'Old lines travel with old data', transform=ex.transAxes, fontsize=7.5, fontweight='bold', color=INK, va='bottom', ha='left')
    ex.text(-.27, 1.015, 'MaxQuant deposits 2023–25, by time from newest raw file to submission', transform=ex.transAxes, fontsize=6,
            color=MUTED, va='bottom', ha='left')
    ex.set_xlabel('Share on MaxQuant 1.x (≤1.6), 95% CI', fontsize=6.5, labelpad=2)
    panel(ex, 'e', x=-.29, y=1.115)
    footnote(fig, .012, .004, 'a–d: submission (deposit) year on one time scale; 2025 = submissions public by 31 Dec 2025. '
                         'e: 95% group-bootstrap intervals over hashed contributor groups.', fontsize=6, color=MUTED)
    save(fig, '31-version-streams')
    return dict(mq16=mq16, mq1=mq1, dn18=dn18, pd2=pd2, fp22=fp22, mq16_rel=mq16_rel, fs16_first=fs16_first, fs16=fs16, stats=stats)


def mix(hexc, i, k):
    """Sequential ramp within one hue: oldest line darkest, newest lightest."""
    c = np.array(matplotlib.colors.to_rgb(hexc))
    dark = c * .55; light = c + (1 - c) * .78
    f = i / max(k - 1, 1)
    return tuple(dark * (1 - f) + light * f)


def darken(c, f=.8):
    return tuple(np.array(matplotlib.colors.to_rgb(c)) * f)


# --------------------------------------------------------------------------------------------
# 33. Software-release age at public release
def fig_tool_age():
    v = load_versions()
    v = v[v.software.isin(['MaxQuant', 'Mascot', 'DIA-NN', 'FragPipe', 'Skyline'])].copy()
    v['rel'] = pd.to_datetime(v.release_date); v['fs'] = pd.to_datetime(v.first_seen)
    v['sd'] = pd.to_datetime(v.submission_date)
    v['age'] = (v.rel - v.fs).dt.days / 365.25
    # sensitivity anchor: first submission by a second distinct contributor group recording the exact release
    a2 = v.groupby(['software', 'release']).apply(lambda x: x.sort_values('sd').drop_duplicates('grp').sd.iloc[min(1, x.grp.nunique() - 1)]).rename('fs2')
    v = v.join(a2, on=['software', 'release'])
    v['age2'] = (v.rel - v.fs2).dt.days / 365.25
    neg = int((v.age < 0).sum())
    v = v[v.age >= 0]                                   # release date before first-seen submission (PRIDE date quirk)
    ceil = {t: (pd.Timestamp('2025-12-31') - v[v.software == t].fs.min()).days / 365.25 for t in v.software.unique()}
    res = []
    for t, x in v.groupby('software'):
        lo, hi = group_boot_median(x.age, x.grp)
        r = x[x.release_year >= 2023]
        lo2, hi2 = group_boot_median(r.age, r.grp)
        res.append({'tool': t, 'datasets': len(x), 'groups': x.grp.nunique(), 'median_years': x.age.median(), 'ci_lo': lo, 'ci_hi': hi,
                    'q10': x.age.quantile(.1), 'q25': x.age.quantile(.25), 'q75': x.age.quantile(.75), 'q90': x.age.quantile(.9),
                    'pct_gt2y': 100 * (x.age > 2).mean(), 'median_2nd_group_anchor': x.age2.median(),
                    'median_at_submission': ((x.sd - x.fs).dt.days / 365.25).median(),
                    'datasets_2023_25': len(r), 'median_2023_25': r.age.median(), 'ci_lo_2023_25': lo2, 'ci_hi_2023_25': hi2,
                    'ceiling_years_end_2025': ceil[t], 'tool_first_seen': str(x.fs.min().date())})
    res = pd.DataFrame(res).sort_values('median_years').reset_index(drop=True)
    res.round(3).to_csv(D / '33-tool-age-at-release.csv', index=False)
    qd = []
    for t, x in v.groupby('software'):
        q = np.quantile(x.age, (np.arange(50) + .5) / 50)
        qd += [{'tool': t, 'dot': i + 1, 'quantile': (i + .5) / 50, 'age_years': a} for i, a in enumerate(q)]
    qd = pd.DataFrame(qd); qd.round(3).to_csv(D / '33-tool-age-at-release-quantile-dots.csv', index=False)

    H = 108
    fig = plt.figure(figsize=(W, H * MM))
    m = res.set_index('tool')
    headmm(fig, f'When a dataset goes public, its Mascot release is typically ≥{m.median_years["Mascot"]:.1f} years old; '
                f'its Skyline release ≥{m.median_years["Skyline"]:.1f}',
           'PRIDE datasets released 2013–2025 whose deposited files record an exact software release (one value per dataset and tool). Age = public release date\n'
           'minus the first PRIDE submission recording the same release. Software exists on or before that date, so every age is a lower bound.', H)
    ax = fig.add_axes([.165, .155, .6, .685])
    order = list(res.tool)                                    # youngest at top
    ypos = {t: len(order) - 1 - i for i, t in enumerate(order)}
    XMAX = 13.6; bw = .3
    stacks = {t: np.unique(np.floor(qd[qd.tool == t].age_years / bw), return_counts=True)[1].max() for t in order}
    step = min(.075, .6 / max(stacks.values()))
    for t in order:
        y = ypos[t]; c = TOOL_COL[t]; r = m.loc[t]
        q = qd[qd.tool == t].age_years.to_numpy()
        bins = np.floor(q / bw); cnt = {}
        xs, ys = [], []
        for a_, b_ in zip(q, bins):
            j = cnt.get(b_, 0); cnt[b_] = j + 1
            xs.append((b_ + .5) * bw); ys.append(y + .03 + j * step)
        ax.scatter(xs, ys, s=6, color=c, lw=0, zorder=3)
        yb = y - .12
        ax.plot([r.q10, r.q90], [yb, yb], color=INK, lw=.6, ls=(0, (1, 1.5)), zorder=3)
        ax.plot([r.q25, r.q75], [yb, yb], color=INK, lw=.7, zorder=3)
        ax.plot([r.ci_lo, r.ci_hi], [yb, yb], color=INK, lw=2.4, solid_capstyle='butt', zorder=4)
        ax.plot(r.median_years, yb, marker='|', ms=7, mew=1.4, color=INK, zorder=5)
        lx = (r.ceiling_years_end_2025 if r.ceiling_years_end_2025 < 8 else r.q90) + .2
        ax.text(lx, yb, f'{r.median_years:.1f} y', fontsize=7, fontweight='bold', va='center', ha='left', color=INK, zorder=6)
        ax.text(lx + .62, yb, f'  ({r.median_at_submission:.1f} y at submission)', fontsize=6.2, va='center', ha='left', color=MUTED)
        yd = y - .27
        ax.plot([r.ci_lo_2023_25, r.ci_hi_2023_25], [yd, yd], color=INK, lw=.7, zorder=4)
        ax.plot(r.median_2023_25, yd, marker='D', ms=3, mfc='white', mec=INK, mew=.7, zorder=5)
        ax.plot([r.ceiling_years_end_2025] * 2, [y - .32, y + .62], color=INK, lw=.8, ls=(0, (2.5, 1.5)), zorder=2)
    ax.set_yticks([ypos[t] for t in order])
    ax.set_yticklabels([f'{t}\n{int(m.datasets[t]):,} datasets\n{int(m.groups[t]):,} groups' for t in order], fontsize=6.5, linespacing=1.2)
    for lab, t in zip(ax.get_yticklabels(), order):
        lab.set_color(darken(TOOL_COL[t], .75))
    ax.tick_params(axis='y', length=0)
    ax.set_ylim(-.45, len(order) - .3); ax.set_xlim(0, XMAX)
    ax.set_xticks(range(0, 14, 2)); ax.set_xlabel('Years from first PRIDE appearance of the release to public release of the dataset')
    ax.grid(axis='x', color=HAIR, lw=.4); ax.set_axisbelow(True)
    kx = fig.add_axes([.79, .155, .2, .685]); kx.axis('off'); kx.set_xlim(0, 1); kx.set_ylim(0, 1)
    yk = .99
    kx.text(0, yk, 'How to read', fontsize=7, fontweight='bold', va='top', color=INK)
    kx.scatter([.03, .07, .11, .07], [yk - .1, yk - .1, yk - .1, yk - .1 + .03], s=6, color=MUTED, lw=0)
    kx.text(.2, yk - .095, 'each dot = 2% of the\ntool’s datasets', fontsize=6.2, va='center', color=INK, linespacing=1.15)
    kx.plot([.02, .16], [yk - .21] * 2, color=INK, lw=2.4, solid_capstyle='butt'); kx.plot(.09, yk - .21, marker='|', ms=7, mew=1.4, color=INK)
    kx.text(.2, yk - .21, 'median, 95% CI\n(group bootstrap)', fontsize=6.2, va='center', color=INK, linespacing=1.15)
    kx.plot([.0, .18], [yk - .31] * 2, color=INK, lw=.6, ls=(0, (1, 1.5))); kx.plot([.04, .14], [yk - .31] * 2, color=INK, lw=.7)
    kx.text(.2, yk - .31, 'IQR; 10th–90th\npercentile', fontsize=6.2, va='center', color=INK, linespacing=1.15)
    kx.plot([.02, .16], [yk - .41] * 2, color=INK, lw=.7); kx.plot(.09, yk - .41, marker='D', ms=3, mfc='white', mec=INK, mew=.7)
    kx.text(.2, yk - .41, 'median, datasets\nreleased 2023–2025', fontsize=6.2, va='center', color=INK, linespacing=1.15)
    kx.plot([.09, .09], [yk - .56, yk - .47], color=INK, lw=.8, ls=(0, (2.5, 1.5)))
    kx.text(.2, yk - .515, 'ceiling: oldest release\nof the tool in PRIDE,\nto end of 2025', fontsize=6.2, va='center', color=INK, linespacing=1.15)
    kx.text(0, yk - .63, f'FragPipe and DIA-NN ages are\ncapped by their short history\nin PRIDE (ceilings {m.ceiling_years_end_2025["FragPipe"]:.1f} and '
                        f'{m.ceiling_years_end_2025["DIA-NN"]:.1f} y).\n\nDiamonds are a composition\ncheck, not a trend: later\nwindows admit older releases.',
            fontsize=6.2, va='top', color=MUTED, linespacing=1.2)
    a2s = ', '.join(f'{t} {m.median_2nd_group_anchor[t]:.1f}' for t in order)
    footnote(fig, .012, .004, f'Sensitivity: anchoring each release at the second contributor group to record it gives medians {a2s} y.\n'
                         f'Excluded: {neg} datasets whose PRIDE release date precedes the first-seen submission; 4 Skyline files whose "version" was "64-bit".',
             fontsize=6, color=MUTED, linespacing=1.25)
    save(fig, '33-tool-age-at-release')
    return res, neg


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['31', '33']
    if '31' in which:
        s = fig_version_streams()
        for k in ['mq16', 'mq1', 'dn18', 'pd2', 'fp22']:
            print(k, np.round(s[k][:3], 1), s[k][3:])
        print('MQ 1.6 by release year 2025', s['mq16_rel'], 'first seen', s['fs16_first'].date(), '2nd group', s['fs16'].date())
        for t, st in s['stats'].items():
            print(t, st['n25'], st['named25'], st['groups25'])
    if '33' in which:
        res, neg = fig_tool_age()
        print(res.round(2).to_string()); print('negative', neg)
