"""Headline figures: one striking, verified finding per figure (descriptive).

Every number drawn here comes from frozen outputs under outputs/Stack-and-reuse (ledger IDs in each
figure's note) or from Deposit-features.parquet. Figures are 183 mm wide, text >= 6 pt; PNG (300 dpi) and PDF.
Outputs: outputs/Stack-and-reuse/Gallery/.
"""
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from common import O
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3, VEND


# 1. Instrument kings: each leading model peaks lower than the last; the fleet diversifies
def fig_instruments():
    m = pd.read_csv(A / 'H1-H3-scale-instruments/H3-model-by-year.csv')
    div = pd.read_csv(A / 'H1-H3-scale-instruments/H3-fleet-diversity-by-year.csv').set_index('release_year')
    kings = ['LTQ Orbitrap', 'LTQ Orbitrap Velos', 'Q Exactive', 'Q Exactive HF', 'Orbitrap Fusion Lumos', 'Orbitrap Exploris 480', 'timsTOF Pro', 'Orbitrap Astral']
    p = m.pivot_table(index='release_year', columns='model', values='credit_share_labelled', aggfunc='sum').reindex(YEARS).fillna(0) * 100
    fig = plt.figure(figsize=(W, 72 * MM))
    headline(fig, 'New instrument models joined a persistent record of earlier technologies',
             'a, Share of each year’s model-labelled dataset credit assigned to selected instrument models, with equal fractional credit when a dataset names several models. Every model named in ≥50 datasets still appeared in 2025 releases; deposits include delayed releases and reanalyses. b, Effective number of models (inverse Simpson index of all fractional model shares); shading: 95% intervals from 300 bootstrap resamples of contributor groups.')
    ax = fig.add_axes([.06, .14, .6, .66])
    cols = ['#9db4c8', '#6f93b3', '#2e6f9e', '#1d4f75', '#5a9e6f', '#e8a33d', '#d1495b', '#8d6cab']
    for k, c in zip(kings, cols):
        y = p[k]; ax.plot(YEARS, y, color=c, lw=1.6)
        pk = y.idxmax(); ax.plot(pk, y[pk], 'o', ms=3, color=c)
        lab = k.replace('Orbitrap ', '').replace('LTQ ', 'LTQ ')
        if k in ['Orbitrap Exploris 480', 'timsTOF Pro', 'Orbitrap Astral']:
            yy = {'Orbitrap Exploris 480': 12, 'timsTOF Pro': 6.6, 'Orbitrap Astral': 2.6}[k]
            ax.text(2025.25, yy, f'{lab} {y[2025]:.0f}%', ha='left', va='center', fontsize=6, color=c, fontweight='bold')
            continue
        off = {'Q Exactive': (0, 1.6), 'LTQ Orbitrap': (.35, 1.2), 'LTQ Orbitrap Velos': (0, 1.4), 'Q Exactive HF': (-1.2, 4.2),
               'Orbitrap Fusion Lumos': (-2.6, 3.2)}[k]
        ax.text(pk + off[0], y[pk] + off[1], f'{lab}\n{y[pk]:.0f}%', ha='center', fontsize=6, color=c, fontweight='bold', linespacing=.95)
    ax.set_xlim(2012.6, 2027.3); ax.set_xticks(range(2014, 2026, 2)); ax.set_ylim(0, 47); yfmt(ax); ax.set_ylabel('Share of model-labelled credit')
    ax2 = fig.add_axes([.74, .14, .24, .66])
    e = div.loc[YEARS]
    ax2.fill_between(YEARS, e.eff_ci_lo, e.eff_ci_hi, color=ACC2, alpha=.18, lw=0)
    ax2.plot(YEARS, e.effective_models_inv_simpson, color=ACC2, lw=1.8)
    for yr in [2013, 2025]:
        v = e.effective_models_inv_simpson[yr]; ax2.text(yr + (.9 if yr == 2013 else -.6), v + (2.6 if yr == 2013 else 1), f'{v:.1f}', ha='center', fontsize=8, fontweight='bold', color=ACC2)
    ax2.set_ylim(0, 15); ax2.set_xlim(2012.5, 2025.5); ax2.set_xticks([2013, 2019, 2025])
    ax2.set_title('Effective number of models', fontsize=7, loc='left', color=INK)
    panel(ax, 'a'); panel(ax2, 'b', x=-.16)
    save(fig, '01-instruments')


# 2. Runs got bigger with every generation
def fig_run_size():
    r = pd.read_csv(A / 'H1-H3-scale-instruments/H1-per-run-size-by-model.csv')
    vend = {'QSTAR': 'SCIEX', 'TripleTOF 5600': 'SCIEX', 'TripleTOF 6600': 'SCIEX', 'maXis': 'Bruker', 'impact II': 'Bruker',
            'timsTOF HT': 'Bruker', 'timsTOF Pro 2': 'Bruker', 'timsTOF Pro': 'Bruker', 'Synapt MS': 'Waters'}
    r['vendor'] = r.model.map(vend).fillna('Thermo Fisher')
    first = m_first()
    r['launch'] = r.model.map(first).fillna(r.first_year)
    fig = plt.figure(figsize=(W, 70 * MM))
    headline(fig, 'One Orbitrap Astral run holds 9\u00d7 the data of an LTQ Orbitrap run',
             'Median size of one raw run by instrument model (bars: interquartile range), placed at the model’s first year in PRIDE.')
    ax = fig.add_axes([.07, .14, .9, .66])
    rng = np.random.default_rng(3)
    for _, x in r.iterrows():
        xx = x.launch + rng.uniform(-.25, .25)
        c = VEND[x.vendor]
        ax.plot([xx, xx], [x.p25, x.p75], color=c, lw=.8, alpha=.5)
        ax.scatter(xx, x.median_per_run_GB, s=8 + np.sqrt(x.datasets) * 1.6, color=c, edgecolor='white', lw=.5, zorder=3)
        if x.model in ['QSTAR', 'LTQ', 'LTQ Orbitrap', 'Q Exactive', 'Orbitrap Fusion Lumos', 'Q Exactive HF-X', 'timsTOF Pro',
                       'Orbitrap Astral', 'Orbitrap Exploris 480', 'Synapt MS', 'timsTOF HT']:
            dx, dy = {'Q Exactive': (.3, 1.45), 'Orbitrap Fusion Lumos': (.1, .55), 'Orbitrap Exploris 480': (.2, .74), 'LTQ Orbitrap': (.3, .95),
                      'TripleTOF 5600': (.3, 1.25), 'Synapt MS': (.3, 1), 'LTQ': (.3, 1), 'QSTAR': (.3, 1), 'Q Exactive HF-X': (.25, 1.15),
                      'timsTOF Pro': (.25, 1), 'Orbitrap Astral': (.25, 1), 'timsTOF HT': (.25, .85)}[x.model]
            ax.text(xx + dx, x.median_per_run_GB * dy, x.model, fontsize=6, va='center', color=c, ha='right' if dx < 0 else 'left')
    ax.set_yscale('log'); ax.set_ylim(.04, 12)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g} GB'))
    ax.set_xlim(2012.3, 2026.2); ax.set_xlabel('First year the model appears in PRIDE')
    ax.set_ylabel('Median size of one raw run')
    for v, c in VEND.items():
        if v in r.vendor.values:
            ax.scatter([], [], color=c, s=20, label=v)
    ax.legend(loc='upper left', ncol=4, handletextpad=.2)
    save(fig, '02-run-size')


def m_first():
    m = pd.read_csv(A / 'H1-H3-scale-instruments/H3-model-by-year.csv')
    m = m[m.datasets >= 3]
    return m.groupby('model').release_year.min()


# 3. More big studies, smaller typical study
def fig_scale():
    s = pd.read_csv(A / 'H1-H3-scale-instruments/H1-scale-by-year.csv').set_index('release_year')
    ci = pd.read_csv(A / 'H1-H3-scale-instruments/H1-year-medians-with-CI.csv').set_index('release_year')
    big = (s.share_ge100_runs * s.datasets_with_raw).round()
    fig = plt.figure(figsize=(W, 68 * MM))
    headline(fig, 'Nine times as many large studies are deposited, yet the typical study shrank by half',
             'PRIDE datasets by release year. Runs = raw MS acquisitions deposited; datasets with vendor raw files only.')
    ax = fig.add_axes([.06, .14, .41, .64])
    ax.bar(YEARS, big[YEARS], color=ACC2, width=.72)
    for yr in [2013, 2025]:
        ax.text(yr, big[yr] + 12, f'{big[yr]:.0f}', ha='center', fontsize=8, fontweight='bold', color=ACC2)
    ax.set_title('Datasets with ≥100 raw runs', loc='left', fontsize=7.5)
    ax.set_xticks([2013, 2016, 2019, 2022, 2025]); ax.set_ylim(0, 600)
    ax2 = fig.add_axes([.57, .14, .41, .64])
    ax2.fill_between(YEARS, ci.ci_lo[YEARS], ci.ci_hi[YEARS], color=ACC, alpha=.15, lw=0)
    ax2.plot(YEARS, ci.median_runs_lb[YEARS], color=ACC, lw=1.8, marker='o', ms=2.5)
    for yr, dy in [(2013, 4), (2016, 4), (2025, 4)]:
        ax2.text(yr, ci.median_runs_lb[yr] + dy, f'{ci.median_runs_lb[yr]:.0f}', ha='center', fontsize=8, fontweight='bold', color=ACC)
    ax2.set_ylim(0, 70); ax2.set_xticks([2013, 2016, 2019, 2022, 2025])
    ax2.set_title('Median raw runs per dataset (95% CI)', loc='left', fontsize=7.5)
    panel(ax, 'a', y=1.06); panel(ax2, 'b', y=1.06)
    save(fig, '03-scale')


# 4. Newcomers: each year's deposits, layered by the year the group first deposited
def fig_newcomers():
    f = F[F.release_year.between(2013, 2025)].copy()
    f['entry'] = f.group_first_year.fillna(f.release_year).clip(2013, 2025).astype(int)
    t = f.pivot_table(index='release_year', columns='entry', values='accession', aggfunc='count').reindex(YEARS).fillna(0)
    first = np.array([t.loc[y, y] for y in YEARS]); tot = t.sum(1).to_numpy()
    fig = plt.figure(figsize=(W, 72 * MM))
    headline(fig, 'Four in ten public datasets come from groups depositing for the first time',
             'PRIDE datasets per release year, layered by the year the contributing group (hashed PI contact) first deposited.')
    ax = fig.add_axes([.07, .13, .6, .67])
    cmap = plt.get_cmap('viridis')
    bottom = np.zeros(len(YEARS))
    for i, e in enumerate(YEARS):
        v = t[e].to_numpy() if e in t else np.zeros(len(YEARS))
        ax.bar(YEARS, v, bottom=bottom, color=cmap(i / 12), width=.82, lw=0)
        bottom += v
    ax.set_xlim(2012.4, 2025.6); ax.set_ylabel('Datasets released')
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
    for e in [2013, 2016, 2019, 2022, 2025]:
        pass
    ax.text(2013, tot[-1] * .93, 'Top band each year: groups\ndepositing for the first time', fontsize=6.5, color=MUTED, va='top')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(2013, 2025))
    cb = fig.colorbar(sm, ax=ax, fraction=.03, pad=.01); cb.set_label('First-deposit year', fontsize=6.5); cb.outline.set_visible(False)
    ax2 = fig.add_axes([.8, .13, .18, .67])
    share = 100 * first / tot
    ax2.plot(YEARS[3:], share[3:], color=ACC, lw=1.8)
    ax2.set_ylim(0, 70); yfmt(ax2); ax2.set_xticks([2016, 2020, 2025])
    ax2.text(2024.2, share[-1] + 4, f'{share[-1]:.0f}%', ha='center', fontsize=8, fontweight='bold', color=ACC)
    ax2.set_title('From first-year groups', loc='left', fontsize=7)
    panel(ax, 'a'); panel(ax2, 'b', x=-.28, y=1.06)
    save(fig, '04-newcomers')
    return share


# 5. Software race
def fig_software():
    t = pd.read_csv(O / 'Figure-data/Tool-year.csv')
    t = t[t.release_year.astype(str).str.fullmatch(r'\d{4}')].astype({'release_year': int, 'deposits': float})
    n = F.groupby('release_year').size()
    t['share'] = 100 * t.deposits / t.release_year.map(n)
    p = t.pivot_table(index='release_year', columns='tool', values='share').reindex(YEARS).loc[2015:]
    fig = plt.figure(figsize=(W, 70 * MM))
    headline(fig, 'DIA-NN went from a handful of public datasets to the third most-named tool in five years',
             'Share of PRIDE datasets naming each analysis tool (protocol text or software field), by release year.')
    ax = fig.add_axes([.06, .13, .78, .67])
    cols = {'MaxQuant': '#2e6f9e', 'Proteome Discoverer': '#8d6cab', 'Mascot': '#9aa5ad', 'DIA-NN': ACC, 'Spectronaut': ACC3, 'FragPipe': '#5a9e6f'}
    pos = repel({k: p[k][2025] for k in cols}, 2.6)
    for k, c in cols.items():
        y = p[k]; lw = 2.2 if k in ['DIA-NN'] else 1.4
        ax.plot(p.index, y, color=c, lw=lw)
        ax.text(2025.2, pos[k], f'{k}  {y[2025]:.0f}%', color=c, fontsize=6.8, va='center', fontweight='bold' if k == 'DIA-NN' else 'normal')
    n_d = t[(t.tool == 'DIA-NN')].set_index('release_year').deposits
    ax.annotate(f'{n_d[2020]:.0f} datasets in 2020', (2020, p['DIA-NN'][2020]), (2018.6, 7), fontsize=6.3, color=ACC,
                arrowprops=dict(arrowstyle='-', color=ACC, lw=.6))
    ax.annotate(f'{n_d[2025]:,.0f} in 2025', (2025, p['DIA-NN'][2025]), (2023.1, 17), fontsize=6.3, color=ACC,
                arrowprops=dict(arrowstyle='-', color=ACC, lw=.6))
    ax.set_xlim(2015, 2025); ax.set_ylim(0, 48); yfmt(ax)
    save(fig, '05-software')


# 6. Standard result-file deposition varies by software
def fig_standards():
    s = pd.read_csv(O / 'Tool-outputs/Standard-rate-by-tool-single-tool-datasets.csv')
    s = s[s.datasets >= 60].sort_values('standard_pct')
    col = {'built-in': ACC2, 'native': '#5a9e6f', 'none': ACC, 'unknown': '#b9b9b9'}
    lab = {'native': 'Writes it by default', 'built-in': 'Can export (an option)', 'none': 'Cannot export', 'unknown': 'Export not documented'}
    fig = plt.figure(figsize=(W, 82 * MM))
    headline(fig, f'Standard result-file deposition varies with the software named: from {s.standard_pct.min():.1f}% to {s.standard_pct.max():.0f}%',
             'Single-tool PRIDE datasets with an mzIdentML or mzTab file (corrected definition), by analysis tool; tools with ≥60 datasets.')
    ax = fig.add_axes([.2, .1, .55, .73])
    y = np.arange(len(s))
    ax.barh(y, s.standard_pct, color=[col[e] for e in s.export], height=.7)
    for i, (v, nn) in enumerate(zip(s.standard_pct, s.datasets)):
        ax.text(v + .8, i, f'{v:.1f}%', va='center', fontsize=6.3)
    ax.set_yticks(y); ax.set_yticklabels(s.tool, fontsize=6.5); ax.set_xlim(0, 86); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%'))
    ax.tick_params(axis='y', length=0)
    c = pd.read_csv(O / 'Tool-outputs/Standard-rate-by-export-capability.csv').set_index('export_class')
    ax2 = fig.add_axes([.82, .1, .16, .73])
    order = [('native', 'Default'), ('built-in', 'Option'), ('none', 'None')]
    for i, (k, l) in enumerate(order):
        v = c.standard_pct[k]
        ax2.bar(i, v, color=col[k], width=.7)
        ax2.errorbar(i, v, yerr=[[v - c.ci_low[k]], [c.ci_high[k] - v]], color=INK, lw=.7, capsize=2)
        ax2.text(i, c.ci_high[k] + 2, f'{v:.0f}%', ha='center', fontsize=8, fontweight='bold', color=col[k])
    ax2.set_xticks(range(3)); ax2.set_xticklabels([l for _, l in order]); ax2.set_ylim(0, 75); yfmt(ax2)
    ax2.set_title('By export setting', loc='left', fontsize=7)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=col[k], label=lab[k]) for k in lab], loc='lower right', title='Standard-format export', title_fontsize=6.5)
    panel(ax, 'a', x=-.3); panel(ax2, 'b', x=-.3, y=1.0)
    save(fig, '06-standards-by-software')


# 7. Software hometowns
def fig_hometowns():
    h = pd.read_csv(O / 'Software-geography/Home-country-lift.csv')
    h = h.copy()
    nm = {'maxquant': 'MaxQuant', 'spectronaut': 'Spectronaut', 'diann': 'DIA-NN', 'msfragger': 'MSFragger', 'peaks': 'PEAKS', 'mascot': 'Mascot',
          'scaffold': 'Scaffold', 'skyline': 'Skyline', 'pfind': 'pFind', 'openms': 'OpenMS', 'comet': 'Comet', 'msgf': 'MS-GF+', 'byonic': 'Byonic',
          'metamorpheus': 'MetaMorpheus', 'xtandem': 'X!Tandem', 'progenesis': 'Progenesis', 'peptideshaker': 'PeptideShaker'}
    h['label'] = h.tool.map(nm).fillna(h.tool) + ' · ' + h.country
    h = h.sort_values('lift')
    fig = plt.figure(figsize=(W, 84 * MM))
    over = int((h.lift > 1).sum()); sig = int((h.lift_ci_low > 1).sum())
    headline(fig, f'Software is named more often in deposits from its developers’ country: {over} of {len(h)} tool\u2013country pairs',
             'Use of each tool by groups in its developer’s country, relative to its use everywhere (1 = no home advantage); 95% group-bootstrap CI. Median 1.8\u00d7; permutation p = 0.0002.')
    ax = fig.add_axes([.25, .09, .72, .75])
    y = np.arange(len(h))
    ax.axvline(1, color=MUTED, lw=.7, ls='--')
    ax.hlines(y, h.lift_ci_low, h.lift_ci_high, color=HAIR, lw=1.6)
    ax.scatter(h.lift, y, color=np.where(h.lift > 1, ACC, ACC2), s=22, zorder=3)
    for i, v in enumerate(h.lift):
        ax.text(max(v, h.lift_ci_high.iloc[i]) + .08, i, f'{v:.1f}×', va='center', fontsize=6.2)
    ax.set_yticks(y); ax.set_yticklabels(h.label, fontsize=6.4); ax.set_xscale('log')
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: '')); ax.set_xticks([.5, 1, 2, 4, 8, 16]); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}×'))
    ax.tick_params(axis='y', length=0); ax.set_xlim(.4, 24)
    save(fig, '07-hometowns')


# 8. Geography: who deposits
def fig_countries():
    f = F[F.release_year.between(2013, 2025) & F.country.notna()]
    p = pd.crosstab(f.release_year, f.country, normalize='index') * 100
    top = ['United States', 'Germany', 'China', 'United Kingdom', 'France']
    fig = plt.figure(figsize=(W, 66 * MM))
    headline(fig, 'China’s share of public PRIDE datasets grew fivefold',
             'Share of country-resolved PRIDE datasets by the contributing group’s country, by release year (28,225 datasets).')
    ax = fig.add_axes([.06, .13, .74, .67])
    cols = {'United States': '#2e6f9e', 'Germany': '#9aa5ad', 'China': ACC, 'United Kingdom': '#8d6cab', 'France': ACC3}
    pos = repel({k: p[k][2025] for k in top}, 1.6)
    for k in top:
        y = p[k].reindex(YEARS); ax.plot(YEARS, y, color=cols[k], lw=2.2 if k == 'China' else 1.3)
        ax.text(2025.2, pos[k], f'{k}  {y[2025]:.0f}%', color=cols[k], fontsize=6.8, va='center', fontweight='bold' if k == 'China' else 'normal')
    ax.text(2013, p['China'][2013] + 1.2, f"{p['China'][2013]:.0f}%", color=ACC, fontsize=7.5, fontweight='bold')
    ax.set_xlim(2013, 2025); ax.set_ylim(0, 25); yfmt(ax)
    save(fig, '08-countries')


# 9. Clocks: what you download is already old
def fig_clocks():
    d = pd.read_csv(A / 'H11-clocks/Dataset-clocks.csv.gz', usecols=['acq_to_rel', 'db_to_rel', 'sw_to_rel'])
    rows = [('acq_to_rel', 'Data (latest acquisition)', ACC2), ('db_to_rel', 'Sequence database', '#5a9e6f'), ('sw_to_rel', 'Software release (lower bound)', ACC)]
    fig = plt.figure(figsize=(W, 66 * MM))
    headline(fig, 'Dated files trace the path to public release: median acquisition age 2.2 years, software-release age at least 3.4',
             'Age at public release of the latest raw-file acquisition date (datasets with dated raw-file names, 34% of those with raw files), the sequence database date and the software release (first seen in PRIDE; lower bound). The plot excludes one negative database interval and two negative software intervals; ages above 10 years are pooled in the final bin.')
    ax = fig.add_axes([.21, .13, .76, .67])
    bins = np.linspace(0, 10, 81)
    for i, (c, lab, col) in enumerate(rows):
        v = d[c].dropna() / 365.25; v = v[v >= 0]
        h, e = np.histogram(v.clip(upper=10), bins=bins, density=True)
        h = h / h.max() * .85
        ax.fill_between(e[:-1], i, i + h, step='post', color=col, alpha=.75, lw=0)
        md = v.median(); ax.plot([md, md], [i, i + .9], color=INK, lw=.9)
        ax.text(md + .12, i + .62, f'median {md:.1f} y', fontsize=7, fontweight='bold')
        ax.text(-.15, i + .3, f'{lab}\n(n = {len(v):,})', ha='right', va='center', fontsize=6.6)
    ax.set_yticks([]); ax.set_xlim(0, 10); ax.set_ylim(0, 3)
    ax.set_xlabel('Years before public release'); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}'))
    save(fig, '09-clocks')


# 10. Releases arrive in waves
def fig_release_waves():
    b = pd.read_csv(O / 'Attacks/release-timing/B1-bulk-days.csv', parse_dates=['release_date'])
    f = F.assign(d=pd.to_datetime(F.release_date).dt.normalize())
    f = f[f.d >= '2018-01-01']
    daily = f.groupby('d').size()
    bulk = set(b.release_date.dt.normalize())
    share = daily[daily.index.isin(bulk)].sum() / daily.sum() * 100
    fig = plt.figure(figsize=(W, 62 * MM))
    headline(fig, f'Public releases come in waves: {share:.0f}% of datasets since 2018 went public on {len(bulk)} days',
             'PRIDE datasets made public per day, 2018\u20132025. Red: days with \u226550 releases. Datasets released on those days went public a median 95 days after their paper (other days: 11).')
    ax = fig.add_axes([.06, .14, .92, .64])
    isb = daily.index.isin(bulk)
    ax.vlines(daily.index[~isb], 0, daily[~isb], color='#9fb4c6', lw=.5)
    ax.vlines(daily.index[isb], 0, daily[isb], color=ACC, lw=1.1)
    mx = daily.idxmax(); ax.text(mx, daily[mx] + 5, f'{daily[mx]} on {mx:%d %b %Y}', fontsize=6.8, color=ACC, ha='center', fontweight='bold')
    ax.set_ylabel('Datasets released per day'); ax.set_ylim(0, daily.max() * 1.15)
    save(fig, '10-release-waves')
    return share, len(bulk)


# 11. Long tail of bytes
def fig_bytes():
    s = pd.read_csv(A / 'H1-H3-scale-instruments/H1-dataset-scale.csv.gz', usecols=['bytes_all'])
    v = np.sort(s.bytes_all.fillna(0).to_numpy())[::-1]
    cum = np.cumsum(v) / v.sum() * 100; x = np.arange(1, len(v) + 1) / len(v) * 100
    tot = v.sum() / 1e15
    fig = plt.figure(figsize=(W, 66 * MM))
    t1 = cum[int(len(v) * .01) - 1]; t10 = cum[int(len(v) * .10) - 1]
    headline(fig, f'1% of datasets hold {t1:.0f}% of the archive’s {tot:.2f} petabytes',
             'Cumulative share of PRIDE bytes held by the largest datasets (all deposited files, 2013–2025).')
    ax = fig.add_axes([.07, .14, .55, .64])
    ax.plot(x, cum, color=ACC2, lw=2); ax.fill_between(x, 0, cum, color=ACC2, alpha=.1)
    for q, c in [(1, t1), (10, t10)]:
        ax.plot([q, q], [0, c], color=ACC, lw=.8, ls='--'); ax.plot([0, q], [c, c], color=ACC, lw=.8, ls='--')
        ax.text(q * 1.15, c - 7, f'top {q}% → {c:.0f}%', color=ACC, fontsize=7.5, fontweight='bold')
    ax.set_xscale('log'); ax.set_xlim(.01, 100); ax.set_ylim(0, 102)
    ax.set_xticks([.01, .1, 1, 10, 100]); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}%')); yfmt(ax)
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ''))
    ax.set_xlabel('Largest datasets (share of all datasets)'); ax.set_ylabel('Share of all bytes')
    ax2 = fig.add_axes([.72, .14, .26, .64])
    lay = pd.read_csv(A / 'H4-H7-databases-storage/H7-top-level-bytes-by-layer.csv')
    ax2.axis('off')
    ax2.text(0, .95, 'What the bytes are', fontsize=7.5, fontweight='bold', va='top')
    g = lay.groupby('layer').pct_bytes.sum().sort_values(ascending=False).head(5)
    names = {'vendor_raw': 'Vendor raw files', 'tool_results': 'Analysis results', 'open_spectra': 'Open spectra (mzML, MGF)',
             'archived_results': 'Zipped results', 'other': 'Other'}
    for i, (k, val) in enumerate(g.items()):
        ax2.barh(.7 - i * .16, val / 100, height=.06, color=ACC2 if i == 0 else '#9fb4c6', left=0)
        ax2.text(0, .7 - i * .16 + .045, f'{names.get(k, k)}  {val:.0f}%', fontsize=6.5, va='bottom')
    ax2.set_xlim(0, 1); ax2.set_ylim(0, 1)
    panel(ax, 'a'); panel(ax2, 'b', x=0, y=.98)
    save(fig, '11-bytes')
    return t1, t10, tot


# 12. Results are there when you open the files
def fig_hidden_results():
    b = pd.read_csv(O / 'Archive-contents/Before-after-by-year.csv').set_index('release_year')
    h = pd.read_csv(O / 'Header-sniff/Recognisable-by-year-three-ways.csv').set_index('release_year')
    yr = 2025
    steps = [('Standard format\n(mzIdentML/mzTab)', b.has_standard_results[yr]), ('Recognisable\nfrom file names', h.recognisable_before[yr]),
             ('+ files inside\nzip archives', h.recognisable_after[yr]), ('+ table headers\nread', h.recognisable_with_headers[yr])]
    fig = plt.figure(figsize=(W, 66 * MM))
    headline(fig, f'Reading deposited outputs reveals identifiable results in {steps[-1][1]:.0f}% of 2025 datasets; {steps[0][1]:.0f}% contain standard result files',
             f'PRIDE datasets released in {yr} (n = {int(b.datasets[yr]):,}) with identification/quantification results, by how deeply the files are read.')
    ax = fig.add_axes([.06, .14, .55, .64])
    prev = 0
    for i, (lab, v) in enumerate(steps):
        ax.bar(i, v - prev, bottom=prev, color=ACC if i == 0 else ACC2, alpha=1 if i in (0, 3) else .6, width=.65)
        ax.text(i, v + 2, f'{v:.0f}%', ha='center', fontsize=8.5, fontweight='bold', color=ACC if i == 0 else ACC2)
        prev = v
    ax.set_xticks(range(4)); ax.set_xticklabels([s for s, _ in steps], fontsize=6.4); ax.set_ylim(0, 100); yfmt(ax)
    ax2 = fig.add_axes([.7, .14, .28, .64])
    yy = [y for y in YEARS if y in h.index and y >= 2015]
    ax2.plot(yy, b.has_standard_results[yy], color=ACC, lw=1.6, label='Standard format')
    ax2.plot(yy, h.recognisable_with_headers[yy], color=ACC2, lw=1.6, label='Identifiable, all files read')
    ax2.set_ylim(0, 100); yfmt(ax2); ax2.legend(loc='center left', fontsize=6); ax2.set_xticks([2015, 2020, 2025])
    ax2.set_title('By release year', loc='left', fontsize=7)
    panel(ax, 'a', x=-.1, y=1.04); panel(ax2, 'b', x=-.2, y=1.06)
    save(fig, '12-hidden-results')
    return steps


# 13. Documentation habits persist
def fig_habits():
    p = pd.read_csv(A / 'H5-H6-H8-strategies-software-docs/H8-documentation-persistence-within-groups.csv')
    nm = {'doc_readme': 'README', 'doc_sdrf': 'SDRF sample table', 'doc_design_author': 'Own design table', 'doc_params': 'Parameter/log files'}
    t = p[p.previous_dataset_had_it].set_index('documentation')
    fig = plt.figure(figsize=(W, 62 * MM))
    headline(fig, 'README and SDRF files recur within contributor groups far above their base rates',
             'Share of a group’s next PRIDE dataset carrying each file when the previous one did, vs the rate expected for datasets of the same years.')
    ax = fig.add_axes([.2, .13, .78, .68])
    y = np.arange(len(nm))[::-1]
    for yi, k in zip(y, nm):
        ax.barh(yi + .18, t.pct_now[k], height=.34, color=ACC)
        ax.barh(yi - .18, t.pct_expected_same_years[k], height=.34, color='#c9cdd1')
        ax.text(t.pct_now[k] + 1, yi + .18, f'{t.pct_now[k]:.0f}%', va='center', fontsize=7, fontweight='bold', color=ACC)
        ax.text(t.pct_expected_same_years[k] + 1, yi - .18, f'{t.pct_expected_same_years[k]:.1f}%', va='center', fontsize=6.3, color=MUTED)
        ax.text(78, yi, f'{t.pct_now[k] / t.pct_expected_same_years[k]:.0f}×', va='center', fontsize=10, fontweight='bold', color=INK)
    ax.set_yticks(y); ax.set_yticklabels(list(nm.values())); ax.set_xlim(0, 85); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%'))
    ax.tick_params(axis='y', length=0)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=ACC, label='Previous dataset had it'), Patch(color='#c9cdd1', label='Expected for the same years')],
              loc='lower center', bbox_to_anchor=(.5, .98), ncol=2)
    save(fig, '13-habits')


# 26. Later literature: which datasets are named
def fig_reuse():
    import json
    t = pd.read_csv(O / 'Reuse-mentions/Reuse-rates.csv')
    t = t[t.outcome == 'D1']
    orr = json.load(open(O / 'Attacks/reuse-mentions/Results.json'))['A3_adjusted_odds_ratios']
    allr = t[t.feature == 'All'].iloc[0]
    rows = [('SDRF sample table', 'SDRF', 'sdrf'), ('DIA acquisition', 'DIA', 'dia'), ('Human samples', 'Human', 'human'),
            ('Thermo Fisher only', 'Thermo raw', 'thermo'), ('Standard result files', 'Standard files', 'standard')]
    fig = plt.figure(figsize=(W, 72 * MM))
    headline(fig, 'Qualifying papers more often name DIA and human datasets; standard result files show no detectable association',
             f'PRIDE datasets released 2013\u20132022 (n = {int(allr.n):,}) with a qualifying non-primary accession mention during the release calendar year '
             f'or the next three calendar years ({allr.std_pct:.1f}% overall). Qualifying papers name datasets from \u22652 contributor groups and no more than 100 datasets; '
             'publication and release order within a year is unresolved. Thermo Fisher only means a single-vendor assignment from instrument labels. '
             'Rates are standardised to release year and human status; odds ratios use the adjusted model in Methods. '
             'Europe PMC text mining provides an incomplete proxy that does not establish reanalysis.')
    ax = fig.add_axes([.17, .19, .42, .58]); ax2 = fig.add_axes([.71, .19, .27, .58])
    y = np.arange(len(rows))[::-1]
    for yi, (lab, feat, k) in zip(y, rows):
        a = t[(t.feature == feat) & (t.level == 'True')].iloc[0]; b = t[(t.feature == feat) & (t.level == 'False')].iloc[0]
        ax.plot([b.std_pct, a.std_pct], [yi, yi], color=HAIR, lw=2, zorder=1)
        ax.errorbar(b.std_pct, yi, xerr=[[b.std_pct - b.lo], [b.hi - b.std_pct]], fmt='o', color='#9aa5ad', ms=4, lw=.8, capsize=0)
        c = ACC if k != 'standard' else ACC2
        mfc = 'white' if k == 'sdrf' else c
        ax.errorbar(a.std_pct, yi, xerr=[[a.std_pct - a.lo], [a.hi - a.std_pct]], fmt='o', color=c, mfc=mfc, ms=5, lw=.9, capsize=0)
        ax.text(a.hi + .4, yi + .18, f'{a.std_pct:.1f}%', fontsize=6.5, color=c, fontweight='bold')
        e, lo, hi = orr[k]
        ax2.errorbar(e, yi, xerr=[[e - lo], [hi - e]], fmt='s', color=c, mfc=mfc, ms=4.5, lw=.9, capsize=0)
        ax2.text(hi * 1.08, yi, f'{e:.2f}', fontsize=6.5, va='center')
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows]); ax.tick_params(axis='y', length=0)
    ax.set_xscale('log'); ax.set_xlim(1, 40); ax.set_xticks([1, 2, 5, 10, 20])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}%')); ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ''))
    ax.set_xlabel('Share named (standardised to year \u00d7 human status)')
    ax.text(1.05, y[0] + .55, 'grey: without the feature', fontsize=6, color=MUTED)
    ax.set_ylim(-.6, len(rows) - .2)
    ax2.axvline(1, color=MUTED, lw=.7, ls='--'); ax2.set_xscale('log'); ax2.set_xlim(.5, 12)
    ax2.set_xticks([.5, 1, 2, 4, 8]); ax2.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}')); ax2.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ''))
    ax2.set_yticks(y); ax2.set_yticklabels([]); ax2.tick_params(axis='y', length=0); ax2.set_ylim(-.6, len(rows) - .2)
    ax2.set_xlabel('Adjusted odds ratio (95% CI)')
    ax2.set_title('Adjusted for year, field, vendor, software,\ncitations of the original paper, size', loc='left', fontsize=6.3, color=MUTED)
    ax.set_title('Qualifying mentions by dataset feature', loc='left', fontsize=7)
    panel(ax, 'a', x=-.36, y=1.03); panel(ax2, 'b', x=-.08, y=1.03)
    footnote(fig, .17, -.005, 'Open circle: SDRF files were often added after release by reuse projects (91 of 205 predate the format), so direction cannot be read.',
             fontsize=6, color=MUTED)
    save(fig, '26-reuse')


if __name__ == '__main__':
    fig_instruments(); fig_run_size(); fig_scale(); s = fig_newcomers(); fig_software(); fig_standards(); fig_hometowns()
    fig_countries(); fig_clocks(); w = fig_release_waves(); b = fig_bytes(); h = fig_hidden_results(); fig_habits(); fig_reuse()
    print('newcomer share 2025', s[-1], 'waves', w, 'bytes', b, 'hidden', h)
