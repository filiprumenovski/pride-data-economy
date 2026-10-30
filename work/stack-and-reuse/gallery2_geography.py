"""Gallery figures 20-22: geography of deposits, software and practices (descriptive).

20-world-growth             proportional symbols at each country's label point (area = datasets released
                            2023-2025), coloured by the growth multiple versus 2016-2018 (binned around the
                            archive-wide multiple), with a Europe zoom and ranked dumbbells (2016-18 -> 2023-25
                            datasets, log scale) carrying 95% group-bootstrap intervals for each multiple.
21-tool-country-heatmap     tool x country lift (Software-geography/Tool-country-lift.csv, presence metric,
                            equal contributor-group weighting, 2019-2025, countries with >= 60 groups), log
                            colour centred on 1; square area = evidence (expected naming groups); rows and
                            columns ordered by average-linkage clustering of evidence-shrunk log lift; home
                            cells outlined; dots mark cells whose interval excludes 1.
22-country-small-multiples  DIA share and standard-format (mzIdentML/mzTab, corrected nongenerated_expanded)
                            share by release year, 2015-2025, for the 20 countries with most datasets;
                            3-year count-weighted window; 95% group-bootstrap bands; archive line in each panel.

Definitions: `standards` is the corrected nongenerated_expanded flag (the legacy inclusive flag is never
used). `dia_expanded` False means "DIA not stated". Country = the contributing group's country; 22% of
datasets have no resolved country and are excluded from country rates (denominators recorded in the CSVs).
A deposited instrument or tool is not an ownership record; first deposit is not adoption.
Uncertainty: 95% percentile intervals resampling contributor groups (pi_group; deposits with no group are
their own cluster), seed 20261001. Claims are descriptive only.
"""
import sys
import numpy as np, pandas as pd
import geopandas as gpd
from pyproj import Transformer
from matplotlib.patches import Rectangle, Circle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3, panel
from scipy.cluster.hierarchy import linkage, leaves_list, optimal_leaf_ordering, dendrogram
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from common import features
from common import O, SEED

D = G / 'data'; D.mkdir(exist_ok=True)
HERE = __import__('pathlib').Path(__file__).resolve().parent
LAND = '#eeece7'
MINBASE = 15   # minimum 2016-2018 datasets for a coloured growth class
B = 2000


def gid(d):
    return d.pi_group.fillna('solo:' + d.accession)


def boot_ratio(d, num_mask, den_mask, rng, reps=B):
    """Group bootstrap of sum(num)/sum(den) where num/den are boolean deposit masks."""
    z = pd.DataFrame({'g': gid(d), 'a': num_mask.astype(int), 'b': den_mask.astype(int)}).groupby('g')[['a', 'b']].sum().to_numpy()
    idx = rng.integers(0, len(z), (reps, len(z)))
    a = z[:, 0][idx].sum(1); b = z[:, 1][idx].sum(1)
    with np.errstate(divide='ignore', invalid='ignore'):
        r = a / b
    return np.nanquantile(r, .025), np.nanquantile(r, .975)


def world_frame():
    w = gpd.read_file(HERE / 'ne_50m_admin_0_countries.geojson')
    return w[w.CONTINENT != 'Antarctica'].copy()


def match_names(w, index):
    names = {}
    for i, row in w.iterrows():
        for col in ('NAME', 'ADMIN', 'NAME_LONG'):
            if row[col] in index:
                names[i] = row[col]
    return names


# ============================================================================ 20 world growth
def fig_world_growth():
    rng = np.random.default_rng(SEED)
    f = F[F.release_year.between(2013, 2025)].copy()
    f['per'] = np.select([f.release_year.between(2016, 2018), f.release_year.between(2023, 2025)], ['early', 'late'], 'other')
    tot = f[f.per != 'other'].groupby('per').size()
    res = f[(f.per != 'other') & f.country.notna()]
    rtot = res.groupby('per').size()
    t = res.groupby(['country', 'per']).size().unstack(fill_value=0).rename(columns={'early': 'd_2016_2018', 'late': 'd_2023_2025'})
    t['added'] = t.d_2023_2025 - t.d_2016_2018
    t['multiple'] = t.d_2023_2025 / t.d_2016_2018.replace(0, np.nan)
    m0 = rtot['late'] / rtot['early']
    # group-bootstrap interval for each multiple (groups pooled across both periods)
    lo, hi = {}, {}
    for c, x in res.groupby('country'):
        if t.loc[c, 'd_2016_2018'] >= 5:
            lo[c], hi[c] = boot_ratio(x, x.per == 'late', x.per == 'early', rng)
    t['multiple_lo'] = pd.Series(lo); t['multiple_hi'] = pd.Series(hi)
    if FROZEN_DATA:
        saved = pd.read_csv(D / '20-world-growth.csv', index_col='country')
        t[['multiple_lo', 'multiple_hi']] = saved[['multiple_lo', 'multiple_hi']]
    w = world_frame()
    nm = match_names(w, t.index)
    w['country'] = w.index.map(nm)
    assert set(t.index) <= set(nm.values()), set(t.index) - set(nm.values())
    meta = w.dropna(subset=['country']).drop_duplicates('country').set_index('country')
    t['continent'] = meta.CONTINENT.reindex(t.index)
    tr = Transformer.from_crs('EPSG:4326', 'EPSG:8857', always_xy=True)
    # label points; nudge a few so symbols sit on the territory people expect
    lx, ly = meta.LABEL_X.copy(), meta.LABEL_Y.copy()
    t['x'], t['y'] = tr.transform(lx.reindex(t.index).to_numpy(), ly.reindex(t.index).to_numpy())
    # regional decomposition of growth
    reg = t.groupby('continent')[['d_2016_2018', 'd_2023_2025', 'added']].sum()
    reg['share_of_added'] = reg.added / reg.added.sum()
    eu = reg.loc['Europe', 'share_of_added']; asia = reg.loc['Asia', 'share_of_added']; na = reg.loc['North America', 'share_of_added']
    cn = t.loc['China']; de = t.loc['Germany']; us = t.loc['United States']
    # colour classes around the archive multiple
    bins = [0, 2, 2.75, 4, 6, np.inf]
    cols = ['#2e6f9e', '#93b6d2', '#b3aea4', '#ec9d90', '#c23b4e']
    labs = ['under 2×', '2–2.75×', '2.75–4×', '4–6×', '6× or more']
    t['cls'] = np.where(t.d_2016_2018 >= MINBASE, np.digitize(t.multiple.fillna(99), bins[1:-1]), -1)
    K = 0.30   # marker area (pt^2) per dataset

    def colour(r):
        return cols[int(r.cls)] if r.cls >= 0 else 'white'

    DIRS = {'r': (1, 0, 'left', 'center'), 'l': (-1, 0, 'right', 'center'), 't': (0, 1, 'center', 'bottom'), 'b': (0, -1, 'center', 'top'),
            'tr': (.75, .75, 'left', 'bottom'), 'tl': (-.75, .75, 'right', 'bottom'), 'br': (.75, -.75, 'left', 'top'), 'bl': (-.75, -.75, 'right', 'top')}

    def symbols(ax, sub, labels, fs=6.2):
        """labels: country -> direction key (optionally with extra offset in points)."""
        sub = sub.sort_values('d_2023_2025', ascending=False)
        ax.scatter(sub.x, sub.y, s=sub.d_2023_2025 * K, c=[colour(r) for _, r in sub.iterrows()],
                   edgecolor=[('white' if r.cls >= 0 else '#8a867e') for _, r in sub.iterrows()],
                   lw=[.5 if r.cls >= 0 else .6 for _, r in sub.iterrows()], zorder=3)
        for c, spec in labels.items():
            if c not in sub.index:
                continue
            d, extra = (spec, 0) if isinstance(spec, str) else spec
            r = sub.loc[c]; ux, uy, ha, va = DIRS[d]; rad = np.sqrt(r.d_2023_2025 * K) / 2 + 1.2 + extra
            ax.annotate(c, (r.x, r.y), xytext=(ux * rad, uy * rad), textcoords='offset points', ha=ha, va=va, fontsize=fs, color=INK,
                        zorder=5, path_effects=[pe.withStroke(linewidth=1.8, foreground='white')])

    fig = plt.figure(figsize=(W, 126 * MM))
    headline(fig, f'Half of PRIDE\u2019s {m0:.1f}-fold growth came from Europe, yet China ({cn.multiple:.1f}\u00d7) added more datasets than Germany ({de.multiple:.1f}\u00d7)',
             f'Circle area: datasets released 2023\u20132025, by the contributing group\u2019s country; colour: multiple of the 2016\u20132018 count (archive-wide {m0:.1f}\u00d7). '
             f'{rtot["early"]:,} and {rtot["late"]:,} country-resolved datasets\n'
             f'({rtot["early"] / tot["early"]:.0%} of each period; the rest have no resolved country). Datasets added: Europe {reg.loc["Europe", "added"]:,.0f} ({eu:.0%}), '
             f'North America {reg.loc["North America", "added"]:,.0f} ({na:.0%}), Asia {reg.loc["Asia", "added"]:,.0f} ({asia:.0%}; China {cn.added:,.0f}). Deposits only: not ownership or national output.')
    wp = w.to_crs('EPSG:8857')
    # world map
    ax = fig.add_axes([.0, .355, .62, .5], anchor='N')
    wp.plot(ax=ax, color=LAND, edgecolor='white', linewidth=.2)
    ax.set_xlim(-1.40e7, 1.62e7); ax.set_ylim(-5.0e6, 7.6e6); ax.set_aspect('equal'); ax.axis('off')
    symbols(ax, t, {'United States': 'bl', 'China': 'b', 'Canada': 't', 'Australia': 'b', 'India': 'b', 'Brazil': 'b',
                    'Japan': 'r', 'South Korea': 'tl', 'Israel': 'l', 'Singapore': 'r', 'Taiwan': 'r'})
    eb = [-1.05e6, 4.25e6, 2.25e6, 7.75e6]  # Europe zoom window, EPSG:8857 metres
    ax.add_patch(Rectangle((eb[0], eb[1]), eb[2] - eb[0], eb[3] - eb[1], fill=False, lw=.5, ec='#8a867e', zorder=6))
    # Europe zoom (same symbol scale)
    ax2 = fig.add_axes([.625, .355, .37, .5], anchor='N')
    wp.plot(ax=ax2, color=LAND, edgecolor='white', linewidth=.4)
    ax2.set_xlim(eb[0], eb[2]); ax2.set_ylim(eb[1], eb[3]); ax2.set_aspect('equal')
    ax2.set_xticks([]); ax2.set_yticks([]); ax2.set_xlabel(''); ax2.set_ylabel('')
    for sp in ax2.spines.values():
        sp.set_visible(True); sp.set_linewidth(.5); sp.set_color('#8a867e')
    inwin = t[(t.x > eb[0]) & (t.x < eb[2]) & (t.y > eb[1]) & (t.y < eb[3])]
    symbols(ax2, inwin, {'Germany': 'tr', 'United Kingdom': 'tl', 'France': 'bl', 'Netherlands': 't', 'Switzerland': 'b',
                         'Spain': 'b', 'Italy': 'br', 'Denmark': 'tl', 'Sweden': 'r', 'Austria': 'br', 'Czechia': 'r', 'Belgium': 'l',
                         'Poland': 'r', 'Norway': 'l', 'Portugal': 'l', 'Greece': 'b', 'Finland': 'r', 'Ireland': 'l', 'Hungary': 'r'})
    panel(ax, 'a', x=.03, y=.97); panel(ax2, 'b', x=-.015, y=.955)
    ax2.text(.025, .975, 'Europe', transform=ax2.transAxes, fontsize=6.6, color=MUTED, va='top', fontweight='bold')
    # keys over the South Pacific
    kx = fig.add_axes([.035, .405, .2, .028])
    for i, c in enumerate(cols):
        kx.add_patch(Rectangle((i, 0), 1, 1, fc=c, ec='white', lw=.6))
    for i, v in enumerate(bins[1:-1]):
        kx.text(i + 1, -.3, f'{v:g}\u00d7', ha='center', va='top', fontsize=6)
    kx.add_patch(Rectangle((5.35, 0), 1, 1, fc='white', ec='#8a867e', lw=.6))
    kx.text(5.85, -.3, f'<{MINBASE} in\n2016\u201318', ha='center', va='top', fontsize=6, linespacing=.95)
    pos = 2 + (m0 - bins[2]) / (bins[3] - bins[2])
    kx.plot([pos, pos], [-.08, 1.08], color=INK, lw=.9, solid_capstyle='butt')
    kx.text(pos, 1.25, f'archive {m0:.1f}\u00d7', ha='center', va='bottom', fontsize=6, color=INK)
    kx.text(0, 2.5, 'Multiple, 2016\u201318 \u2192 2023\u201325', fontsize=6.4, color=INK, va='bottom', fontweight='bold')
    kx.set_xlim(0, 6.4); kx.set_ylim(0, 1); kx.axis('off')
    sx = fig.add_axes([.27, .392, .12, .07]); sx.axis('off'); sx.set_xlim(0, 1); sx.set_ylim(0, 1)
    for i, v in enumerate([100, 500, 2000]):
        sx.scatter(.1 + i * .34, .35, s=v * K, facecolor='none', edgecolor='#8a867e', lw=.6, clip_on=False)
        sx.text(.1 + i * .34, -.25, f'{v:,}', ha='center', va='top', fontsize=6)
    sx.text(0, 1.2, 'Datasets 2023\u201325', fontsize=6.4, fontweight='bold', va='bottom')
    # dumbbells: 20 largest contributors
    top = t.sort_values('d_2023_2025', ascending=False).head(20)
    fig_letter(fig, .008, .305, 'c', va='bottom')
    fig.text(.03, .305, 'The 20 largest contributors, 2023\u201325: datasets in each period (log scale) and growth multiple with 95% group-bootstrap interval',
             fontsize=7, fontweight='bold', color=INK)
    for k, part in enumerate([top.iloc[:10], top.iloc[10:]]):
        a = fig.add_axes([.115 + k * .5, .06, .29, .21])
        y = np.arange(len(part))[::-1]
        for yi, (c, r) in zip(y, part.iterrows()):
            col = colour(r)
            a.plot([r.d_2016_2018, r.d_2023_2025], [yi, yi], color=col if r.cls >= 0 else '#b9b5ad', lw=1.7, solid_capstyle='butt', zorder=2)
            a.scatter(r.d_2016_2018, yi, s=10, facecolor='white', edgecolor='#6d6a64', lw=.6, zorder=3)
            a.scatter(r.d_2023_2025, yi, s=18, color=col if r.cls >= 0 else 'white', edgecolor='white' if r.cls >= 0 else '#6d6a64', lw=.5, zorder=3)
            a.text(1.035, yi, f'{r.multiple:.1f}\u00d7', transform=a.get_yaxis_transform(), va='center', ha='left', fontsize=6.4, color=INK,
                   fontweight='bold' if c in ('China', 'Germany') else 'normal')
            a.text(1.2, yi, f'{r.multiple_lo:.1f}\u2013{r.multiple_hi:.1f}', transform=a.get_yaxis_transform(), va='center', ha='left', fontsize=6, color=MUTED)
        a.set_yticks(y); a.set_yticklabels(part.index, fontsize=6.4)
        a.set_xscale('log'); a.set_xlim(12, 3000); a.set_ylim(-.7, len(part) - .2)
        a.set_xticks([30, 100, 300, 1000, 3000]); a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:,.0f}'))
        a.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ''))
        a.tick_params(axis='y', length=0)
        a.grid(axis='x', color=HAIR, lw=.4); a.set_axisbelow(True)
        a.set_xlabel('Datasets (log scale)', fontsize=6.4, labelpad=1.5)
        a.text(1.035, len(part) - .1, 'multiple', transform=a.get_yaxis_transform(), fontsize=6, color=MUTED, va='bottom')
        a.text(1.2, len(part) - .1, '95% CI', transform=a.get_yaxis_transform(), fontsize=6, color=MUTED, va='bottom')
        if k == 0:
            r = part.iloc[0]; yy = len(part) - 1
            a.annotate('2016\u201318', (r.d_2016_2018, yy), xytext=(0, 4.5), textcoords='offset points', ha='center', va='bottom', fontsize=6, color=MUTED)
            a.annotate('2023\u201325', (r.d_2023_2025, yy), xytext=(0, 4.5), textcoords='offset points', ha='center', va='bottom', fontsize=6, color=MUTED)
    save(fig, '20-world-growth')
    out = t.drop(columns=['x', 'y', 'cls']).sort_values('d_2023_2025', ascending=False)
    out['archive_multiple_resolved'] = m0
    out['denominator_resolved_2016_2018'] = rtot['early']; out['denominator_resolved_2023_2025'] = rtot['late']
    out['denominator_all_2016_2018'] = tot['early']; out['denominator_all_2023_2025'] = tot['late']
    out.round(4).to_csv(D / '20-world-growth.csv', index_label='country')
    reg.round(4).to_csv(D / '20-world-growth-continents.csv')
    print(reg); print(top[['d_2016_2018', 'd_2023_2025', 'multiple', 'multiple_lo', 'multiple_hi']].round(2))


# ============================================================================ 21 tool x country heatmap
TOOL_NAME = {'maxquant': 'MaxQuant', 'spectronaut': 'Spectronaut', 'diann': 'DIA-NN', 'msfragger': 'MSFragger', 'peaks': 'PEAKS',
             'mascot': 'Mascot', 'scaffold': 'Scaffold', 'skyline': 'Skyline', 'pfind': 'pFind', 'openms': 'OpenMS', 'comet': 'Comet',
             'msgf': 'MS-GF+', 'byonic': 'Byonic', 'metamorpheus': 'MetaMorpheus', 'xtandem': 'X!Tandem', 'progenesis': 'Progenesis',
             'peptideshaker': 'PeptideShaker'}
DIV = LinearSegmentedColormap.from_list('div', ['#1f5a85', '#5f95c0', '#b9d1e4', '#f6f5f2', '#f1bcb3', '#d9707a', '#a8283d'])
FULL_E = 5      # expected naming groups at which a cell is drawn full size
PRIOR_E = 2     # pseudo-groups for the evidence-shrunk lift used only to order rows and columns


def tool_coverage():
    """Per country: equal-group share of 2019-2025 deposits naming any of the 17 tools, same filter as the lift table."""
    d = features()
    d = d[d.release_year.between(2019, 2025) & d.country.notna() & d.pi_group.notna()].copy()
    ct = d.groupby('country').pi_group.nunique(); d = d[d.country.isin(ct[ct >= 60].index)].copy()
    d['any17'] = d[['sw_' + t for t in TOOL_NAME]].any(axis=1)
    g = d.groupby(['country', 'pi_group']).any17.mean()
    return g.groupby('country').mean(), g.mean()


def fig_tool_heatmap():
    L = pd.read_csv(O / 'Software-geography/Tool-country-lift.csv')
    L['expected_groups'] = L.groups * L.global_equal_group_share
    L['observed_group_equiv'] = L.equal_group_share * L.groups
    nz = L.groups_naming_tool > 0
    # interval rule: bootstrap interval for non-zero cells; for zero cells the percentile interval is degenerate,
    # so a shortfall counts only when >= 3 naming groups were expected (Poisson P(0) < 0.05)
    L['above'] = nz & (L.lift_ci_low > 1)
    L['below'] = (nz & (L.lift_ci_high < 1)) | (~nz & (L.expected_groups >= 3))
    L['shrunk_log2_lift'] = np.log2((L.observed_group_equiv + PRIOR_E) / (L.expected_groups + PRIOR_E))
    M = L.pivot(index='tool', columns='country', values='shrunk_log2_lift')
    Zr = optimal_leaf_ordering(linkage(M.values, 'average', metric='correlation'), M.values, metric='correlation')
    Zc = optimal_leaf_ordering(linkage(M.T.values, 'average', metric='correlation'), M.T.values, metric='correlation')
    rows = list(M.index[leaves_list(Zr)]); cols = list(M.columns[leaves_list(Zc)])
    cov, cov_all = tool_coverage()
    gl = L.groupby('tool').global_equal_group_share.first()
    n_groups = L.groupby('country')['groups'].first()
    N = len(L); sig = int(L.above.sum() + L.below.sum())
    h = L[L.home]; h_up = int(h.above.sum())
    o = L[~L.home]; o_up = int(o.above.sum())
    cn = L[L.country == 'China']
    fig = plt.figure(figsize=(W, 140 * MM))
    headline(fig, f'Software habits are national: {sig} of {N} tool\u2013country rates differ clearly from the archive, home cells most often',
             f'Lift = share of contributor groups (equal weight) naming a tool in a country \u00f7 its share across all groups; 2019\u20132025 deposits, 25 countries with \u226560 groups '
             f'({int(n_groups.sum()):,} groups).\nHome cells (outlined) lie above 1 in {h_up} of {len(h)} cases, other cells in {o_up} of {len(o)}. '
             f'At 95% about {0.05 * N:.0f} cells would differ by chance. Rows and columns are ordered by clustering, a display aid rather than natural groups.')
    x0, y0, cw, chh = .175, .2, .62 / len(cols), .5 / len(rows)
    ax = fig.add_axes([x0, y0, cw * len(cols), chh * len(rows)])
    norm = LogNorm(.25, 4)
    Li = L.set_index(['tool', 'country'])
    for i, t in enumerate(rows):
        for j, c in enumerate(cols):
            r = Li.loc[(t, c)]
            side = float(np.clip(np.sqrt(r.expected_groups / FULL_E), .38, 1))
            v = float(np.clip(r.lift, .25, 4))
            col = DIV(norm(v))
            ax.add_patch(Rectangle((j + (1 - side) / 2, i + (1 - side) / 2), side, side, fc=col, ec='none'))
            if r.home:
                ax.add_patch(Rectangle((j + .04, i + .04), .92, .92, fill=False, ec=INK, lw=.9, zorder=4))
            if r.above or r.below:
                lum = np.dot(col[:3], [.299, .587, .114])
                ax.scatter(j + .5, i + .5, s=2.6, color='white' if lum < .55 else INK, lw=0, zorder=5)
    ax.set_xlim(0, len(cols)); ax.set_ylim(len(rows), 0)
    ax.set_xticks(np.arange(len(cols)) + .5); ax.set_xticklabels(cols, rotation=90, fontsize=6.4)
    ax.set_yticks(np.arange(len(rows)) + .5); ax.set_yticklabels([TOOL_NAME[t] for t in rows], fontsize=6.6)
    ax.tick_params(length=0, pad=2)
    for sp in ax.spines.values():
        sp.set_linewidth(.5); sp.set_color('#8a867e')
    # column dendrogram + coverage bars
    axc = fig.add_axes([x0, y0 + chh * len(rows) + .085, cw * len(cols), .055])
    dc = dendrogram(Zc, no_plot=True)
    for xs, ys in zip(dc['icoord'], dc['dcoord']):
        axc.plot(np.array(xs) / 10, ys, color='#9b978f', lw=.5)
    axc.set_xlim(0, len(cols)); axc.set_ylim(0, None); axc.axis('off')
    axb = fig.add_axes([x0, y0 + chh * len(rows) + .012, cw * len(cols), .06])
    vals = cov.reindex(cols) * 100
    axb.bar(np.arange(len(cols)) + .5, vals, width=.72, color='#b3aea4', lw=0)
    axb.axhline(cov_all * 100, color=INK, lw=.6, ls=(0, (2, 1.5)))
    axb.set_xlim(0, len(cols)); axb.set_ylim(0, 100); axb.set_xticks([]); axb.set_yticks([0, 50, 100])
    axb.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%')); axb.tick_params(axis='y', labelsize=6, length=1.5, pad=1)
    axb.text(len(cols) + .3, cov_all * 100, f'all {cov_all:.0%}', fontsize=6, va='center', color=INK)
    axb.text(-.6, 108, 'Groups naming any of these 17 tools', fontsize=6.3, ha='right', va='bottom', color=INK, transform=axb.transData)
    # row dendrogram
    axr = fig.add_axes([.012, y0, .07, chh * len(rows)])
    dr = dendrogram(Zr, no_plot=True, orientation='left')
    for xs, ys in zip(dr['icoord'], dr['dcoord']):
        axr.plot(ys, np.array(xs) / 10, color='#9b978f', lw=.5)
    axr.set_ylim(len(rows), 0); axr.set_xlim(None, 0); axr.invert_xaxis(); axr.axis('off')
    # right margin: archive-wide share of groups naming each tool
    axg = fig.add_axes([x0 + cw * len(cols) + .012, y0, .085, chh * len(rows)])
    gv = gl.reindex(rows) * 100
    axg.barh(np.arange(len(rows)) + .5, gv, height=.62, color='#b3aea4', lw=0)
    for i, v in enumerate(gv):
        axg.text(v + 1, i + .5, f'{v:.1f}%' if v < 10 else f'{v:.0f}%', va='center', fontsize=6, color=INK)
    axg.set_ylim(len(rows), 0); axg.set_xlim(0, 52); axg.set_yticks([]); axg.set_xticks([])
    axg.set_title('Archive-wide\nshare of groups', fontsize=6.3, loc='left', pad=2, color=INK)
    for sp in axg.spines.values():
        sp.set_visible(False)
    # China note
    jc = cols.index('China')
    nb = int(cn.below.sum())
    axb.annotate(f'China: below the archive for {nb} of 17 tools', (jc + .5, vals['China']), xytext=(0, 14), textcoords='offset points',
                 ha='center', fontsize=6.2, color=INK, arrowprops=dict(arrowstyle='-', lw=.5, color=INK))
    # legends
    lx = fig.add_axes([.06, .045, .2, .016])
    grad = np.logspace(np.log10(.25), np.log10(4), 256)
    lx.imshow(DIV(norm(grad))[None, :, :], aspect='auto', extent=(np.log2(.25), np.log2(4), 0, 1))
    lx.set_yticks([]); lx.set_xticks(np.log2([.25, .5, 1, 2, 4])); lx.set_xticklabels(['\u2264\u00bc\u00d7', '\u00bd\u00d7', '1', '2\u00d7', '\u22654\u00d7'], fontsize=6)
    lx.tick_params(length=1.5, pad=1.5)
    lx.set_title('Lift (log scale)', fontsize=6.4, loc='left', pad=2)
    k2 = fig.add_axes([.32, .03, .4, .045]); k2.axis('off'); k2.set_xlim(0, 40); k2.set_ylim(0, 4.5)
    for i, e in enumerate([.5, 2, 5]):
        side = float(np.clip(np.sqrt(e / FULL_E), .38, 1)) * 1.6
        k2.add_patch(Rectangle((i * 3.2 + (1.6 - side) / 2, .9 + (1.6 - side) / 2), side, side, fc='#b3aea4', ec='none'))
        k2.text(i * 3.2 + .8, .5, f'{e:g}', ha='center', va='top', fontsize=6)
    k2.text(0, 3.4, 'Expected naming groups (evidence)', fontsize=6.4, va='bottom')
    k2.add_patch(Rectangle((18, .9), 1.6, 1.6, fill=False, ec=INK, lw=.9)); k2.text(20.2, 1.7, 'developer\u2019s home country', va='center', fontsize=6.2)
    k2.add_patch(Rectangle((18, -1.4), 1.6, 1.6, fc='#b9d1e4', ec='none')); k2.scatter(18.8, -.6, s=2.6, color=INK, lw=0)
    k2.text(20.2, -.6, '95% interval excludes 1', va='center', fontsize=6.2)
    save(fig, '21-tool-country-heatmap')
    out = L.assign(tool_name=L.tool.map(TOOL_NAME), row_order=L.tool.map({t: i for i, t in enumerate(rows)}),
                   col_order=L.country.map({c: i for i, c in enumerate(cols)}), country_any17_group_share=L.country.map(cov))
    out.sort_values(['row_order', 'col_order']).round(4).to_csv(D / '21-tool-country-heatmap.csv', index=False)
    print(sig, h_up, len(h), o_up, len(o), nb, rows, cols)


if __name__ == '__main__':
    which = sys.argv[1:] or ['20', '21', '22']
    if '20' in which: fig_world_growth()
    if '21' in which: fig_tool_heatmap()
    if '22' in which: fig_small_multiples()
