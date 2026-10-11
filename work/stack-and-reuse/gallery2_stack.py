"""Gallery figures 17-19: the analysis stack, the tool ecosystem and lab habits (descriptive).

17-stack-alluvial  datasets released 2023-2025: instrument vendor -> acquisition -> main software
                   family -> results deposited. Ribbons coloured by the deposited-results outcome.
18-tool-chord      ordered co-occurrence matrix of analysis tools named in the same dataset,
                   2019-2025; colour = lift over independence, area = co-naming datasets,
                   hollow = 95% group-bootstrap interval for lift includes 1.
19-lab-barcodes    the 200 contributor groups (hashed pi_group) with most datasets 2013-2025: one tick
                   per dataset at its release date, coloured by main software family; marginal =
                   per-group share of consecutive datasets that switch family. Persistence compared
                   with a within-release-year permutation of family labels.

Definitions: `standards` is the corrected nongenerated_expanded flag (top-level mzIdentML/mzTab,
/generated/ paths excluded); the legacy inclusive flag is never used. Software = named in the
structured software field or the data-processing protocol (absence = not detected). `dia_expanded`
False means "DIA not stated", not verified DDA. Instruments in a deposit are not ownership records.
Uncertainty: 95% percentile intervals resampling contributor groups (pi_group; deposits with no
group are their own cluster), seed 20261001. Claims are descriptive only.
"""
import numpy as np, pandas as pd
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle
import matplotlib.patheffects as pe
from galstyle import *
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3, VEND
from common import features, O, SEED, cluster_boot, SOFTWARE

D = G / 'data'; D.mkdir(exist_ok=True)
FAM_COL = {'MaxQuant': '#2e6f9e', 'Proteome Discoverer': '#8d6cab', 'Mascot': '#7f8c96', 'FragPipe': '#5a9e6f',
           'DIA-NN': '#d1495b', 'Spectronaut': '#e8a33d', 'Several of these': '#3d3d3d',
           'Other tools only': '#b7ad99', 'No tool named': '#e2dfd8'}
HALO = [pe.withStroke(linewidth=2.2, foreground='white')]


def headmm(fig, title, sub, h_mm):
    """galstyle.headline with offsets fixed in millimetres, so tall figures keep the house spacing."""
    if record_head(title, sub):
        return
    fig.text(.012, 1 - 1.6 / h_mm, title, fontsize=10.5, fontweight='bold', va='top', color=INK)
    fig.text(.012, 1 - 6.2 / h_mm, sub, fontsize=7, va='top', color=MUTED, linespacing=1.25)
    r = fig.canvas.get_renderer()
    for tx in fig.texts[-2:]:
        wmm = tx.get_window_extent(r).width / fig.dpi * 25.4
        if wmm > 181:
            print(f'WARNING headline text {wmm:.0f} mm wide: {tx.get_text()[:60]}')


def load():
    d = features()
    h = pd.read_csv(O / 'Habits/Dataset-form.csv.gz', usecols=['accession', 'form', 'recognisable'])
    d = d.merge(h, on='accession', how='left', validate='one_to_one')
    d['cluster'] = d.pi_group.fillna('acc:' + d.accession)
    return d


def ci(values, groups):
    lo, hi = cluster_boot(np.asarray(values, float), np.asarray(groups), reps=2000, seed=SEED)
    return 100 * lo, 100 * hi


# --------------------------------------------------------------------------------------------
# 17. Alluvial of the analysis stack
# --------------------------------------------------------------------------------------------
OUT_ORDER = ['Standard file', 'Tool-native results', 'No recognisable results']
OUT_COL = {'Standard file': '#1f5f8b', 'Tool-native results': '#c4c8cc', 'No recognisable results': '#e8a33d'}


def ribbon(ax, x0, x1, a0, a1, b0, b1, color, alpha=.82):
    """Ribbon from left edge [a0,a1] at x0 to right edge [b0,b1] at x1 (cubic Bezier sides)."""
    xm = (x0 + x1) / 2
    verts = [(x0, a0), (xm, a0), (xm, b0), (x1, b0), (x1, b1), (xm, b1), (xm, a1), (x0, a1), (x0, a0)]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4, MPath.LINETO, MPath.CURVE4, MPath.CURVE4,
             MPath.CURVE4, MPath.CLOSEPOLY]
    ax.add_patch(PathPatch(MPath(verts, codes), facecolor=color, edgecolor='none', alpha=alpha, lw=0))


def fig_alluvial(d):
    s = d[d.release_year.between(2023, 2025)].copy()
    N = len(s)
    vmap = {'Thermo Fisher': 'Thermo Fisher', 'Bruker': 'Bruker', 'SCIEX': 'SCIEX'}
    s['c_vendor'] = s.vendor_single.map(vmap).fillna('Other or several')
    s['c_acq'] = np.where(s.dia_expanded, 'DIA', 'DDA or not stated')
    s['c_sw'] = s.family
    # results deposited: corrected top-level standard flag first; everything else from the frozen form labels
    s['c_out'] = np.where(s.standards, 'Standard file',
                          np.where(s.form.isin(['tool output, loose', 'tool output, zipped only', 'tool output, renamed', 'standard']),
                                   'Tool-native results', 'No recognisable results'))
    n_std_in_archive = int(((s.form == 'standard') & ~s.standards).sum())
    n_adhoc = int((s.form == 'ad hoc tables only').sum()); n_none = int((s.form == 'no results files').sum())
    cols = ['c_vendor', 'c_acq', 'c_sw', 'c_out']
    order = {'c_vendor': ['Thermo Fisher', 'SCIEX', 'Other or several', 'Bruker'],
             'c_acq': ['DDA or not stated', 'DIA'],
             'c_sw': ['Mascot', 'Other tools only', 'Proteome Discoverer', 'Several of these', 'MaxQuant', 'No tool named',
                      'FragPipe', 'Spectronaut', 'DIA-NN'],
             'c_out': OUT_ORDER}
    head = {'c_vendor': 'Instrument vendor', 'c_acq': 'Acquisition', 'c_sw': 'Main software named', 'c_out': 'Results deposited'}
    gap = .022 * N
    H = N + gap * (len(order['c_sw']) - 1)
    xs = {'c_vendor': 0, 'c_acq': 1.05, 'c_sw': 2.1, 'c_out': 3.15}; nw = .045
    # node extents (y increases downward: we invert the axis)
    nodes = {}
    for c in cols:
        cnt = s[c].value_counts().reindex(order[c]).fillna(0)
        tot = cnt.sum() + gap * (len(cnt) - 1); y = (H - tot) / 2
        for k, v in cnt.items():
            nodes[(c, k)] = [y, y + v]; y += v + gap

    HM = 128; fig = plt.figure(figsize=(W, HM * MM))
    # key numbers (computed)
    dt = s.c_sw.isin(['DIA-NN', 'Spectronaut'])
    r_dt, r_ot = 100 * s.standards[dt].mean(), 100 * s.standards[~dt].mean()
    ci_dt, ci_ot = ci(s.standards[dt], s.cluster[dt]), ci(s.standards[~dt], s.cluster[~dt])
    k_dt = int(s.standards[dt].sum())
    headmm(fig, f'DIA-NN/Spectronaut datasets almost never include a standard results file: {r_dt:.1f}% vs {r_ot:.0f}%',
             f'PRIDE datasets released 2023–2025 (n = {N:,}). Width = datasets; colour = results deposited; std = share with a standard file.\n'
             f'Standard file = top-level mzIdentML/mzTab (corrected definition). DIA-NN or Spectronaut as only main tool: {k_dt} of {int(dt.sum()):,} '
             f'(95% CI {ci_dt[0]:.1f}–{ci_dt[1]:.1f}%);\nall others {ci_ot[0]:.0f}–{ci_ot[1]:.0f}% (group bootstrap). '
             f'DIA = stated in metadata or protocol. Vendor = instrument named in the deposit, not ownership.', HM)
    ax = fig.add_axes([.0, .0, 1, .84])
    # ribbons: for each adjacent column pair, split by final outcome so the outcome colour runs through
    for c0, c1 in zip(cols[:-1], cols[1:]):
        g = s.groupby(list(dict.fromkeys([c0, c1, 'c_out']))).size().rename('n').reset_index()
        g['o0'] = g[c0].map({k: i for i, k in enumerate(order[c0])}); g['o1'] = g[c1].map({k: i for i, k in enumerate(order[c1])})
        g['oo'] = g.c_out.map({k: i for i, k in enumerate(OUT_ORDER)})
        out_y = {k: nodes[(c0, k)][0] for k in order[c0]}; in_y = {k: nodes[(c1, k)][0] for k in order[c1]}
        g = g.sort_values(['o0', 'o1', 'oo']); top0 = {}
        for i, r in g.iterrows():
            top0[i] = out_y[r[c0]]; out_y[r[c0]] += r.n
        g = g.sort_values(['o1', 'o0', 'oo']); top1 = {}
        for i, r in g.iterrows():
            top1[i] = in_y[r[c1]]; in_y[r[c1]] += r.n
        for i, r in g.iterrows():
            ribbon(ax, xs[c0] + nw, xs[c1], top0[i], top0[i] + r.n, top1[i], top1[i] + r.n, OUT_COL[r.c_out],
                   alpha=.9 if r.c_out == 'Standard file' else .72)
    # nodes + labels
    std_by = {c: s.groupby(c).standards.mean() * 100 for c in cols}
    for (c, k), (y0, y1) in nodes.items():
        fc = OUT_COL[k] if c == 'c_out' else '#2b2b2b'
        ax.add_patch(Rectangle((xs[c], y0), nw, y1 - y0, facecolor=fc, edgecolor='white', lw=.4, zorder=3))
        n = y1 - y0; pct = 100 * n / N; ym = (y0 + y1) / 2
        lab = f'{k}'; num = f'{n:,.0f} · {pct:.0f}%' if pct >= 1 else f'{n:,.0f} · {pct:.1f}%'
        if c == 'c_vendor':
            sr = std_by[c][k]
            ax.text(xs[c] - .03, ym, lab, ha='right', va='bottom', fontsize=7, fontweight='bold', color=INK)
            t2 = ax.text(xs[c] - .03, ym, f'{sr:.0f}% std', ha='right', va='top', fontsize=6.3, color=OUT_COL['Standard file'],
                         fontweight='bold')
            ax.annotate(f'{num} \u00b7 ', xy=(0, 0), xycoords=t2, ha='right', va='bottom', fontsize=6.3, color=MUTED)
        elif c == 'c_out':
            extra = {'Tool-native results': f'incl. {n_std_in_archive} with mzIdentML/mzTab\nonly inside archives',
                     'No recognisable results': f'ad hoc tables {n_adhoc:,}; no results files {n_none:,}'}.get(k, 'mzIdentML or mzTab')
            ax.text(xs[c] + nw + .03, ym, f'{lab}\n\n', ha='left', va='center', fontsize=7, fontweight='bold',
                    color=OUT_COL[k] if k != 'Tool-native results' else '#6b7177', linespacing=1.05)
            ax.text(xs[c] + nw + .03, ym, f'\n{num}\n', ha='left', va='center', fontsize=6.3, color=INK, linespacing=1.05)
            ax.text(xs[c] + nw + .03, ym + .0 * N, f'\n\n\n{extra}', ha='left', va='center', fontsize=6, color=MUTED, linespacing=1.05)
        else:
            sr = std_by[c][k]
            ax.text(xs[c] + nw + .025, ym, lab, ha='left', va='bottom', fontsize=6.6, fontweight='bold', color=INK,
                    path_effects=HALO, zorder=5)
            t2 = ax.text(xs[c] + nw + .025, ym, f'{num} · ', ha='left', va='top', fontsize=6, color=INK, path_effects=HALO, zorder=5)
            # standard-file rate of the node, in the standard colour
            ax.annotate(f'{sr:.1f}% std' if sr < 1 else f'{sr:.0f}% std', xy=(1, 0), xycoords=t2, va='bottom', ha='left',
                        fontsize=6, color=OUT_COL['Standard file'], fontweight='bold', path_effects=HALO, zorder=5)
    # vendor case: how much of Bruker's lower rate follows its software/acquisition mix (Thermo rates applied)
    th, br = s[s.vendor_single == 'Thermo Fisher'], s[s.vendor_single == 'Bruker']
    rate = th.groupby(['family', 'dia_expanded']).standards.mean()
    br_exp = 100 * np.nanmean(br.set_index(['family', 'dia_expanded']).index.map(rate))
    br_obs, th_obs = 100 * br.standards.mean(), 100 * th.standards.mean()
    share_mix = (th_obs - br_exp) / (th_obs - br_obs)
    yb = nodes[('c_vendor', 'Bruker')][1]
    oth = s[s.c_vendor == 'Other or several'].vendor_single.value_counts()
    ax.text(-.55, yb + .03 * N, 'Other or several: ' + ', '.join(f'{k} {v}' for k, v in oth.items() if v >= 50) + ', other/unknown '
            + str(int(oth[oth < 50].sum())) + '.\n\n' + f'Bruker vs Thermo Fisher: {br_obs:.1f}% vs {th_obs:.1f}%.\nWith Thermo Fisher rates per software \u00d7 acquisition,\n'
            f'Bruker\u2019s mix would give {br_exp:.1f}%: about {share_mix:.0%} of the gap\nfollows software and acquisition.',
            ha='left', va='top', fontsize=6, color=MUTED, linespacing=1.2)
    for c in cols:
        ax.text(xs[c] + nw / 2, -0.035 * N, head[c], ha='center', va='bottom', fontsize=7, color=MUTED, fontweight='bold')
    ax.set_xlim(-.58, 3.92); ax.set_ylim(max(H, yb + .17 * N), -.05 * N); ax.axis('off')
    save(fig, '17-stack-alluvial')
    # plotted table + node table
    flows = s.groupby(cols).size().rename('datasets').reset_index()
    flows['pct_of_all'] = 100 * flows.datasets / N
    flows.columns = ['vendor', 'acquisition', 'software_family', 'results_deposited', 'datasets', 'pct_of_all']
    flows.to_csv(D / '17-stack-alluvial.csv', index=False)
    nd = pd.DataFrame([(c[2:], k, y1 - y0, 100 * (y1 - y0) / N, std_by[c].get(k, np.nan)) for (c, k), (y0, y1) in nodes.items()],
                      columns=['axis', 'node', 'datasets', 'pct_of_all', 'standard_file_pct'])
    nd.to_csv(D / '17-stack-alluvial-nodes.csv', index=False)
    return dict(br_obs=br_obs, th_obs=th_obs, br_exp=br_exp, N=N, r_dt=r_dt, ci_dt=ci_dt, r_ot=r_ot, ci_ot=ci_ot, n_dt=int(dt.sum()), k_dt=k_dt)


# --------------------------------------------------------------------------------------------
# 18. Which tools are named together (ordered co-occurrence matrix)
# --------------------------------------------------------------------------------------------
TOOL_NAME = {'mascot': 'Mascot', 'sequest': 'SEQUEST', 'peaks': 'PEAKS', 'scaffold': 'Scaffold', 'progenesis': 'Progenesis',
             'msgf': 'MS-GF+', 'peptideshaker': 'PeptideShaker/SearchGUI', 'comet': 'Comet', 'xtandem': 'X!Tandem',
             'byonic': 'Byonic', 'proteinpilot': 'ProteinPilot', 'plgs': 'PLGS', 'skyline': 'Skyline',
             'msfragger': 'MSFragger/FragPipe', 'pfind': 'pFind/pLink', 'openms': 'OpenMS/OpenSWATH', 'percolator': 'Percolator',
             'metamorpheus': 'MetaMorpheus', 'maxquant': 'MaxQuant', 'pd': 'Proteome Discoverer', 'diann': 'DIA-NN',
             'spectronaut': 'Spectronaut'}
SHRINK = 5          # pseudo-count for the clustering similarity log2((n_ab + 5) / (E_ab + 5))
N_CLUST = 8


def _tree(C, n, N):
    from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform
    E = np.outer(n, n) / N; S = np.log2((C + SHRINK) / (E + SHRINK)); np.fill_diagonal(S, S.max())
    Dm = S.max() - S; np.fill_diagonal(Dm, 0); v = squareform(Dm, checks=False)
    return optimal_leaf_ordering(linkage(v, 'average'), v)


def fig_cooccurrence(d):
    from scipy.cluster.hierarchy import leaves_list, fcluster, dendrogram
    s = d[d.release_year.between(2019, 2025)]
    keys = list(SOFTWARE); X = s[['sw_' + k for k in keys]].to_numpy(float); N = len(X)
    n = X.sum(0); C = X.T @ X; E = np.outer(n, n) / N
    with np.errstate(divide='ignore'):
        L = C / E
    # group bootstrap of lift (ratio of sums over resampled contributor groups)
    cl = pd.factorize(s.cluster)[0]; Gn = cl.max() + 1
    ng = np.bincount(cl, minlength=Gn).astype(float)
    Xg = np.zeros((Gn, len(keys))); np.add.at(Xg, cl, X)
    Cg = np.zeros((Gn, len(keys) ** 2)); np.add.at(Cg, cl, (X[:, :, None] * X[:, None, :]).reshape(len(X), -1))
    rng = np.random.default_rng(SEED); reps = 2000
    Lb = np.empty((reps, len(keys), len(keys))); same = np.zeros((len(keys), len(keys))); kb = 0
    Z0 = _tree(C, n, N); lab0 = fcluster(Z0, N_CLUST, 'maxclust')
    for r in range(reps):
        w = rng.multinomial(Gn, np.full(Gn, 1 / Gn)).astype(float)
        Nb = w @ ng; nb = w @ Xg; Cb = (w @ Cg).reshape(len(keys), len(keys))
        with np.errstate(divide='ignore', invalid='ignore'):
            Lb[r] = Cb * Nb / np.outer(nb, nb)
        if r < 500:   # cluster stability: how often each pair is co-assigned when the tree is rebuilt
            lb = fcluster(_tree(Cb, nb, Nb), N_CLUST, 'maxclust'); same += lb[:, None] == lb[None, :]; kb += 1
    lo, hi = np.nanquantile(Lb, .025, axis=0), np.nanquantile(Lb, .975, axis=0)
    co = same / kb
    o = leaves_list(Z0); K = [keys[i] for i in o]; P = len(K)
    # tidy table of pairs
    rows = []
    for a in range(P):
        for b in range(a):
            i, j = o[a], o[b]
            rows.append(dict(tool_a=TOOL_NAME[K[a]], tool_b=TOOL_NAME[K[b]], datasets_both=int(C[i, j]), datasets_a=int(n[i]),
                             datasets_b=int(n[j]), expected_if_independent=E[i, j], lift=L[i, j], lift_ci_low=lo[i, j],
                             lift_ci_high=hi[i, j], jaccard=C[i, j] / (n[i] + n[j] - C[i, j]),
                             cluster_a=int(lab0[i]), cluster_b=int(lab0[j]), bootstrap_coassignment=co[i, j]))
    T = pd.DataFrame(rows); T['N_datasets'] = N
    T.to_csv(D / '18-tool-chord.csv', index=False)
    alone = {k: 100 * ((X[:, keys.index(k)] == 1) & (X.sum(1) == 1)).sum() / n[keys.index(k)] for k in keys}

    # finding (computed): strongest reliable affinity and MaxQuant vs Proteome Discoverer avoidance
    rel = T[(T.datasets_both >= 30) & (T.lift_ci_low > 1)].sort_values('lift', ascending=False).iloc[0]
    mp = T[((T.tool_a == 'MaxQuant') & (T.tool_b == 'Proteome Discoverer')) | ((T.tool_b == 'MaxQuant') & (T.tool_a == 'Proteome Discoverer'))].iloc[0]
    mq_alone = alone['maxquant']
    HM = 128; fig = plt.figure(figsize=(W, HM * MM))
    headmm(fig, f'Analysis tools come in camps: some pairs are named together up to {rel.lift:.0f}\u00d7 more often than chance',
           f'Tools named in the same PRIDE dataset (software field or processing protocol), released 2019\u20132025 (N = {N:,}; '
           f'{100 * (X.sum(1) >= 2).mean():.0f}% name \u22652 of\nthe 22 tools). Colour = lift = datasets naming both / number expected '
           f'if named independently. Grey = 95% group-bootstrap interval includes 1,\nor <5 datasets. '
           f'Order: average-linkage tree on log lift. Brackets = camps, with the share of 500 bootstrap rebuilds keeping their pairs together.', HM)
    # geometry in data units: cells 1x1, rows top->bottom
    ax = fig.add_axes([.01, .0, .98, .83])
    from matplotlib.colors import TwoSlopeNorm, LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list('div', ['#2166ac', '#67a9cf', '#d1e5f0', '#f7f7f7', '#fddbc7', '#ef8a62', '#b2182b'])
    norm = TwoSlopeNorm(0, -4, 5)
    x0 = 0.0
    for a in range(P):
        for b in range(a):
            i, j = o[a], o[b]
            ok = C[i, j] >= 5 and (lo[i, j] > 1 or hi[i, j] < 1)
            if ok:
                ax.add_patch(Rectangle((x0 + b + .05, a + .05), .9, .9, facecolor=cmap(norm(np.clip(np.log2(L[i, j]), -4, 5))), lw=0))
                v = L[i, j]
                if v >= 8 or v <= 1 / 4:
                    txt = f'{v:.0f}' if v >= 1 else f'¹⁄{1 / v:.0f}'
                    ax.text(x0 + b + .5, a + .52, txt, ha='center', va='center', fontsize=6,
                            color='white' if abs(np.log2(v)) > 3 else INK)
            else:
                ax.add_patch(Rectangle((x0 + b + .3, a + .3), .4, .4, facecolor='#d9d9d9', lw=0))
        # tool label on the diagonal
        ax.text(x0 + a + .1, a + .5, TOOL_NAME[K[a]], ha='left', va='center', fontsize=6.6, color=INK)
    # cluster brackets on the left
    labs = lab0[o]
    def cname(members):
        m = set(members)
        if 'pd' in m: return 'Proteome Discoverer\u2013\nMascot pipeline'
        if 'diann' in m and 'msfragger' in m: return 'FragPipe with\nDIA-NN'
        if 'xtandem' in m and 'msgf' in m: return 'Open-source\nsearch engines'
        return None
    start = 0
    for a in range(1, P + 1):
        if a == P or labs[a] != labs[start]:
            members = K[start:a]; nm = cname(members)
            if len(members) > 1 and nm:
                ax.plot([x0 - .35, x0 - .35], [start + .12, a - .12], color=INK, lw=.8)
                wc = [co[keys.index(p), keys.index(q)] for ii, p in enumerate(members) for q in members[ii + 1:]]
                ax.text(x0 - .6, (start + a) / 2, f'{nm}\n{100 * np.mean(wc):.0f}% stable', ha='right', va='center',
                        fontsize=6.3, color=INK, linespacing=1.15)
                stair = [(x0 + start, start + 1)]
                for r in range(start + 1, a):
                    stair += [(x0 + r, r), (x0 + r, r + 1)] if r > start + 1 else [(x0 + r, r + 1)]
                stair = [(x0 + start, start + 1), (x0 + start + 1, start + 1)]
                for r in range(start + 2, a):
                    stair += [(x0 + r - 1, r), (x0 + r, r)]
                stair += [(x0 + a - 1, a), (x0 + start, a)]
                ax.add_patch(plt.Polygon(stair, closed=True, fill=False, ec=INK, lw=.8))
            start = a
    # MaxQuant call-out
    iq = K.index('maxquant')
    mq_rel = T[((T.tool_a == 'MaxQuant') | (T.tool_b == 'MaxQuant')) & (T.datasets_both >= 5) & ((T.lift_ci_low > 1) | (T.lift_ci_high < 1))]
    ax.text(P - 9.6, iq + 4.6, f'MaxQuant ({int(n[keys.index("maxquant")]):,} datasets) is named with every other tool less often than chance\n'
            f'({(mq_rel.lift < 1).sum()} of {len(mq_rel)} reliable pairs below 1; Proteome Discoverer \u00b9\u2044{1 / mp.lift:.0f}, '
            f'95% CI \u00b9\u2044{1 / mp.lift_ci_high:.0f}\u2013\u00b9\u2044{1 / mp.lift_ci_low:.0f}); {mq_alone:.0f}% of its datasets name no other tool.\n'
            f'DIA-NN is co-named with MSFragger/FragPipe, not with Spectronaut.',
            fontsize=6.2, color=INK, va='top', ha='left', linespacing=1.25)
    # marginal: named with no other of the 22 tools, and datasets
    mx = P + 4.4; bw = 5.0
    ax.text(mx, -.4, 'Named with none\nof the other 21 tools', ha='left', va='bottom', fontsize=6.3, color=MUTED, linespacing=1.1)
    for a, k in enumerate(K):
        v = alone[k]
        ax.add_patch(Rectangle((mx, a + .2), bw * v / 100, .6, facecolor='#7f8c96', lw=0))
        ax.text(mx + bw * v / 100 + .15, a + .5, f'{v:.0f}%', va='center', fontsize=6, color=INK)
        ax.text(mx + bw + 1.6, a + .5, f'{int(n[keys.index(k)]):,}', va='center', ha='right', fontsize=6, color=MUTED)
    ax.text(mx + bw + 1.6, -.4, 'Datasets', ha='right', va='bottom', fontsize=6.3, color=MUTED)
    # colour key
    kx, ky = P - 9.6, 1.0
    for ii, lv in enumerate(np.linspace(-4, 5, 37)):
        ax.add_patch(Rectangle((kx + ii * .2, ky), .2, .5, facecolor=cmap(norm(lv)), lw=0))
    for lv, tl in [(-4, '¹⁄₁₆'), (-2, '¼'), (0, '1'), (2, '4'), (4, '16')]:
        xx = kx + (lv + 4) / 9 * 7.4
        ax.plot([xx, xx], [ky + .5, ky + .7], color=INK, lw=.5); ax.text(xx, ky + .8, tl if lv != -4 else '1/16', ha='center', va='top', fontsize=6)
    ax.text(kx, ky - .2, 'Lift (× expected if independent)', fontsize=6.3, va='bottom', color=MUTED)
    ax.set_xlim(-6.5, mx + bw + 2); ax.set_ylim(P + .2, -1.6); ax.axis('off'); ax.set_aspect('equal')
    save(fig, '18-tool-chord')
    return dict(N=N, top=(rel.tool_a, rel.tool_b, rel.lift, rel.lift_ci_low, rel.lift_ci_high, rel.datasets_both),
                mq_pd=(mp.lift, mp.lift_ci_low, mp.lift_ci_high, mp.datasets_both), mq_alone=mq_alone)


if __name__ == '__main__':
    import sys
    d = load()
    which = sys.argv[1:] or ['17', '18', '19']
    if '17' in which:
        print(fig_alluvial(d))
    if '18' in which:
        print(fig_cooccurrence(d))
