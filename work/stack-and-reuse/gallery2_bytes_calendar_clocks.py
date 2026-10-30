"""Gallery 2, figures 23-25: where the bytes are, when datasets go public, and the clocks of a dataset (descriptive).

23 23-petabyte-treemap   stored bytes by file layer (squarified treemap, own implementation) with a linked
                         treemap of uncompressed member bytes inside the listed results zips.
24 24-release-calendar   datasets made public per day, 2018-2025 (calendar heatmap, log colour), bulk days
                         (>= 50 releases/day) and, per release year, median paper -> release interval.
25 25-life-of-a-dataset  the clocks of a dataset relative to public release, in event order (half-eye plots).

Inputs (frozen): Arsenal/H4-H7-databases-storage/H7-*.csv, Archive-contents/Denominators.csv,
Attacks/release-timing/B1-bulk-days.csv, Arsenal/H11-clocks/Dataset-clocks.csv.gz, Deposit-features.parquet.
Uncertainty: 95% percentile interval over 2000 resamples of whole pi_group clusters (datasets without a
group hash are their own cluster), seed 20261001. Plotted tables -> Gallery/data/<NN>-<slug>.csv.
Run from work/stack-and-reuse/: python gallery2_bytes_calendar_clocks.py
"""
import json
import numpy as np, pandas as pd
from matplotlib.patches import Rectangle, Polygon, ConnectionPatch
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm, to_rgb
from galstyle import *
from matplotlib.ticker import FuncFormatter
from galstyle import INK, MUTED, HAIR, ACC, ACC2, ACC3
from common import O, SEED

D = G / 'data'; D.mkdir(exist_ok=True)
REPS = 2000


# ------------------------------------------------------------------ helpers ----------------------------------
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


def squarify(values, x, y, w, h):
    """Squarified treemap (Bruls, Huizing & van Wijk 2000). values sorted descending, > 0.
    Returns list of (x, y, w, h) in the same order."""
    vals = np.asarray(values, float); vals = vals * (w * h) / vals.sum()
    rects = []

    def worst(row, side):
        s = sum(row); return max(max(side ** 2 * r / s ** 2, s ** 2 / (side ** 2 * r)) for r in row)

    def layout(row, x, y, w, h):
        s = sum(row); out = []
        if w >= h:                       # lay the row as a column on the left
            cw = s / h; yy = y
            for r in row:
                out.append((x, yy, cw, r / cw)); yy += r / cw
            return out, (x + cw, y, w - cw, h)
        rh = s / w; xx = x                 # lay the row along the top
        for r in row:
            out.append((xx, y + h - rh, r / rh, rh)); xx += r / rh
        return out, (x, y, w, h - rh)

    i, row = 0, []
    while i < len(vals):
        side = min(w, h); r = vals[i]
        if not row or worst(row + [r], side) <= worst(row, side):
            row.append(r); i += 1
        else:
            out, (x, y, w, h) = layout(row, x, y, w, h); rects += out; row = []
    if row:
        out, _ = layout(row, x, y, w, h); rects += out
    return rects


def tb(v):
    return f'{v:,.0f} TB' if v >= 10 else f'{v:.1f} TB'




# ------------------------------------------------------------------ 23 treemap --------------------------------
BL, BL2, BL3, TOOLC, SCANC, GRY = '#2e6f9e', '#6b98bf', '#b3cadd', '#8f8aa3', '#a39e92', '#d3d0c8'


def fig_treemap():
    H7 = A / 'H4-H7-databases-storage'
    lay = pd.read_csv(H7 / 'H7-top-level-bytes-by-layer.csv')
    S = json.load(open(H7 / 'Summary.json'))['H7_anatomy']
    den = pd.read_csv(O / 'Archive-contents/Denominators.csv').set_index('why')
    tot = lay.TB.sum(); nfiles = lay.files.sum()
    lay['block'] = np.select(
        [lay.layer.eq('vendor_raw') & ~lay.is_archive, lay.layer.eq('vendor_raw') & lay.is_archive, lay.is_archive,
         lay.layer.eq('tool_results'), lay.layer.eq('open_spectra'), lay.layer.eq('standard_results')],
        ['raw_loose', 'raw_zip', 'arch', 'tool', 'spectra', 'std'], 'rest')
    top = lay.groupby('block').agg(TB=('TB', 'sum'), files=('files', 'sum')).sort_values('TB', ascending=False)
    top['pct'] = 100 * top.TB / tot
    col = {'raw_loose': BL, 'raw_zip': BL2, 'spectra': BL3, 'tool': TOOLC, 'arch': '#3d3b37', 'std': ACC, 'rest': GRY}

    z = pd.read_csv(H7 / 'H7-inside-zips-bytes-by-role.csv')
    ztot = z.TB_uncompressed.sum()
    keep = ['scan-level and intermediate search files', 'raw data', 'tool project and result databases', 'final result tables',
            'open spectra (mzML, MGF...)', 'other tables', 'standard result files (mzIdentML/mzTab)']
    z['cell'] = np.where(z.role.isin(keep), z.role, 'everything else')
    zc = z.groupby('cell').agg(TB=('TB_uncompressed', 'sum'), members=('members', 'sum')).sort_values('TB', ascending=False)
    zc['pct'] = 100 * zc.TB / ztot
    zcol = {'scan-level and intermediate search files': SCANC, 'raw data': BL, 'tool project and result databases': TOOLC,
            'final result tables': ACC, 'open spectra (mzML, MGF...)': BL3, 'other tables': GRY,
            'standard result files (mzIdentML/mzTab)': ACC, 'everything else': '#e6e3dc'}

    raw_pct = top.loc[['raw_loose', 'raw_zip'], 'pct'].sum()
    fin, std, scan = (zc.loc[k] for k in ['final result tables', 'standard result files (mzIdentML/mzTab)', 'scan-level and intermediate search files'])
    listed_tb, targ_tb = den.tb['listed'], den.tb.sum()
    rest_parts = lay[lay.block.eq('rest')].groupby('layer').TB.sum().sort_values(ascending=False)

    pd.concat([
        lay.assign(panel='a', bytes_kind='stored', pct=100 * lay.TB / tot)[['panel', 'block', 'layer', 'is_archive', 'files', 'TB', 'pct', 'bytes_kind']],
        z.assign(panel='b', bytes_kind='uncompressed (inside listed zips)', pct=100 * z.TB_uncompressed / ztot, block=z.cell,
                 layer=z.role, is_archive=np.nan, files=z.members, TB=z.TB_uncompressed)[['panel', 'block', 'layer', 'is_archive', 'files', 'TB', 'pct', 'bytes_kind']]
    ]).round(3).to_csv(D / '23-petabyte-treemap.csv', index=False)

    FH = 114
    fig = plt.figure(figsize=(W, FH * MM))
    headline(fig, f'Raw data fill {raw_pct:.0f}% of PRIDE’s {tot / 1000:.2f} petabytes; final result tables are '
                  f'{fin.pct:.0f}% of what results zips hold',
             f'(a) Stored bytes of all {nfiles / 1e6:.2f} million top-level files in the {len(F):,} study datasets, by file layer. (b) Uncompressed bytes inside the\n'
             f'{S["listed_zips"]:,} results zips whose member list could be read ({listed_tb:.1f} of {targ_tb:.1f} TB of non-raw archives; rar, 7z and tar cannot be listed),\n'
             f'by role of each member file (assigned from file names). The two panels use different scales.')

    def tmap(rect):
        a_ = fig.add_axes(rect); wmm, hmm = rect[2] * 183, rect[3] * FH
        a_.set_xlim(0, wmm); a_.set_ylim(0, hmm); a_.axis('off'); return a_, wmm, hmm
    # ---- panel a -----------------------------------------------------------------------------------
    ax, wa, ha = tmap([.012, .115, .488, .63])
    order = ['raw_loose', 'raw_zip', 'tool', 'spectra', 'std', 'rest', 'arch']   # descending except the three slivers,
    assert set(order) == set(top.index)                                             # ordered so the archives sit at the edge facing (b)
    R = dict(zip(order, squarify(top.TB[order].to_numpy(), 0, 0, wa, ha)))
    for k, (x, y, w, h) in R.items():
        ax.add_patch(Rectangle((x, y), w, h, facecolor=col[k], edgecolor='white', lw=1.4 if top.pct[k] > 1 else .5))

    def lab(A_, RR, k, txt, big=False, color='white', fs=6.6):
        x, y, w, h = RR[k]
        A_.text(x + 1.5, y + h - 1.6, txt, color=color, fontsize=9 if big else fs, va='top', ha='left',
                fontweight='bold' if big else 'normal', linespacing=1.12)
    lab(ax, R, 'raw_loose', f'Vendor raw files\n{top.TB.raw_loose:,.0f} TB · {top.pct.raw_loose:.0f}%', big=True)
    lab(ax, R, 'raw_zip', f'Vendor raw files,\nzipped\n{top.TB.raw_zip:,.0f} TB · {top.pct.raw_zip:.0f}%', big=True)
    lab(ax, R, 'tool', f'Analysis-tool\noutputs\n{top.TB.tool:,.0f} TB · {top.pct.tool:.1f}%')
    lab(ax, R, 'spectra', f'Open\nspectra\n{top.TB.spectra:,.0f} TB\n{top.pct.spectra:.1f}%', color=INK)
    # callouts below panel a for the two slivers (leaders leave from the cells' lower edges)
    x, y, w, h = R['std']
    ax.annotate(f'Standard result files\n(mzIdentML, mzTab)  {top.TB["std"]:.1f} TB \u00b7 {top.pct["std"]:.2f}%', (x + w / 2, y + .3), (x - 1.5, -3.5),
                fontsize=6.3, color=ACC, ha='right', va='top', linespacing=1.1, fontweight='bold', annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=ACC, lw=.5, shrinkA=0, shrinkB=0))
    x, y, w, h = R['rest']
    ax.annotate(f'Everything else (tables, FASTA,\nlibraries, documents)  {top.TB.rest:.0f} TB \u00b7 {top.pct.rest:.1f}%', (x + w / 2, y + .3), (x + w / 2 + 1.5, -3.5),
                fontsize=6.2, color=MUTED, ha='left', va='top', linespacing=1.1, annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=MUTED, lw=.5, shrinkA=0, shrinkB=0))
    xa, ya, wa_, ha_ = R['arch']
    lab(ax, R, 'arch', f'Non-raw archives\n{top.TB.arch:.0f} TB \u00b7 {top.pct.arch:.1f}%', fs=6.2)
    panel(ax, 'a', x=3.2 / wa, y=(ha + 1.2) / ha)
    ax.text(5, ha + 1.5, f'Stored bytes, {tot / 1000:.2f} PB', fontsize=7, va='bottom', color=INK)
    # ---- panel b -----------------------------------------------------------------------------------
    bx, wb, hb = tmap([.635, .115, .35, .63])
    ZR = dict(zip(zc.index, squarify(zc.TB.to_numpy(), 0, 0, wb, hb)))
    for k, (x, y, w, h) in ZR.items():
        bx.add_patch(Rectangle((x, y), w, h, facecolor=zcol[k], edgecolor='white', lw=1.2 if zc.pct[k] > 2 else .5))
    lab(bx, ZR, 'scan-level and intermediate search files', f'Scan-level and intermediate\nsearch files\n{scan.TB:.0f} TB · {scan.pct:.0f}%', big=True)
    lab(bx, ZR, 'raw data', f'Raw data\n{zc.TB["raw data"]:.0f} TB\n{zc.pct["raw data"]:.0f}%')
    lab(bx, ZR, 'tool project and result databases', f'Tool\nproject\ndatabases\n{zc.pct["tool project and result databases"]:.1f}%')
    lab(bx, ZR, 'final result tables', f'Final result\ntables\n{fin.TB:.1f} TB\n{fin.pct:.1f}%', fs=6.3)
    lab(bx, ZR, 'open spectra (mzML, MGF...)', f'Open\nspectra\n{zc.pct["open spectra (mzML, MGF...)"]:.1f}%', color=INK, fs=6.2)
    lab(bx, ZR, 'other tables', f'Other\ntables\n{zc.pct["other tables"]:.1f}%', color=INK, fs=6.2)
    lab(bx, ZR, 'everything else', f'Logs, docs,\nlibraries,\nother {zc.pct["everything else"]:.1f}%', color=INK, fs=6.2)
    x, y, w, h = ZR['standard result files (mzIdentML/mzTab)']
    bx.annotate(f'Standard result files {std.TB:.2f} TB · {std.pct:.1f}%', (x + w / 2, y + h * .4), (wb, -4), fontsize=6.3, color=ACC,
                fontweight='bold', ha='right', va='top', linespacing=1.1, annotation_clip=False,
                arrowprops=dict(arrowstyle='-', color=ACC, lw=.5, shrinkA=0, shrinkB=0))
    panel(bx, 'b', x=3.2 / wb, y=(hb + 1.2) / hb)
    bx.text(5, hb + 1.5, f'Inside listed results zips, {ztot:.1f} TB uncompressed', fontsize=7, va='bottom', color=INK)
    # zoom wedge from the archive block to panel b, drawn behind both treemaps
    from matplotlib.patches import Polygon as Poly
    inv = fig.transFigure.inverted()
    pts = [inv.transform(ax.transData.transform(p_)) for p_ in [(xa + wa_, ya + ha_), (xa + wa_, ya)]] + \
          [inv.transform(bx.transData.transform(p_)) for p_ in [(0, 0), (0, hb)]]
    fig.patches.append(Poly([pts[0], pts[3], pts[2], pts[1]], closed=True, transform=fig.transFigure, facecolor='#eeece6', edgecolor='none', zorder=-1))
    ax.add_patch(Rectangle((xa, ya), wa_, ha_, facecolor='none', edgecolor=INK, lw=.9))
    fig.text((pts[0][0] + pts[3][0]) / 2 + .004, pts[1][1] + .16, f'{listed_tb:.0f} TB of\nlistable zips\nunpack to\n{ztot:.0f} TB', fontsize=6.3,
             color='#5e5b55', ha='center', va='center', linespacing=1.15)
    save(fig, '23-petabyte-treemap')
    return dict(total_PB=tot / 1000, raw_pct=raw_pct, final_pct=fin.pct, std_pct_zip=std.pct, std_pct_top=top.pct['std'],
                scan_pct=scan.pct, listed_TB=listed_tb, target_TB=targ_tb, rest=rest_parts.to_dict())


# ------------------------------------------------------------------ 24 release calendar -----------------------
def fig_calendar():
    Y0, Y1 = 2018, 2025
    b = pd.read_csv(O / 'Attacks/release-timing/B1-bulk-days.csv', parse_dates=['release_date'])
    bulk = set(b.release_date)
    rd = pd.to_datetime(F.release_date).dt.normalize()
    cnt_all = rd.value_counts()                                   # bulk definition counts the full 36,249 cohort
    days = pd.date_range(f'{Y0}-01-01', f'{Y1}-12-31', freq='D')
    cnt = cnt_all.reindex(days).fillna(0).astype(int)
    assert set(cnt.index[cnt >= 50]) == bulk, 'bulk-day definition does not reproduce B1'
    n_period = int(cnt.sum()); n_bulk = int(cnt[cnt >= 50].sum())
    # paper -> release interval (H11 clocks; same population as the release-timing attack: lag >= 0, dated OpenAlex paper)
    c = pd.read_csv(A / 'H11-clocks/Dataset-clocks.csv.gz', usecols=['accession', 'release_date', 'release_year', 'grp', 'paper_to_rel'])
    c['release_date'] = pd.to_datetime(c.release_date); c['bulk'] = c.release_date.isin(bulk)
    p = c[c.paper_to_rel.notna()].copy()
    rows = []
    for yr in ['all'] + list(range(Y0, Y1 + 1)):
        z = p if yr == 'all' else p[p.release_year == yr]
        for bk in [True, False]:
            q = z[z.bulk == bk]
            if len(q) < 20:
                continue
            lo, hi = group_boot_median(q.paper_to_rel, q.grp)
            rows.append({'release_year': yr, 'bulk_day': bk, 'datasets': len(q), 'groups': q.grp.nunique(),
                         'median_days': q.paper_to_rel.median(), 'ci_lo': lo, 'ci_hi': hi,
                         'share_gt30d_after': (q.paper_to_rel > 30).mean()})
    T = pd.DataFrame(rows)
    ov = T[T.release_year == 'all'].set_index('bulk_day')
    assert ov.median_days[True] == 95 and ov.median_days[False] == 11
    per = pd.DataFrame({'year': range(Y0, Y1 + 1)})
    per['bulk_days'] = [sum(1 for d in bulk if d.year == y) for y in per.year]
    per['released'] = [int(cnt[str(y)].sum()) for y in per.year]
    per['on_bulk_days'] = [int(cnt[str(y)][cnt[str(y)] >= 50].sum()) for y in per.year]
    per['share_on_bulk'] = per.on_bulk_days / per.released
    cal = pd.DataFrame({'date': days, 'datasets': cnt.to_numpy(), 'bulk_day': cnt.to_numpy() >= 50})
    cal.to_csv(D / '24-release-calendar.csv', index=False)
    T.round(4).to_csv(D / '24-release-calendar-paper-to-release.csv', index=False)
    per.round(4).to_csv(D / '24-release-calendar-per-year.csv', index=False)

    # colour: 0 = paper; 1-49 blue ramp (log); >= 50 red ramp (log)
    blues = LinearSegmentedColormap.from_list('b', ['#e3ebf2', '#9fb9d0', '#4f7fa8', '#24527a'])
    reds = LinearSegmentedColormap.from_list('r', ['#f0a5a0', '#d1495b', '#7a1424'])
    vmax = cnt.max()
    def colour(n):
        if n == 0:
            return '#f3f2ee'
        if n < 50:
            return blues(np.log(n) / np.log(49))
        return reds((np.log(n) - np.log(50)) / (np.log(vmax) - np.log(50)))

    FH = 128
    fig = plt.figure(figsize=(W, FH * MM))
    headline(fig, f'{100 * n_bulk / n_period:.0f}% of datasets released since 2018 went public on {len(bulk)} days, typically months after their paper',
             f'(a) PRIDE datasets made public per day, {Y0}–{Y1} (n = {n_period:,}); red: the {len(bulk)} days with ≥50 releases. '
             f'(b) Median days from the linked paper (OpenAlex date)\nto public release, '
             f'for datasets released on those days vs other days (n = {len(p[p.release_year >= Y0]):,} with a dated paper); bars: 95% CI over contributor groups.')
    cell = 1.78                     # mm
    x0, y0 = 9, 6                   # mm offsets inside the calendar axes
    blk_h = 7 * cell + 3.2
    cw, ch = x0 + 54 * cell + 1, y0 + (Y1 - Y0 + 1) * blk_h
    axc = fig.add_axes([.012, .02, cw / 183, ch / FH]); axc.set_xlim(0, cw); axc.set_ylim(ch, 0); axc.axis('off')
    patches, cols = [], []
    for yi, yr in enumerate(range(Y0, Y1 + 1)):
        top_ = y0 + yi * blk_h
        start = pd.Timestamp(f'{yr}-01-01'); wk0 = start - pd.Timedelta(days=start.weekday())
        yd = pd.date_range(start, f'{yr}-12-31')
        for d in yd:
            wk = (d - wk0).days // 7; wd = d.weekday()
            patches.append(Rectangle((x0 + wk * cell, top_ + wd * cell), cell * .88, cell * .88)); cols.append(colour(cnt[d]))
        # month boundaries (step line between months)
        for m in range(2, 13):
            d = pd.Timestamp(yr, m, 1); wk = (d - wk0).days // 7; wd = d.weekday()
            xs = x0 + wk * cell - cell * .06
            pts = [(xs + cell, top_ - .1), (xs + cell, top_ + wd * cell - cell * .06), (xs, top_ + wd * cell - cell * .06), (xs, top_ + 7 * cell)] if wd else \
                  [(xs, top_ - .1), (xs, top_ + 7 * cell)]
            axc.plot(*zip(*pts), color='white', lw=1.1, solid_joinstyle='miter', zorder=3)
        axc.text(x0 - 1.2, top_ + 3.5 * cell, str(yr), ha='right', va='center', fontsize=7, fontweight='bold', color=INK)
    axc.add_collection(PatchCollection(patches, facecolor=cols, edgecolor='none'))
    for i, m in enumerate('JFMAMJJASOND'):
        d = pd.Timestamp(Y0, i + 1, 1); wk0 = pd.Timestamp(f'{Y0}-01-01') - pd.Timedelta(days=pd.Timestamp(f'{Y0}-01-01').weekday())
        axc.text(x0 + ((d - wk0).days // 7 + 2) * cell, y0 - 1, m, fontsize=6, color=MUTED, ha='center', va='bottom')
    for wd, l in [(0, 'M'), (2, 'W'), (4, 'F'), (6, 'S')]:
        axc.text(x0 + 54 * cell - .2, y0 + wd * cell + cell * .44, l, fontsize=6, color=MUTED, va='center', ha='left')
    # label the largest day
    mx = cnt.idxmax(); yi = mx.year - Y0; wk0 = pd.Timestamp(mx.year, 1, 1) - pd.Timedelta(days=pd.Timestamp(mx.year, 1, 1).weekday())
    px, py = x0 + ((mx - wk0).days // 7 + .44) * cell, y0 + yi * blk_h + (mx.weekday() + .44) * cell
    axc.annotate(f'{cnt[mx]} on {mx.day} {mx:%b}', (px, py), (px + 9, y0 + yi * blk_h - 1.2), fontsize=6.2, color=ACC, fontweight='bold',
                 va='bottom', arrowprops=dict(arrowstyle='-', color=ACC, lw=.5, shrinkA=0, shrinkB=1))
    panel(axc, 'a', x=3 / cw, y=1 - 1.5 / ch)
    # per-year text: bulk days and share of that year's releases on them
    xt = x0 + 54 * cell + 6.5
    axc.text(x0 + 54 * cell + .4, y0 + .44 * cell, 'Mon', fontsize=6, color=MUTED, va='center')
    axc.text(x0 + 54 * cell + .4, y0 + 6.44 * cell, 'Sun', fontsize=6, color=MUTED, va='center')
    for yi, r in per.iterrows():
        top_ = y0 + yi * blk_h
        if r.bulk_days:
            axc.text(xt, top_ + 2.2 * cell, f'{r.bulk_days} bulk day{"s" if r.bulk_days > 1 else ""}', fontsize=6.3, color=ACC, va='center', fontweight='bold')
            axc.text(xt, top_ + 4.6 * cell, f'{100 * r.share_on_bulk:.0f}% of {r.released:,}', fontsize=6.2, color=MUTED, va='center')
        else:
            axc.text(xt, top_ + 3.4 * cell, f'none · {r.released:,}', fontsize=6.2, color=MUTED, va='center')
    # legend in the spare bottom row
    ly = y0 + (Y1 - Y0 + 1) * blk_h + 1.5
    steps = [(0, '0'), (1, '1'), (3, '3'), (10, '10'), (30, '30'), (50, '50'), (100, '100'), (200, '200'), (int(vmax), f'{int(vmax)}')]
    lx = x0
    axc.text(x0 - 1.2, ly + cell * .44, 'per day', fontsize=6, color=MUTED, ha='right', va='center')
    for i, (v, l) in enumerate(steps):
        xx = lx + i * cell * 2.3 + (2.2 if v >= 50 else 0)
        axc.add_patch(Rectangle((xx, ly), cell * 1.6, cell * .88, facecolor=colour(v), edgecolor='none'))
        axc.text(xx + cell * .8, ly + cell * 1.1, l, fontsize=6, color=MUTED, ha='center', va='top')
    axc.text(lx + 9 * cell * 2.3 + 3.4, ly + cell * .44, '≥50 = bulk day (red scale)', fontsize=6, color=ACC, va='center')
    axc.set_ylim(ly + 5, 0)
    # ---- panel b: median paper -> release, bulk vs other days, aligned to each year block -------------------
    bxl = .735
    bx = fig.add_axes([bxl, .02 + (ch - (ly + 5)) / FH * 0, .985 - bxl, (ly + 5) / FH])
    bx.set_ylim(ly + 5, 0); bx.set_xscale('log'); bx.set_xlim(3, 1000)
    yc = {yr: y0 + i * blk_h + 3.5 * cell for i, yr in enumerate(range(Y0, Y1 + 1))}
    yc['all'] = ly + 1.5
    for yr, yy in yc.items():
        q = T[T.release_year == yr].set_index('bulk_day')
        if True in q.index and False in q.index:
            bx.plot([q.median_days[False], q.median_days[True]], [yy, yy], color=HAIR, lw=1.6, zorder=1, solid_capstyle='butt')
        for bk, cc in [(False, ACC2), (True, ACC)]:
            if bk not in q.index:
                continue
            m, lo, hi = q.median_days[bk], q.ci_lo[bk], q.ci_hi[bk]
            big = yr == 'all'
            bx.plot([max(lo, 3), hi], [yy, yy], color=cc, lw=2.6 if big else 1.8, alpha=.35, solid_capstyle='butt', zorder=2)
            bx.plot(m, yy, 'o', ms=5 if big else 4, color=cc, mec='white', mew=.6, zorder=3)
            off = 1.35 if bk else 1 / 1.35
            bx.text(m * off if bk else lo / 1.12, yy, f'{m:.0f}', fontsize=7.5 if big else 6.4, color=cc, va='center',
                    ha='left' if bk else 'right', fontweight='bold' if big else 'normal')
    bx.axhline(ly - 1.6, color=HAIR, lw=.6)
    bx.text(3.3, yc['all'], '2018–25', fontsize=6.4, fontweight='bold', va='center', color=INK)
    bx.set_yticks([]); bx.set_xticks([3, 10, 30, 100, 300, 1000]); bx.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:g}'))
    bx.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: '')); bx.tick_params(axis='x', labelsize=6.2)
    bx.xaxis.tick_top(); bx.xaxis.set_label_position('top')
    bx.set_xlabel('Median days from paper to public release', fontsize=6.4, labelpad=2)
    bx.text(1.0, -0.005, '', transform=bx.transAxes)
    bx.text(np.sqrt(11 * 95), y0 + 7 * blk_h + 3.5 * cell + 6.3, '', fontsize=6)
    bx.text(.97, .995, 'bulk days', transform=bx.transAxes, fontsize=6.4, color=ACC, ha='right', va='top', fontweight='bold')
    bx.text(.03, .995, 'other days', transform=bx.transAxes, fontsize=6.4, color=ACC2, ha='left', va='top', fontweight='bold')
    panel(bx, 'b', x=-.03, y=1.0)
    save(fig, '24-release-calendar')
    return dict(n_period=n_period, n_bulk=n_bulk, share=n_bulk / n_period, bulk_days=len(bulk),
                overall=ov[['datasets', 'groups', 'median_days', 'ci_lo', 'ci_hi', 'share_gt30d_after']].to_dict('index'))


if __name__ == '__main__':
    import sys
    which = sys.argv[1:] or ['23', '24', '25']
    if '23' in which:
        print(fig_treemap())
    if '24' in which:
        print(fig_calendar())
