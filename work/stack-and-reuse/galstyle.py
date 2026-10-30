"""Shared style for the headline figures (gallery.py, gallery2_*.py).

Import: from galstyle import *  -> G (output dir), A (Arsenal dir), F (Deposit-features), YEARS, palette,
headline(fig, title, subtitle), save(fig, name), repel(dict, gap), yfmt(ax), W/MM sizes.

Original docstring:
Headline figures: one striking, verified finding per figure (descriptive).

Every number drawn here comes from frozen outputs under outputs/Stack-and-reuse (ledger IDs in each
figure's note) or from Deposit-features.parquet. Figures are 183 mm wide, text >= 6 pt; PNG (300 dpi) and PDF.
Outputs: outputs/Stack-and-reuse/Gallery/.
"""
import json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.transforms import Bbox
from common import O

G = O / 'Gallery'; G.mkdir(exist_ok=True)
A = O / 'Arsenal'
MM = 1 / 25.4; W = 183 * MM
INK, MUTED, HAIR, ACC, ACC2, ACC3 = '#111111', '#7a7873', '#d8d6cf', '#d1495b', '#2e6f9e', '#e8a33d'
VEND = {'Thermo Fisher': '#2e6f9e', 'Bruker': '#d1495b', 'SCIEX': '#e8a33d', 'Waters': '#5a9e6f', 'Agilent': '#8d6cab'}
plt.rcParams.update({
    'font.family': 'Liberation Sans', 'font.size': 7, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
    'axes.linewidth': .6, 'axes.edgecolor': '#333333', 'xtick.major.width': .6, 'ytick.major.width': .6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': True, 'axes.spines.right': True,
    'legend.frameon': False, 'legend.fontsize': 6.5, 'figure.facecolor': 'white', 'savefig.dpi': 300, 'pdf.fonttype': 42})
# The public release keeps its reduced analytical features under data/early.
_feature_path = O / 'Deposit-features.parquet'
if not _feature_path.exists():
    _feature_path = O.parents[1] / 'data/early/Deposit-features.parquet'
F = pd.read_parquet(_feature_path)
YEARS = list(range(2013, 2026))


# Journal mode (GAL_JOURNAL=1): titles/subtitles are recorded for the figure legend instead of drawn, panel
# letters are shifted by GAL_PANEL_OFFSET so blocks can be stacked into one multi-panel figure, and figures are
# saved to Gallery/blocks/ with a JSON sidecar (title, subtitle, letters used).
JOURNAL = os.environ.get('GAL_JOURNAL') == '1'
OFFSET = int(os.environ.get('GAL_PANEL_OFFSET', '0'))
FROZEN_DATA = os.environ.get('GAL_FROZEN_DATA') == '1'
_META = {'title': None, 'sub': None, 'letters': []}


def record_head(title, sub):
    """Store the headline for the legend; True in journal mode (caller should then draw nothing)."""
    _META['title'], _META['sub'] = title, sub
    return JOURNAL


def headline(fig, title, sub):
    if record_head(title, sub):
        return
    fig.text(.012, .975, title, fontsize=10.5, fontweight='bold', va='top', color=INK)
    fig.text(.012, .905, sub, fontsize=7, va='top', color=MUTED)


def save(fig, name):
    if JOURNAL:
        B = G / 'blocks'; B.mkdir(exist_ok=True)
        bounds = journal_layout(fig)
        _META['bounds_inches'] = list(bounds.extents)
        fig.savefig(B / f'{name}.png', dpi=600, bbox_inches=bounds, pad_inches=0)
        fig.savefig(B / f'{name}.pdf', bbox_inches=bounds, pad_inches=0)
        (B / f'{name}.json').write_text(json.dumps(_META))
        _META.update(title=None, sub=None, letters=[], notes=[])
        plt.close(fig); return
    fig.savefig(G / f'{name}.png', bbox_inches='tight', pad_inches=.04)
    fig.savefig(G / f'{name}.pdf', bbox_inches='tight', pad_inches=.04)
    plt.close(fig)



def journal_layout(fig):
    """Align panel letters by row and retain a common horizontal figure coordinate system."""
    labels = fig.findobj(lambda a: a.get_gid() == 'panel-letter')
    if not labels:
        panel(fig.axes[0], 'a')
        labels = [fig.axes[0].texts[-1]]
    for label in labels:
        label.set_visible(False)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px_mm = fig.dpi / 25.4
    positions = []
    for label in labels:
        if label.axes is not None:
            ax = label.axes
            box = ax.get_tightbbox(renderer)
            top = max(box.y1, ax.bbox.y1 + 4 * px_mm)
            left = box.x0 - 2 * px_mm
            row = ax.bbox.y1
        else:
            left, bottom = label.get_transform().transform(label.get_position())
            top = bottom + label.get_fontsize() * fig.dpi / 72
            row = top - 4 * px_mm
        positions.append([label, left, top, row])
    rows = []
    for item in sorted(positions, key=lambda p: -p[3]):
        if not rows or abs(rows[-1][0][3] - item[3]) > 3 * px_mm:
            rows.append([])
        rows[-1].append(item)
    for row in rows:
        top = max(p[2] for p in row)
        first = min(row, key=lambda p: p[1])
        for label, left, _, _ in row:
            # Insets keep their own gutter; outer rows share the same left edge.
            if label is first[0] and left < .3 * fig.bbox.width and not (label.axes and label.axes in [c for a in fig.axes for c in a.child_axes]):
                left = 1.5 * px_mm
            label.set_transform(fig.transFigure)
            label.set_position(fig.transFigure.inverted().transform((left, top)))
            label.set_ha('left' if label is first[0] else 'right')
            label.set_va('top')
            label.set_visible(True)
    fig.canvas.draw()
    box = fig.get_tightbbox(fig.canvas.get_renderer())
    pad = .8 * MM
    return Bbox.from_extents(min(0, box.x0 - pad), box.y0 - pad,
                             max(fig.get_figwidth(), box.x1 + pad), box.y1 + pad)


def repel(vals, gap):
    """Spread label positions (dict label -> y) so neighbours are at least `gap` apart."""
    items = sorted(vals.items(), key=lambda kv: kv[1]); out = {}
    ys = [v for _, v in items]
    for _ in range(200):
        moved = False
        for i in range(1, len(ys)):
            if ys[i] - ys[i - 1] < gap:
                d = (gap - (ys[i] - ys[i - 1])) / 2; ys[i] += d; ys[i - 1] -= d; moved = True
        if not moved:
            break
    return {k: y for (k, _), y in zip(items, ys)}


def panel(ax, letter, x=-0.075, y=1.0):
    """Bold lowercase panel letter at the top-left corner of an axes (Nature style)."""
    letter = chr(ord(letter) + OFFSET); _META['letters'].append(letter)
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9.5, fontweight='bold', va='bottom', ha='right', color=INK, gid='panel-letter')


def fig_letter(fig, x, y, letter, **kw):
    """Panel letter placed in figure coordinates (for panels without their own axes); shifted like panel()."""
    letter = chr(ord(letter) + OFFSET); _META['letters'].append(letter)
    fig.text(x, y, letter, fontsize=9.5, fontweight='bold', color=INK, gid='panel-letter', **kw)


def footnote(fig, x, y, text, **kw):
    """A note drawn under the figure; in journal mode it is appended to the legend text instead."""
    if JOURNAL:
        _META['notes'] = _META.get('notes', []) + [text.replace('\n', ' ')]
        return
    fig.text(x, y, text, **kw)


def yfmt(ax, pct=True):
    if pct:
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v:.0f}%'))


