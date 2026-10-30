"""Shared figure style, palette, layout fitting and publication export (183 mm, >=5 pt, PDF/TIFF/PNG)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm, TwoSlopeNorm
from matplotlib.patches import Rectangle, Patch
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from common import O

F = O / 'Figures'; F.mkdir(exist_ok=True); D = O / 'Figure-data'
FS = F / 'Supplementary'; FS.mkdir(exist_ok=True)
MM = 1 / 25.4; W = 183 * MM
INK, INK2, MUTED, HAIR = '#0b0b0b', '#52514e', '#898781', '#c3c2b7'
# Soft palette. Six tool hues at modest chroma (OKLCH C=0.105) with alternating lightness (L 0.60/0.73)
# in stacking order; this order passes the dataviz validator (adjacent CVD dE 15.0, normal-vision 19.6).
FAM = {'Mascot': '#5abd9e', 'Proteome Discoverer': '#8573b8', 'MaxQuant': '#caa056', 'FragPipe': '#af6589',
       'Spectronaut': '#73ade7', 'DIA-NN': '#b66854', 'Several or other tools': '#d3cfc7', 'No tool named': '#eeece7'}
SEQ = LinearSegmentedColormap.from_list('softblue', ['#f7f9fb', '#e3ebf4', '#c9d8ea', '#aac2de', '#8aa9cf', '#6d90bd', '#5679a6', '#45628a'])
DIV = LinearSegmentedColormap.from_list('softdiv', ['#45628a', '#7d9cc4', '#c3d3e8', '#f2f0eb', '#ebc8bc', '#c98a76', '#9c5b48'])
BLUE_D, BLUE_L, RED_D, RED_L = '#5679a6', '#d6e2f0', '#9c5b48', '#efd3c9'   # listed / not listed
plt.rcParams.update({
    'font.family': 'Liberation Sans', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5, 'xtick.labelsize': 6,
    'ytick.labelsize': 6, 'legend.fontsize': 6, 'axes.linewidth': .5, 'xtick.major.width': .5, 'ytick.major.width': .5,
    'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.direction': 'out', 'ytick.direction': 'out', 'axes.edgecolor': '#3a3a38',
    'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.grid': False, 'legend.frameon': False,
    'axes.titleweight': 'normal', 'axes.titlelocation': 'left', 'axes.titlepad': 3, 'figure.facecolor': 'white',
    'axes.facecolor': 'white', 'savefig.dpi': 600, 'pdf.fonttype': 42, 'hatch.linewidth': .4, 'hatch.color': '#b9b7b0'})


MIN_PT = 5.0        # Nature: no text below 5 pt
MAX_H = 247 * MM    # Nature: maximum page depth
SPECS = []


def fit(fig, pad=.06, iters=10):
    """Keep the canvas size fixed and move all axes so every artist (labels, legends,
    colour bars) lies inside it with `pad` inches to spare. Text keeps its point size."""
    W_, H_ = fig.get_size_inches()
    fig.set_dpi(300)  # measure text at export resolution (hinting changes extents at low dpi)
    for _ in range(iters):
        fig.canvas.draw(); r = fig.canvas.get_renderer()
        tb = fig.get_tightbbox(r)
        pos = [ax.get_position() for ax in fig.axes]
        x0 = min(p.x0 for p in pos) * W_; x1 = max(p.x1 for p in pos) * W_
        y0 = min(p.y0 for p in pos) * H_; y1 = max(p.y1 for p in pos) * H_
        nx0, nx1 = pad + (x0 - tb.x0), W_ - pad - (tb.x1 - x1)
        ny0, ny1 = pad + (y0 - tb.y0), H_ - pad - (tb.y1 - y1)
        if max(abs(nx0 - x0), abs(nx1 - x1), abs(ny0 - y0), abs(ny1 - y1)) < 2e-3:
            break
        sx, sy = (nx1 - nx0) / (x1 - x0), (ny1 - ny0) / (y1 - y0)
        for ax, p in zip(fig.axes, pos):
            ax.set_position([(nx0 + (p.x0 * W_ - x0) * sx) / W_, (ny0 + (p.y0 * H_ - y0) * sy) / H_,
                             p.width * W_ * sx / W_, p.height * H_ * sy / H_])
    fig.canvas.draw(); tb = fig.get_tightbbox(fig.canvas.get_renderer())
    return tb


def save(fig, name, folder=None):
    folder = folder or F
    W_, H_ = fig.get_size_inches()
    assert abs(W_ - W) < 1e-6 and H_ <= MAX_H + 1e-6, (W_, H_)
    tb = fit(fig)
    inside = tb.x0 >= -1e-3 and tb.y0 >= -1e-3 and tb.x1 <= W_ + 1e-3 and tb.y1 <= H_ + 1e-3
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    assert inside, (name, tb)
    assert min(sizes) >= MIN_PT - 1e-9, (name, min(sizes))
    # Exact canvas: no bbox_inches='tight', which would silently change the dimensions.
    fig.savefig(folder / f'{name}.pdf')                                   # vector, TrueType (Type 42) embedded
    fig.savefig(folder / f'{name}.tiff', dpi=600)
    from PIL import Image as _I  # flatten to RGB (no alpha channel), LZW, 600 dpi tag
    _I.open(folder / f'{name}.tiff').convert('RGB').save(folder / f'{name}.tiff', compression='tiff_lzw', dpi=(600, 600))
    fig.savefig(folder / f'{name}.png', dpi=300)                          # preview
    from PIL import Image
    im = np.asarray(Image.open(folder / f'{name}.png').convert('L'))
    edge = np.concatenate([im[:3].ravel(), im[-3:].ravel(), im[:, :3].ravel(), im[:, -3:].ravel()])
    assert edge.min() > 245, (name, 'ink touches the canvas edge')
    SPECS.append({'figure': name, 'path': folder / f'{name}.pdf', 'width_in': round(W_, 3), 'height_in': round(H_, 3), 'width_mm': round(W_ / MM, 1),
                  'height_mm': round(H_ / MM, 1), 'raster_dpi': 600, 'min_font_pt': min(sizes), 'max_font_pt': max(sizes),
                  'content_inside_canvas': bool(inside), 'formats': ['pdf (vector, Type 42 fonts)', 'tiff (600 dpi, LZW, RGB)', 'png (300 dpi preview)']})
    plt.close(fig)


def letter(ax, s, dx=-.02, dy=1.0):
    ax.annotate(s, (0, 1), xycoords='axes fraction', xytext=(-14 + dx, 6 + dy), textcoords='offset points',
                fontsize=8, fontweight='bold', va='bottom', ha='left')


def ink_for(rgba):
    r, g, b = rgba[:3]
    return 'white' if (0.299 * r + 0.587 * g + 0.114 * b) < .55 else INK


def heat(ax, M, cmap, norm, N=None, nmin=0, fmt='{:.0f}', fs=5.2, annotate=True, gap=.6):
    """Cells as rectangles with a white gap; hatched where N < nmin; values printed inside."""
    M = np.asarray(M, float); ny, nx = M.shape
    for i in range(ny):
        for j in range(nx):
            v = M[i, j]
            low = N is not None and (np.isnan(N[i, j]) or N[i, j] < nmin)
            if np.isnan(v) or low:
                ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor='white', edgecolor='#d8d6cf', hatch='//////', lw=0))
                continue
            c = cmap(norm(v)); ax.add_patch(Rectangle((j - .5, i - .5), 1, 1, facecolor=c, edgecolor='white', lw=gap))
            if annotate:
                ax.text(j, i, fmt.format(v), ha='center', va='center', fontsize=fs, color=ink_for(c))
    ax.set_xlim(-.5, nx - .5); ax.set_ylim(ny - .5, -.5)
    ax.tick_params(length=0)


def cbar(fig, ax, cmap, norm, label, ticks=None, orientation='horizontal'):
    sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, cax=ax, orientation=orientation, ticks=ticks)
    cb.outline.set_linewidth(.4); cb.ax.tick_params(length=1.5, width=.4, labelsize=5.5)
    cb.set_label(label, fontsize=5.8, labelpad=1.5)
    return cb




def write_specs(path):
    import json, subprocess
    for sp in SPECS:
        out = subprocess.run(['pdffonts', str(sp['path'])], capture_output=True, text=True).stdout.splitlines()[2:]
        sp['pdf_fonts'] = sorted({l.split()[0] for l in out}); sp['pdf_fonts_embedded'] = all(' yes ' in l for l in out)
        assert sp['pdf_fonts_embedded'], sp
        sp['path'] = str(sp['path'].relative_to(O))
    path.write_text(json.dumps(SPECS, indent=2))
