"""Assemble the gallery figures into journal-ready multi-panel figures (descriptive; no new numbers).

Each gallery figure is re-rendered as a *block* in journal mode (galstyle GAL_JOURNAL=1): its headline and
subtitle are not drawn but recorded for the legend, and its panel letters are shifted (GAL_PANEL_OFFSET) so
that blocks stacked into one figure carry continuous letters a, b, c ... A block with no lettered panels gets
one letter from this script. Blocks are stacked at full width (7.20 inches) up to 9.72 inches high, rendered at 600 dpi.

  python compose_journal.py            render blocks, assemble Main/ and Supplementary/ figures + Legends.md
  python compose_journal.py --no-render reuse rendered blocks
  python compose_journal.py --frozen-data preserve saved bootstrap tables during layout edits

Outputs: outputs/Stack-and-reuse/Journal-figures/{Main,Supplementary}/Figure-N.{png,pdf}, Legends.md.
"""
import json, os, re, subprocess, sys
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import O

HERE = Path(__file__).resolve().parent
G = O / 'Gallery'; B = G / 'blocks'
OUT = O / 'Journal-figures'
MM = 1 / 25.4; W_IN = 7.2; MAX_H_IN = 9.72; DPI = 600
W_MM = W_IN * 25.4; MAX_H_MM = MAX_H_IN * 25.4

# slug -> command that renders that gallery figure alone
GAL = lambda fn: [sys.executable, '-c', f'import gallery; gallery.{fn}()']
RUN = {
    '01-instruments': GAL('fig_instruments'), '02-run-size': GAL('fig_run_size'), '03-scale': GAL('fig_scale'),
    '04-newcomers': GAL('fig_newcomers'), '05-software': GAL('fig_software'), '06-standards-by-software': GAL('fig_standards'),
    '07-hometowns': GAL('fig_hometowns'), '08-countries': GAL('fig_countries'), '09-clocks': GAL('fig_clocks'),
    '10-release-waves': GAL('fig_release_waves'), '11-bytes': GAL('fig_bytes'), '12-hidden-results': GAL('fig_hidden_results'),
    '13-habits': GAL('fig_habits'), '26-reuse': GAL('fig_reuse'), '42-breadth': [sys.executable, 'gallery4_breadth.py'],
    '14-fleet-history': [sys.executable, 'gallery2_scale.py', '14'], '15-run-ridgelines': [sys.executable, 'gallery2_scale.py', '15'],
    '17-stack-alluvial': [sys.executable, 'gallery2_stack.py', '17'], '20-world-growth': [sys.executable, 'gallery2_geography.py', '20'],
    '23-petabyte-treemap': [sys.executable, 'gallery2_bytes_calendar_clocks.py', '23'],
    '27-isobaric-plex': [sys.executable, 'gallery3_methods.py', '27'], '28-ptm-landscape': [sys.executable, 'gallery3_methods.py', '28'],
    '29-species': [sys.executable, 'gallery3_methods.py', '29'], '30-file-formats': [sys.executable, 'gallery3_methods.py', '30'],
    '31-version-streams': [sys.executable, 'gallery3_versions.py', '31'], '33-tool-age-at-release': [sys.executable, 'gallery3_versions.py', '33'],
    '32-database-age': [sys.executable, 'gallery3_databases.py'],
    '34-contributor-concentration': [sys.executable, 'gallery3_producers.py', '34'], '35-collaboration': [sys.executable, 'gallery3_producers.py', '35'],
    '36-research-fields': [sys.executable, 'gallery3_producers.py', '36'], '37-journals': [sys.executable, 'gallery3_producers.py', '37'],
    '38-instrument-arrival': [sys.executable, 'gallery3_timing.py', '38'], '39-data-age-by-country': [sys.executable, 'gallery3_timing.py', '39'],
    '40-weekly-rhythm': [sys.executable, 'gallery3_timing.py', '40'], '41-paper-and-release': [sys.executable, 'gallery3_timing.py', '41'],
}

# Paper: "The data economy of public proteomics: thirteen years of PRIDE". Blocks per figure, top to bottom.
MAIN = [
    ('An expanding community contributing to PRIDE', ['42-breadth', '04-newcomers', '20-world-growth']),
    ('A growing archive spanning different experimental scales', ['15-run-ridgelines', '02-run-size', '11-bytes']),
    ('A public record of instruments and analysis workflows', ['01-instruments', '05-software', '07-hometowns']),
    ('Analytical content shared in standard and tool-specific formats', ['06-standards-by-software', '12-hidden-results', '13-habits']),
    ('Connecting public deposits to the research literature', ['09-clocks', '37-journals', '26-reuse']),
]
EXTENDED = [
    ('Concentration of supply', ['34-contributor-concentration']),
    ('Fleet history of every common instrument model', ['14-fleet-history']),
    ('Arrival of new instrument models by continent', ['38-instrument-arrival']),
    ('Software release lines in use', ['31-version-streams']),
    ('Software and sequence-database age', ['33-tool-age-at-release', '32-database-age']),
    ('Where the petabytes are', ['23-petabyte-treemap']),
    ('File formats', ['30-file-formats']),
    ('Isobaric labelling and modifications', ['27-isobaric-plex', '28-ptm-landscape']),
    ('Organisms', ['29-species']),
    ('Countries', ['08-countries', '39-data-age-by-country']),
    ('Collaboration', ['35-collaboration']),
    ('Research fields of the linked papers', ['36-research-fields']),
    ('Release timing', ['10-release-waves', '41-paper-and-release']),
    ('The submission clock', ['40-weekly-rhythm']),
]


def shift_letters(text, k):
    """Re-letter panel references written in a block's own subtitle/notes ((a), a, / b: / b–c) by k letters."""
    if not k or not text:
        return text
    sh = lambda c: chr(ord(c) + k)
    pat = r'\(([a-h])\)|(?<![\w’\'.])([a-h])–([a-h])(?![\w])|(?<![\w’\'.–])([a-h])(?=[,:])'

    def f(m):
        if m.group(1):
            return f'({sh(m.group(1))})'
        if m.group(2):
            return f'{sh(m.group(2))}–{sh(m.group(3))}'
        return sh(m.group(4))
    text = re.sub(pat, f, text)
    return text


def render(slug, offset):
    env = dict(os.environ, GAL_JOURNAL='1', GAL_PANEL_OFFSET=str(offset), MPLBACKEND='Agg')
    r = subprocess.run(RUN[slug], cwd=HERE, env=env, capture_output=True, text=True)
    if r.returncode != 0 or not (B / f'{slug}.json').exists():
        raise RuntimeError(f'{slug}: {r.stderr[-800:]}')
    return json.loads((B / f'{slug}.json').read_text())


def block_size_mm(slug):
    im = Image.open(B / f'{slug}.png'); return im.size[0] / DPI * 25.4, im.size[1] / DPI * 25.4


def assemble(kind, n, title, slugs, do_render):
    """Render each block with the right letter offset, stack them, write PNG/PDF; return legend text."""
    offset, metas = 0, []
    for slug in slugs:
        meta = render(slug, offset) if do_render else json.loads((B / f'{slug}.json').read_text())
        own = len(sorted(set(meta['letters'])))
        meta['slug'], meta['first'] = slug, chr(ord('a') + offset)
        meta['add_letter'] = own == 0
        meta['n_letters'] = max(own, 1)
        offset += meta['n_letters']; metas.append(meta)
    sizes = [block_size_mm(s) for s in slugs]
    gap = 5.0
    bounds = [m.get('bounds_inches', [0, 0, w * MM, h * MM])
              for m, (w, h) in zip(metas, sizes)]
    left = min(b[0] for b in bounds); right = max(b[2] for b in bounds)
    canvas_w = (right - left) / MM
    heights = [(b[3] - b[1]) / MM for b in bounds]
    H = sum(heights) + gap * (len(slugs) - 1)
    scale = min(1.0, W_MM / canvas_w, MAX_H_MM / H)
    fig = plt.figure(figsize=(W_IN, H * scale * MM))
    y = H
    for slug, h, b, m in zip(slugs, heights, bounds, metas):
        y -= h
        x = (b[0] - left) / MM * scale / W_MM
        width = (b[2] - b[0]) / MM * scale / W_MM
        ax = fig.add_axes([x, y / H, width, h / H])
        ax.imshow(Image.open(B / f'{slug}.png')); ax.axis('off')
        if m['add_letter']:
            fig.text(.008, (y + h) / H - .002, m['first'], fontsize=9.5 * scale,
                     fontweight='bold', va='top', ha='left')
        y -= gap
    d = OUT / kind; d.mkdir(parents=True, exist_ok=True)
    stem = f'Figure-{n}' if kind == 'Main' else f'Supplementary-Figure-{n}'
    fig.savefig(d / f'{stem}.png', dpi=DPI); fig.savefig(d / f'{stem}.pdf', dpi=DPI); plt.close(fig)
    # legend: figure title, then one sentence per block from its recorded headline and subtitle
    parts = []
    for m in metas:
        last = chr(ord(m['first']) + m['n_letters'] - 1)
        lab = m['first'] if m['first'] == last else f"{m['first']}–{last}"
        k = ord(m['first']) - ord('a') if not m['add_letter'] else 0
        t = (m['title'] or '').rstrip('.'); s = shift_letters((m['sub'] or '').replace('\n', ' ').strip(), k)
        notes = ' '.join(shift_letters(n, k) for n in m.get('notes', []))
        parts.append(f"**{lab},** {t}. {s} {notes}".strip())
    label = 'Figure' if kind == 'Main' else 'Supplementary Fig.'
    return f"**{label} {n} | {title}.** " + ' '.join(parts), round(H * scale, 1), [m['slug'] for m in metas]


if __name__ == '__main__':
    do_render = '--no-render' not in sys.argv
    if '--frozen-data' in sys.argv:
        os.environ['GAL_FROZEN_DATA'] = '1'
    leg, log = ['# Figure legends\n'], {}
    for kind, figs in [('Main', MAIN), ('Supplementary', EXTENDED)]:
        if kind == 'Supplementary':
            leg.append('\n# Supplementary figure legends\n')
        for i, (title, slugs) in enumerate(figs, 1):
            text, h, used = assemble(kind, i, title, slugs, do_render)
            leg.append(text + '\n'); log[f'{kind} {i}'] = {'width_in': W_IN, 'height_mm': h, 'dpi': DPI, 'blocks': used}
            print(kind, i, title, h, 'mm', flush=True)
    (OUT / 'Legends.md').write_text('\n'.join(leg))
    (OUT / 'Layout.json').write_text(json.dumps(log, indent=1))
