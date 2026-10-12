"""Typeset Paper.md as a two-column A4 article (headless Chromium, CSS columns).

An unnumbered graphical overview follows the abstract. Article text flows in two columns;
main and supplementary figures follow the references.
Figures and their legends stay together, with multiple figures per page when they fit.
The original lossless 600-dpi PNG masters are embedded without downsampling.

  python build_paper_2col.py  ->  outputs/Stack-and-reuse/Paper-2col.html / .pdf
"""
import os, re, subprocess, json
from html import escape
from pathlib import Path
import markdown
from PIL import Image
from common import O

# Reuse the house CSS and helpers from build_perspective_2col.py without running that script.
_src = (Path(__file__).resolve().parent / 'build_perspective_2col.py').read_text()
CSS = re.search(r'CSS = """(.*?)"""', _src, re.S).group(1)
CHROME = os.environ.get('CHROME') or re.search(r"CHROME = '([^']+)'", _src).group(1)
exec(re.search(r'(def blank_before_lists\(t\):.*?return .*?\n)', _src, re.S).group(1))

J = O / 'Journal-figures'
md = lambda s: markdown.markdown(s, extensions=['tables'])
SUP = str.maketrans('0123456789,', '⁰¹²³⁴⁵⁶⁷⁸⁹˒')


def number_citations(t):
    """[Key] citations -> superscript numbers in order of first appearance; reference list renumbered."""
    body, refs = t.split('## References', 1)
    order = []
    for k in re.findall(r'\[([A-Za-z][A-Za-z0-9]+)\]', body):
        if k not in order:
            order.append(k)
    num = {k: i + 1 for i, k in enumerate(order)}
    def rep(m):
        keys = re.findall(r'\[([A-Za-z][A-Za-z0-9]+)\]', m.group(0))
        return '<sup>' + ','.join(str(num[k]) for k in keys) + '</sup>'
    body = re.sub(r'(?:\s?\[[A-Za-z][A-Za-z0-9]+\])+', rep, body)
    entries = dict(re.findall(r'^\[([A-Za-z][A-Za-z0-9]+)\] (.+)$', refs, re.M))
    missing = [k for k in order if k not in entries]; unused = [k for k in entries if k not in num]
    assert not missing, f'citations without reference: {missing}'
    if unused:
        print('unused references dropped:', unused)
    head = refs.split('\n', 1)[0]
    tail = re.sub(r'^\[[A-Za-z][A-Za-z0-9]+\] .+\n?', '', refs.split('\n', 1)[1], flags=re.M).strip()
    lst = '\n'.join(f'{num[k]}. {entries[k]}' for k in order)
    return body + '## References' + head + '\n\n' + lst + '\n\n' + tail + '\n'


text = blank_before_lists(number_citations((O / 'Paper.md').read_text()))
title = re.search(r'^# (.+)$', text, re.M).group(1)
secs = {m.group(1).strip(): m.group(2).strip() for m in re.finditer(r'^## (.+?)\n(.*?)(?=^## |\Z)', text, re.M | re.S)}
legends = (J / 'Legends.md').read_text()
LEG = {}
for m in re.finditer(r'^\*\*(Figure|Supplementary Fig\.) (\d+) \|(.*?)$', legends, re.M):
    LEG[(m.group(1), int(m.group(2)))] = m.group(0)


def figure(kind, n):
    path = J / ('Main' if kind == 'Figure' else 'Supplementary') / (f'Figure-{n}.png' if kind == 'Figure' else f'Supplementary-Figure-{n}.png')
    with Image.open(path) as im:
        assert im.width == 4320, f'{path}: expected 7.20 inches at 600 dpi (4320 pixels)'
    return f'<figure class="span"><img src="{path.as_uri()}"/><figcaption>{md(LEG[(kind, n)])}</figcaption></figure>'


body = f'<h2>Introduction</h2>{md(secs["Introduction"])}<h2>Results</h2>'
subs = re.findall(r'^### (.+?)\n(.*?)(?=^### |\Z)', secs['Results'], re.M | re.S)
for h, b in subs:
    body += f'<h3 class="sub">{h.strip()}</h3>' + md(b.strip())
body += f'<h2>Discussion</h2>{md(secs["Discussion"])}'
body += f'<div class="box"><h3>Box 1. Data, definitions and limits</h3>{md(secs["Box 1. Data, definitions and limits"])}</div>'
refs_key = [k for k in secs if k.startswith('References')][0]
refs = re.sub(r'\n\(S\):.*$', '', secs[refs_key], flags=re.S)
end = (f'<h2 class="end">Methods</h2>{md(secs["Methods"])}'
       f'<h2 class="end">Supplementary information</h2>{md(secs["Supplementary information"])}'
       f'<h2 class="end">Data and code availability</h2>{md(secs["Data and code availability"])}'
       f'<h2 class="end">Author declarations</h2>{md(secs["Author declarations"])}'
       f'<h2 class="end">References</h2><div class="refs">{md(refs)}</div>')
supp_methods = f'<section class="figures"><h2>Supplementary methods</h2><div class="cols">{md(secs["Supplementary methods"])}</div></section>'
reference_guide = f'<section class="figures reference-guide"><h2>Supplementary reference guide</h2>{md(secs["Supplementary reference guide"])}</section>'
n_ext = max(n for (k, n) in LEG if k != 'Figure')
main = '<section class="figures"><h2>Main figures</h2>' + ''.join(figure('Figure', i) for i in range(1, 6)) + '</section>'
ext = '<section class="figures"><h2>Supplementary figures</h2>' + ''.join(figure('Supplementary Fig.', i) for i in range(1, n_ext + 1)) + '</section>'
CSS2 = CSS + """
@page { size: A4; margin: 12mm 13.5mm; }
body { line-height: 1.32; }
.authors, a { font-variant-ligatures: none; }
a { color: inherit; text-decoration: none; }
.cols p { margin-bottom: 4pt; }
.figures { break-before: page; }
.figures h2 { margin: 0 0 4pt; font-size: 10pt; break-after: avoid; }
figure.span { margin: 0 0 4pt; break-inside: avoid; }
figure.span img { width: 7.2in; max-width: 100%; }
.cols { orphans: 3; widows: 3; }
.refs { overflow-wrap: anywhere; }
.reference-guide { font-size: 8pt; }
.reference-guide h3 { font-size: 9pt; margin: 8pt 0 3pt; break-after: avoid; }
.reference-guide table { font-size: 7.5pt; break-inside: auto; }
.reference-guide tr { break-inside: avoid; }
.reference-guide sup { font-size: inherit; vertical-align: baseline; }
.reference-guide thead { display: table-header-group; }
h3.sub { font-family: 'Liberation Sans', sans-serif; font-size: 8.8pt; font-weight: bold; margin: 7pt 0 2.5pt; break-after: avoid; }
figure.span figcaption { font-size: 7pt; line-height: 1.2; }
.standfirst ol { margin: 2pt 0 5pt; padding-left: 14pt; }
.graphical-abstract { margin: 0 0 10pt; break-inside: avoid; }
.graphical-abstract img { display: block; width: 7.2in; max-width: 100%; height: auto; }
.graphical-abstract p { margin: 0; }
.graphical-abstract p:last-child { margin-top: 3pt; font-family: sans-serif; font-size: 6pt; color: #555; }
"""
metadata = json.loads((O / 'Submission-metadata.json').read_text())
author_line = escape(metadata['author'] + ' · ' + metadata['affiliation'])
if metadata.get('corresponding_email'):
    author_line += '<br>Correspondence: ' + escape(metadata['corresponding_email'])
html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{title}</title><style>{CSS2}</style></head><body>
<div class="kicker">Analysis</div><h1>{title}</h1>
<div class="authors">{author_line}<br>Pre-submission review · October 2026</div>
<div class="standfirst">{md(secs['Abstract'])}</div>
<figure class="graphical-abstract">{md(secs['Graphical abstract'])}</figure>
<div class="cols">{body}{end}</div>{main}{supp_methods}{reference_guide}{ext}</body></html>"""
out_html = O / 'Paper-2col.html'; out_pdf = O / 'Paper-2col.pdf'
out_html.write_text(html)
subprocess.run([CHROME, '--headless', '--no-sandbox', '--disable-gpu', '--no-pdf-header-footer', '--allow-file-access-from-files',
                f'--print-to-pdf={out_pdf}', out_html.as_uri()], check=True, capture_output=True, timeout=300)
print(out_pdf, out_pdf.stat().st_size // 1024, 'KB')
