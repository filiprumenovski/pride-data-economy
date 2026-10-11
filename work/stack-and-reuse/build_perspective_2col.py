"""Typeset Perspective-draft.md as a two-column journal-style A4 PDF (headless Chromium, CSS columns).

Layout: full-width title block with the abstract as a bold standfirst; two justified serif columns;
figures span both columns and are placed after the section that discusses them, with their legends;
boxes are shaded panels in the column flow; end matter and numbered references follow the body.

  python build_perspective_2col.py  ->  outputs/Stack-and-reuse/Perspective-2col.html / .pdf
"""
import base64, re, subprocess
import markdown
from common import O

FIG = O / 'Figures'
CHROME = '/opt/pw-browsers/chromium-1194/chrome-linux/chrome'
FIGS = ['Figure-1-anatomy-of-a-deposit', 'Figure-2-instruments', 'Figure-3-software-and-identification-files',
        'Figure-4-dia-and-identification-files', 'Figure-5-same-groups-different-software', 'Figure-6-world-maps']


def blank_before_lists(t):
    out, prev = [], ''
    for line in t.splitlines():
        item = bool(re.match(r'\s*([-*]|\d+\.)\s', line)); prev_item = bool(re.match(r'\s*([-*]|\d+\.)\s', prev))
        if item and prev.strip() and not prev_item:
            out.append('')
        out.append(line); prev = line
    return '\n'.join(out)

MD = O / 'Perspective-draft.md'
text = blank_before_lists(MD.read_text())
md = lambda s: markdown.markdown(s, extensions=['tables'])

title = re.search(r'^# (.+)$', text, re.M).group(1)
secs = {m.group(1).strip(): m.group(2).strip() for m in re.finditer(r'^## (.+?)\n(.*?)(?=^## |\Z)', text, re.M | re.S)}
abstract = secs['Abstract']
body_keys = [k for k in secs if re.match(r'\d\.', k)]
legend_blocks = [l for l in re.split(r'\n(?=\*\*Figure \d)', secs['Figure legends']) if l.startswith('**Figure')]
LEG = dict(zip(FIGS, legend_blocks))
# figure placement: after the section that discusses it
AFTER = {'2.': [FIGS[0]], '3.': [FIGS[1]], '4.': [FIGS[2], FIGS[3], FIGS[4]], '5.': [FIGS[5]]}
BOX_AFTER = {'1.': 'Box 1. Data and definitions', '7.': 'Box 2. What these data cannot say'}


def figure(name):
    b64 = base64.b64encode((FIG / f'{name}.png').read_bytes()).decode()
    lines = [l.strip() for l in LEG[name].strip().splitlines() if l.strip()]
    head = re.sub(r'^\*\*(Figure \d) \| ([^*]+)\*\*$', r'**\1 | \2**', lines[0])
    panels = ' '.join(re.sub(r'^-\s+', '', l) for l in lines[1:])   # run-in panel descriptions, journal style
    return (f'<figure class="span"><img src="data:image/png;base64,{b64}"/>'
            f'<figcaption>{md(head + " " + panels)}</figcaption></figure>')


def strip_trailing_rule(s):
    return re.sub(r'\n-{3,}\s*$', '', s.strip())


html_body = ''
for k in body_keys:
    num = k.split(' ', 1)[0]
    html_body += f'<h2>{k}</h2>' + md(strip_trailing_rule(secs[k]))
    if num in BOX_AFTER:
        bk = BOX_AFTER[num]
        html_body += f'<div class="box"><h3>{bk}</h3>{md(strip_trailing_rule(secs[bk]))}</div>'
    for f in AFTER.get(num, []):
        html_body += figure(f)

refs = secs[[k for k in secs if k.startswith('References')][0]]
refs = re.sub(r'^Status:.*$', '', refs, flags=re.M)
end = (f'<h2 class="end">Data and code availability</h2>{md(secs["Data and code availability"])}'
       f'<h2 class="end">Supplementary information</h2>{md(secs["Supplementary information"])}'
       f'<h2 class="end">References</h2><div class="refs">{md(refs)}</div>'
       '<p class="note">References marked (S) were checked only against a search snippet or abstract and must be verified against the full text.</p>')
CSS = """
@page { size: A4; margin: 16mm 15mm 17mm 15mm; }
* { box-sizing: border-box; }
body { margin: 0; font-family: 'Liberation Serif', 'DejaVu Serif', serif; font-size: 9.1pt; line-height: 1.36; color: #111;
       hyphens: auto; -webkit-hyphens: auto; }
.kicker { font-family: 'Liberation Sans', sans-serif; font-size: 7.6pt; letter-spacing: .14em; text-transform: uppercase;
          color: #2f5f8f; font-weight: bold; border-top: 2.2pt solid #2f5f8f; padding-top: 4pt; margin-bottom: 6pt; }
h1 { font-family: 'Liberation Sans', sans-serif; font-size: 21pt; line-height: 1.12; margin: 0 0 7pt; font-weight: bold; letter-spacing: -.01em; }
.authors { font-family: 'Liberation Sans', sans-serif; font-size: 8.6pt; color: #333; margin-bottom: 9pt; }
.standfirst { font-family: 'Liberation Sans', sans-serif; font-size: 9.2pt; line-height: 1.38; font-weight: bold; color: #1a1a1a;
              border-bottom: .6pt solid #999; padding-bottom: 9pt; margin-bottom: 10pt; }
.standfirst p { margin: 0 0 5pt; } .standfirst ul { margin: 2pt 0 5pt; padding-left: 12pt; } .standfirst li { margin: 0 0 1.5pt; }
.cols { column-count: 2; column-gap: 6.5mm; column-fill: balance; }
.cols p { margin: 0 0 5pt; text-align: justify; }
.cols ul, .cols ol { margin: 0 0 5pt; padding-left: 11pt; } .cols li { margin: 0 0 1.5pt; text-align: left; }
h2 { font-family: 'Liberation Sans', sans-serif; font-size: 9.8pt; font-weight: bold; margin: 9pt 0 3pt; color: #111; break-after: avoid; }
h2.end { font-size: 8.4pt; margin-top: 8pt; }
table { width: 100%; border-collapse: collapse; margin: 3pt 0 7pt; font-family: 'Liberation Sans', sans-serif; font-size: 6.9pt;
        line-height: 1.25; break-inside: avoid; }
th { text-align: left; border-top: .8pt solid #222; border-bottom: .5pt solid #222; padding: 2pt 4pt 2pt 0; font-weight: bold; }
td { border-bottom: .3pt solid #c8c8c8; padding: 1.8pt 4pt 1.8pt 0; vertical-align: top; }
tr:last-child td { border-bottom: .8pt solid #222; }
code { font-family: 'DejaVu Sans Mono', monospace; font-size: 7pt; }
.box { background: #eef2f7; border-top: 1.4pt solid #2f5f8f; padding: 6pt 7pt 3pt; margin: 6pt 0 8pt;
       font-family: 'Liberation Sans', sans-serif; font-size: 7.5pt; line-height: 1.33; }
.box h3 { font-size: 8.2pt; margin: 0 0 4pt; color: #2f5f8f; } .box p { text-align: left; margin: 0 0 4pt; } .box ul { padding-left: 10pt; }
figure.span { column-span: all; break-inside: avoid; margin: 6pt 0 10pt; padding: 0; }
figure.span img { width: 100%; display: block; }
figcaption { font-family: 'Liberation Sans', sans-serif; font-size: 7.2pt; line-height: 1.3; margin-top: 3pt; color: #222; }
figcaption p { margin: 0 0 1.5pt; text-align: left; } figcaption ul { margin: 0; padding-left: 0; list-style: none; }
figcaption li { display: inline; margin: 0 2pt 0 0; }
.refs { font-size: 7.4pt; line-height: 1.3; } .refs ol { padding-left: 13pt; }
.note { font-size: 7pt; color: #555; font-style: italic; }
.cols > h2:first-child { margin-top: 0; }
"""
html = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{title}</title><style>{CSS}</style></head><body>
<div class="kicker">Perspective</div><h1>{title}</h1>
<div class="authors">[Authors and affiliations] · Draft, October 2026 · Not submitted</div>
<div class="standfirst">{md(abstract)}</div>
<div class="cols">{html_body}{end}</div></body></html>"""
out_html = O / 'Perspective-2col.html'; out_pdf = O / 'Perspective-2col.pdf'
out_html.write_text(html)
subprocess.run([CHROME, '--headless', '--no-sandbox', '--disable-gpu', '--no-pdf-header-footer',
                f'--print-to-pdf={out_pdf}', out_html.as_uri()], check=True, capture_output=True, timeout=180)
print(out_pdf, out_pdf.stat().st_size // 1024, 'KB')
