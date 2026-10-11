"""Automated check of the file-layer classification against two independent references.

1. PRIDE's own rules: every filename pattern and raw-format extension in
   reference/pride-submission-validation.yaml (v2.2.2) is applied to each file name, giving
   PRIDE's category (RAW, ANALYSIS, STANDARD, PEAK, SPECTRAL_LIBRARY, FASTA, METADATA, OTHER).
   Specific names (e.g. evidence.txt, report.parquet, *.msf) take precedence over bare generic
   extensions (.tsv/.csv/.txt/.xlsx/.json/.xml/.parquet), which are ambiguous by PRIDE's own
   account and are reported separately.
2. The category each submitter chose for the file at upload (RAW, PEAK, SEARCH, RESULT, OTHER).

Agreement is computed per file and per dataset-layer flag. Disagreements are written out for
inspection. Neither reference is ground truth; agreement between independent methods is the
evidence that the classification is not an artefact of our own rules.
"""
import re, json
import pandas as pd, yaml
from common import R, O
from layer_rules import RULES
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Audit'; OUT.mkdir(exist_ok=True)
Y = yaml.safe_load(open(R / 'work/stack-and-reuse/reference/pride-submission-validation.yaml'))
GENERIC = {'tsv', 'csv', 'txt', 'xlsx', 'json', 'xml', 'parquet', 'properties', 'py', 'log.txt', 'pdf'}


def pat2re(p):
    p = p.strip().lower()
    rx = re.escape(p).replace(re.escape('{prefix}'), '.*').replace(r'\*', '.*')
    return '^' + rx + '$'


specific, generic = [], []
entries = list(Y['common_files'].values()) + [f for t in Y['tools'] for f in t.get('files', [])]
for f in entries:
    for p in (f.get('filename_patterns') or []):
        for q in p.split(','):
            q = q.strip()
            ext = q.lower().replace('{prefix}', '').lstrip('.')
            (generic if (q.lower().startswith('{prefix}.') and ext in GENERIC) else specific).append((pat2re(q), f['category']))
for rf in Y['raw_formats']:
    for e in rf['extensions']:
        specific.append(('^.*' + re.escape(e.lower()) + '$', rf['category']))
specific = list(dict.fromkeys(specific)); generic = list(dict.fromkeys(generic))
spec_re = [(re.compile(p), c) for p, c in specific]
gen_re = [(re.compile(p), c) for p, c in generic]


def pride_category(name):
    cats = {c for r, c in spec_re if r.match(name)}
    if cats:
        return '|'.join(sorted(cats)), 'specific'
    cats = {c for r, c in gen_re if r.match(name)}
    return ('|'.join(sorted(cats)), 'generic') if cats else (None, None)


case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
c = connect()
f = c.execute(f"""
with f as (select project_accession a, lower(file_name) n,
  lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
  coalesce(file_category_value, '') cat
  from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
  and project_accession in (select accession from study))
select a accession, n, cat, {case} layer from f""").df()
# PRIDE patterns are written against the uploaded name (e.g. sample.d.zip, report.tsv); try the
# full name first and then the name without its compression suffix.
u = pd.Series(f.n.unique())
pc = {n: pride_category(n) for n in u}
for n, v in list(pc.items()):
    if v[0] is None:
        stripped = re.sub(r'\.(gz|zip|bz2|7z|tar|xz|tgz)$', '', n)
        if stripped != n:
            pc[n] = pride_category(stripped)
f['pride_cat'] = f.n.map(lambda n: pc[n][0]); f['pride_match'] = f.n.map(lambda n: pc[n][1])
EXPECT = {'vendor_raw': {'RAW'}, 'open_spectra': {'PEAK', 'RAW'}, 'standard_results': {'STANDARD'},
          'tool_results': {'ANALYSIS', 'OTHER', 'PEAK'}, 'spectral_library': {'SPECTRAL_LIBRARY'},
          'sequence_db': {'FASTA'}, 'sample_metadata': {'METADATA'}, 'archived_results': {'ANALYSIS'},
          'unidentified_tables': {'ANALYSIS', 'OTHER', 'METADATA'}, 'documents': {'OTHER'}, 'bookkeeping': {'OTHER'},
          'other': {'OTHER'}}
f['agree'] = [bool(set(pcat.split('|')) & EXPECT[l]) if isinstance(pcat, str) else None for pcat, l in zip(f.pride_cat, f.layer)]
res = {}
for m in ['specific', 'generic']:
    x = f[f.pride_match == m]
    res[m] = {'files': int(len(x)), 'agreement_pct': round(100 * x.agree.mean(), 2),
              'by_layer': x.groupby('layer').agree.agg(['size', 'mean']).round(4).to_dict('index')}
res['unmatched_by_pride_rules'] = int(f.pride_cat.isna().sum())
res['unmatched_by_layer'] = f[f.pride_cat.isna()].layer.value_counts().to_dict()
# submitter-chosen category vs our layer
SUB = {'vendor_raw': {'RAW'}, 'open_spectra': {'PEAK', 'RAW'}, 'standard_results': {'RESULT', 'SEARCH'},
       'tool_results': {'SEARCH', 'RESULT', 'OTHER', 'PEAK'}, 'archived_results': {'SEARCH', 'RESULT'}}
s = f[f.layer.isin(SUB)]
s_ag = [cat in SUB[l] for cat, l in zip(s.cat, s.layer)]
res['submitter_category_agreement'] = {'files': int(len(s)), 'agreement_pct': round(100 * pd.Series(s_ag).mean(), 2),
    'by_layer': pd.DataFrame({'layer': s.layer.values, 'ok': s_ag}).groupby('layer').ok.mean().round(4).to_dict()}
(OUT / 'Automated-validation.json').write_text(json.dumps(res, indent=2))
f[(f.agree == False) & (f.pride_match == 'specific')].groupby(['layer', 'pride_cat']).n.agg(['size', lambda x: list(x.head(5))]).reset_index().to_csv(OUT / 'Disagreements-specific.csv', index=False)
print(json.dumps({k: v for k, v in res.items() if k != 'specific'} | {'specific': {k: res['specific'][k] for k in ['files', 'agreement_pct']}}, indent=1))
print(pd.DataFrame(res['specific']['by_layer']).T.to_string())
