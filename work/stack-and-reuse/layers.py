"""What does a public proteomics dataset contain? File-level layer classification.

Every listed file outside PRIDE's /generated/ folders is assigned one layer from its name
(and PRIDE file category where the name is ambiguous). Rules are ordered; the first match wins.
Datasets get one boolean and one byte total per layer. Known tool outputs are matched by
their documented file names; other tables are kept separate as "unidentified tables" so the
result layer can be reported strictly (known outputs only) and loosely (+ tables).
"""
import json
import pandas as pd
from common import R, O
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Layers'; OUT.mkdir(exist_ok=True)
c = connect()
# lower-cased name with one trailing compression/archive suffix removed (b); full lower name (n)
from layer_rules import RULES
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
files = c.execute(f"""
with f as (
  select project_accession a, lower(file_name) n,
         lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
         coalesce(file_category_value, '') cat, coalesce(file_size_bytes, 0) sz, checksum
  from pride_files
  where project_accession in (select accession from study)
    and coalesce(primary_public_location_value, '') not like '%/generated/%')
select a accession, {case} layer, count(*) files, sum(sz) bytes from f group by 1, 2""").df()
w = files.pivot_table(index='accession', columns='layer', values='files', fill_value=0)
b = files.pivot_table(index='accession', columns='layer', values='bytes', fill_value=0)
layers = [l for l, _ in RULES] + ['other']
w = w.reindex(columns=layers, fill_value=0); b = b.reindex(columns=layers, fill_value=0)
d = pd.read_parquet(O / 'Deposit-features.parquet').set_index('accession')
ds = (w > 0).add_prefix('has_').join(b.add_prefix('bytes_')).join(w.add_prefix('n_'))
ds = d[['release_year', 'country', 'pi_group', 'dia_expanded', 'vendor_single', 'standards', 'submission_type']].join(ds, how='left')
for l in layers:
    ds['has_' + l] = ds['has_' + l].fillna(False).astype(bool)
# check: standard_results layer must reproduce the corrected `standards` flag exactly
mismatch = int((ds.has_standard_results != ds.standards.astype(bool)).sum())
ds['any_results_strict'] = ds.has_standard_results | ds.has_tool_results
ds['any_results_loose'] = ds.any_results_strict | ds.has_unidentified_tables | ds.has_archived_results
ds.reset_index().to_parquet(OUT / 'Dataset-layers.parquet', index=False)
summary = {'datasets': len(ds), 'files_classified': int(files.files.sum()), 'bytes_classified_PB': float(files.bytes.sum() / 1e15),
           'standards_flag_mismatches': mismatch,
           'files_by_layer': files.groupby('layer').files.sum().sort_values(ascending=False).astype(int).to_dict(),
           'TB_by_layer': (files.groupby('layer').bytes.sum() / 1e12).round(1).sort_values(ascending=False).to_dict()}
(OUT / 'Layer-summary.json').write_text(json.dumps(summary, indent=2)); print(json.dumps(summary, indent=2))
