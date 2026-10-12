"""What is inside the deposited results archives (members listed by archive_contents.py)?

Every member is classified with the same rules as top-level files (layer_rules.py), by name
only: members have no submitter category, and the parent's is not inherited, so a bare .xml or
.out inside a results zip is not counted as tool output (conservative). Members are also matched
against each registry tool's native-file patterns. Directory entries and macOS resource forks
(__MACOSX/, ._*) are dropped.

Questions:
  1. Content of each archive: the layer holding most of its uncompressed bytes, and whether it
     holds a standard file, known tool results, raw data, spectra or only tables.
  2. Check of the name-based layer rules: archives we classified from their names as tool
     results (e.g. txt.zip, maxquant*.zip) — do they contain tool results?
  3. Dataset level, before vs after opening archives: share with standard files, with
     recognisable results (standard or known tool output), and with results only in ad hoc
     tables; by release year. The primary outcome (`standards`, top-level files) is unchanged;
     'standard inside an archive' is reported as a labelled sensitivity.
Denominators: every candidate archive (Targets.csv.gz) is counted as listed / not listable
(tar.gz, tgz, rar, 7z, tar) / failed, overall and per dataset.

The listing holds ~36 million members; it is converted to Parquet and all work on members is done
in DuckDB (memory-capped).
"""
import json
import numpy as np, pandas as pd, duckdb
from common import R, O
from layer_rules import RULES
from software_registry import TOOLS

AC = O / 'Archive-contents'
LIST, MEMB = AC / 'Listings.jsonl', AC / 'Members.parquet'
L = pd.read_parquet(O / 'Layers/Dataset-layers.parquet').set_index('accession')
T = pd.read_csv(AC / 'Targets.csv.gz')

# ---- stream listings -> archive table + member parquet --------------------------------------
if not MEMB.exists() or (LIST.exists() and MEMB.stat().st_mtime < LIST.stat().st_mtime):
    # DuckDB's native JSON reader: ~15 s for the 2.4 GB listing (a Python json loop took minutes)
    cx = duckdb.connect(); cx.execute("set memory_limit='10GB'; set threads=4; set preserve_insertion_order=false")
    src = f"""read_ndjson('{LIST}', columns={{'url': 'VARCHAR', 'accession': 'VARCHAR', 'file_name': 'VARCHAR', 'layer': 'VARCHAR',
              'cat': 'VARCHAR', 'size': 'BIGINT', 'status': 'VARCHAR', 'n_members': 'BIGINT', 'members': 'JSON'}},
              maximum_object_size=67108864)"""
    cx.execute(f"""copy (select url, m[1] member, m[2]::BIGINT usize from (
                 select url, unnest(json_transform(members, '[["VARCHAR"]]')) m from {src} where status in ('ok', 'partial')))
               to '{MEMB}' (format parquet, compression zstd)""")
    cx.execute(f"""copy (select any_value(accession) accession, url, any_value(file_name) file_name, any_value(layer) layer,
                 any_value(cat) cat, any_value(size) size, arg_max(status, status in ('ok', 'partial')) status, max(n_members) n_members
                 from {src} group by url) to '{AC / 'Archives.csv.gz'}' (header, compression gzip)""")
A = pd.read_csv(AC / 'Archives.csv.gz')
A['listed'] = A.status.isin(['ok', 'partial'])   # partial: member count differs from the header; members kept

con = duckdb.connect(); con.execute("set memory_limit='4GB'; set threads=2")
con.execute(f"set temp_directory='{AC / '.duckdb_tmp'}'")
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
tool_cols = ', '.join(f"regexp_matches(b, '{'|'.join(f'(?:{p})' for p in t['native']).replace(chr(39), chr(39) * 2)}') as \"tool_{k}\""
                      for k, t in TOOLS.items())
con.execute(f"""
create table m as
with raw as (select url, member, usize, lower(member) n,
                    lower(regexp_replace(member, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b, '' cat
             from read_parquet('{MEMB}')
             where not ends_with(member, '/') and not regexp_matches(member, '(^|/)__MACOSX/|(^|/)\\._'))
select url, usize, {case} inner_layer, regexp_matches(n, '\\.(zip|tgz|tar|rar|7z)$|\\.tar\\.gz$') nested, {tool_cols}
from raw""")
TOOLC = [f'tool_{k}' for k in TOOLS]
n_members = con.execute('select count(*) from m').fetchone()[0]

# ---- 1. per archive --------------------------------------------------------------------------
any_tools = ', '.join(f'bool_or("{c}") "{c}"' for c in TOOLC)
per = con.execute(f"""
with lb as (select url, inner_layer, sum(usize) s from m group by all),
dom as (select url, arg_max(inner_layer, s) dominant_layer from lb group by url)
select m.url, count(*) members, sum(usize) usize,
  bool_or(inner_layer = 'standard_results') has_standard, bool_or(inner_layer = 'tool_results') has_tool_results,
  bool_or(inner_layer = 'vendor_raw') has_raw, bool_or(inner_layer = 'open_spectra') has_spectra,
  bool_or(inner_layer = 'unidentified_tables') has_tables, bool_or(inner_layer = 'sequence_db') has_fasta,
  bool_or(inner_layer = 'sample_metadata') has_sdrf, bool_or(nested) has_nested_archive,
  any_value(dominant_layer) dominant_layer, {any_tools}
from m join dom using (url) group by m.url""").df()
per['tools_found'] = per[TOOLC].apply(lambda r: ', '.join(c[5:] for c in TOOLC if r[c]), axis=1)
per['content'] = np.select(
    [per.has_standard, per.has_tool_results, per.has_raw | per.has_spectra, per.has_tables, per.has_nested_archive],
    ['standard file', 'known tool results', 'raw data or spectra (no results)', 'tables only', 'nested archive only'], 'other')
A = A.merge(per, on='url', how='left')
for c in ['has_standard', 'has_tool_results', 'has_raw', 'has_spectra', 'has_tables', 'has_fasta', 'has_sdrf'] + TOOLC:
    A[c] = A[c].fillna(False).astype(bool)
A.drop(columns=TOOLC).to_csv(AC / 'Archive-content.csv.gz', index=False)

# denominators over all candidate archives
T['listed'] = T.url.isin(A.url[A.listed])
T['why'] = np.where(T.listed, 'listed', np.where(T.fmt != 'zip', 'format not listable (' + T.fmt + ')',
           T.url.map(A.set_index('url').status).fillna('not attempted (empty file)')))
den = T.groupby('why').agg(archives=('url', 'size'), datasets=('accession', 'nunique'),
                           tb=('file_size_bytes', lambda x: x.sum() / 1e12)).round(3)
den.to_csv(AC / 'Denominators.csv')

ok = A[A.listed]
content = ok.groupby(['layer', 'content']).size().unstack(fill_value=0)
content.to_csv(AC / 'Content-by-parent-layer.csv')
content_pct = content.div(content.sum(1), axis=0).mul(100).round(1)
dominant = ok.groupby(['layer', 'dominant_layer']).usize.sum().unstack(fill_value=0)
dominant.div(dominant.sum(1), axis=0).mul(100).round(1).to_csv(AC / 'Uncompressed-bytes-by-dominant-layer.csv')

# ---- 2. name-based rule check ---------------------------------------------------------------
nb = ok[ok.layer == 'tool_results']
rule_check = {'archives_named_as_tool_results': int(len(nb)),
              'contain_tool_results_or_standard_pct': round(100 * (nb.has_tool_results | nb.has_standard).mean(), 2),
              'contain_tables_only_pct': round(100 * (nb.content == 'tables only').mean(), 2),
              'contain_raw_or_spectra_without_results_pct': round(100 * (nb.content == 'raw data or spectra (no results)').mean(), 2)}

# ---- 3. dataset level before / after --------------------------------------------------------
dsa = ok.groupby('accession').agg(inner_standard=('has_standard', 'any'), inner_tool=('has_tool_results', 'any'),
                                  inner_tables=('has_tables', 'any'), inner_raw=('has_raw', 'any'),
                                  inner_fasta=('has_fasta', 'any'), inner_sdrf=('has_sdrf', 'any'))
unread = set(T.accession[~T.listed])
D = L[['release_year', 'dia_expanded', 'standards', 'has_standard_results', 'has_tool_results', 'has_unidentified_tables',
       'has_archived_results', 'has_vendor_raw', 'has_open_spectra', 'has_sequence_db', 'has_sample_metadata']].join(dsa)
for c in dsa.columns:
    D[c] = D[c].fillna(False).astype(bool)
D['unread_archive'] = D.index.isin(unread)
D['standard_after'] = D.has_standard_results | D.inner_standard
D['recognisable_before'] = D.has_standard_results | D.has_tool_results
D['recognisable_after'] = D.standard_after | D.has_tool_results | D.inner_tool
D['adhoc_only_before'] = ~D.recognisable_before & (D.has_unidentified_tables | D.has_archived_results)
D['adhoc_only_after'] = ~D.recognisable_after & (D.has_unidentified_tables | D.inner_tables | (D.has_archived_results & D.unread_archive))
D['fasta_after'] = D.has_sequence_db | D.inner_fasta
D['sdrf_after'] = D.has_sample_metadata | D.inner_sdrf
D.to_csv(AC / 'Dataset-before-after.csv.gz')
cols = ['has_standard_results', 'standard_after', 'recognisable_before', 'recognisable_after', 'adhoc_only_before',
        'adhoc_only_after', 'has_sequence_db', 'fasta_after', 'has_sample_metadata', 'sdrf_after', 'unread_archive']
by_year = D.groupby('release_year')[cols].mean().mul(100).round(1)
by_year['datasets'] = D.groupby('release_year').size()
by_year.to_csv(AC / 'Before-after-by-year.csv')
by_mode = D[D.release_year >= 2019].groupby('dia_expanded')[cols].mean().mul(100).round(1)
by_mode['datasets'] = D[D.release_year >= 2019].groupby('dia_expanded').size()
by_mode.index = by_mode.index.map({True: 'DIA', False: 'not DIA'}); by_mode.to_csv(AC / 'Before-after-by-mode-2019-2025.csv')

# which tools' outputs sit inside archives, 2019-2025 datasets
ok19 = ok[ok.accession.map(L.release_year) >= 2019]
tool_inside = ok19.groupby('accession')[TOOLC].any().sum().sort_values(ascending=False).rename(lambda c: c[5:])
tool_inside.to_csv(AC / 'Tool-outputs-inside-archives-2019-2025.csv', header=['datasets'])

summary = {'archives_listed': int(len(ok)), 'members': int(n_members), 'denominators': den.reset_index().to_dict('records'),
           'content_of_listed_archives_pct': ok.content.value_counts(normalize=True).mul(100).round(1).to_dict(),
           'name_based_rule_check': rule_check,
           'datasets_with_standard_only_inside_archives': int((D.inner_standard & ~D.has_standard_results).sum()),
           'datasets_2019_2025': {k: round(float(v), 1) for k, v in D[D.release_year >= 2019][cols].mean().mul(100).items()}}
(AC / 'Summary.json').write_text(json.dumps(summary, indent=2, default=str))
pd.set_option('display.width', 250)
print(json.dumps(summary, indent=1, default=str)); print(content_pct.to_string()); print(by_year.to_string()); print(by_mode.to_string())
print(tool_inside.head(15).to_string())
