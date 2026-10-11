"""Deposit-level feature table for the analytical-stack / reuse study.

New experiment (October 2026). Reads the frozen databases read-only and writes only to
outputs/Stack-and-reuse/. The primary standard-result outcome is the corrected
`standards` (= nongenerated_expanded) flag from the tracked cohort; legacy inclusive flags
are carried only under their explicit legacy name and are never used as an outcome.

File-derived signals classify *listed file names* outside repository generated/ paths.
They are evidence that an output of a given tool was deposited, not proof of the tool's
use, and their absence is "not detected", not verified absence.
"""
from pathlib import Path
import sys, json
import pandas as pd
R = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R / 'scripts'))
from query import connect

O = R / 'outputs/Stack-and-reuse'; O.mkdir(parents=True, exist_ok=True)
c = connect()

# --- file-name evidence, non-generated paths only -------------------------------------
files = c.execute(r"""
with f as (
  select project_accession a,
         lower(file_name) n,
         lower(regexp_replace(file_name, '\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
         coalesce(file_size_bytes,0) sz, coalesce(total_downloads,0) dl, file_category_value cat
  from pride_files
  where project_accession in (select accession from study)
    and coalesce(primary_public_location_value,'') not like '%/generated/%')
select a accession,
  count(*) f_files, sum(sz) f_bytes, sum(dl) f_downloads,
  sum(dl) filter (where cat='RAW') f_raw_downloads,
  count(*) filter (where cat='RAW') f_raw_files,
  -- vendor raw formats
  bool_or(b like '%.raw') fmt_raw_ext,
  bool_or(b like '%.d' or n like '%.d.zip' or b like '%.tdf' or b like '%.baf') fmt_bruker_d,
  bool_or(b like '%.wiff' or b like '%.wiff2' or b like '%.wiff.scan') fmt_sciex_wiff,
  -- open spectrum formats
  bool_or(b like '%.mzml') fmt_mzml,
  bool_or(b like '%.mzxml') fmt_mzxml,
  bool_or(b like '%.mgf' or b like '%.ms2' or b like '%.mzdata' or b like '%.pkl' or b like '%.dta') fmt_peaklist,
  -- native tool outputs
  bool_or(b similar to '.*(proteingroups|evidence|msms|peptides|allpeptides|modificationspecificpeptides)\.txt' or n like '%mqpar%' or b like '%.apl' or (n like '%txt%' and n like '%combined%')) tool_maxquant,
  bool_or(b similar to '.*\.(msf|pdresult|pdstudy|pdresultview|msfview)') tool_pd,
  bool_or(b similar to '.*report(\.pr_matrix|\.pg_matrix|\.gg_matrix|\.stats|\.unique_genes_matrix)?\.(tsv|parquet)' or b like '%.speclib' or n like '%diann%' or n like '%dia-nn%' or n like '%dia_nn%') tool_diann,
  bool_or(b similar to '.*\.(sne|kit|bgsfasta)' or n like '%spectronaut%') tool_spectronaut,
  bool_or(b similar to '.*(combined_protein|combined_peptide|combined_ion|combined_modified_peptide|psm|protein|peptide|ion)\.tsv' or n like '%fragpipe%' or b like '%.workflow' or b like '%.pepxml') tool_fragpipe_tpp,
  bool_or(b similar to '.*\.(dat)' ) tool_mascot,
  bool_or(b similar to '.*\.(sky|skyd|blib)' or n like '%.sky.zip') tool_skyline,
  bool_or(b similar to '.*\.(sf3|sfd)') tool_scaffold,
  bool_or(b similar to '.*\.(group)') tool_proteinpilot,
  -- tabular result tables of any provenance
  bool_or(b similar to '.*\.(tsv|csv|xlsx|xls|txt|parquet)' and cat in ('SEARCH','RESULT','OTHER')) fmt_tables,
  bool_or(n like '%sdrf%') fmt_sdrf,
  bool_or(n like '%readme%') fmt_readme,
  bool_or(b like '%.fasta' or b like '%.fa') fmt_fasta
from f group by 1
""").df()

# --- reuse and attention signals -----------------------------------------------------
mentions = c.execute("""
select project_accession accession,
  count(distinct mentioning_pmid) filter (where not is_primary_reference) sec_mentions,
  count(distinct mentioning_pmid) filter (where not is_primary_reference and has_text_mined_accessions) sec_mentions_tm,
  min(mentioning_pub_year) filter (where not is_primary_reference) first_sec_year
from europepmc_mentions group by 1""").df()
icite = c.execute("""
select project_accession accession, max(citation_count) primary_citations, max(relative_citation_ratio) primary_rcr
from icite_citations where lookup_status='ok' or citation_count is not null group by 1""").df()
proj = c.execute("""
select project_accession accession, submission_date, publication_date release_date,
  date_diff('day', submission_date, publication_date) lag_days, total_file_downloads project_downloads,
  project_license from pride_projects""").df()

study = c.execute('select * from study').df()
vend = pd.read_csv(R / 'outputs/Global-proteomics-study/Vendor-patterns/Deposit-vendor-assignments.csv.gz',
                   usecols=['accession', 'vendor', 'credit'])
vend = vend[vend.accession.isin(study.accession)]
assert abs(vend.groupby('accession').credit.sum() - 1).max() < 1e-9, 'credit conservation'
vw = vend.pivot_table(index='accession', columns='vendor', values='credit', aggfunc='sum', fill_value=0)
vw.columns = ['credit_' + x.replace(' ', '_') for x in vw.columns]
single = vend.groupby('accession').vendor.agg(lambda s: s.iloc[0] if len(s) == 1 else 'Multiple').rename('vendor_single')

d = (study.merge(files, on='accession', how='left').merge(mentions, on='accession', how='left')
     .merge(icite, on='accession', how='left').merge(proj, on='accession', how='left')
     .merge(vw, left_on='accession', right_index=True, how='left')
     .merge(single, left_on='accession', right_index=True, how='left'))
for col in ['sec_mentions', 'sec_mentions_tm']:
    d[col] = d[col].fillna(0).astype(int)
assert len(d) == 36249 and d.accession.is_unique
# repeat-contributor sequence (hashed PI contact; ordered by release date then accession)
d = d.sort_values(['pi_group', 'release_date', 'accession'])
d['group_seq'] = d.groupby('pi_group').cumcount() + 1
d.loc[d.pi_group.isna(), 'group_seq'] = pd.NA
d['group_total'] = d.groupby('pi_group').accession.transform('size')
d['group_first_year'] = d.groupby('pi_group').release_year.transform('min')
d = d.sort_values('accession').reset_index(drop=True)
drop = ['title', 'description', 'keywords', 'processing', 'publication_title', 'experiments', 'software', 'quant']
d.drop(columns=drop).to_parquet(O / 'Deposit-features.parquet', index=False)
summary = {
    'deposits': len(d), 'with_file_rows': int(d.f_files.notna().sum()),
    'primary_outcome': 'standards == nongenerated_expanded (asserted)',
    'standards_equals_nongenerated_expanded': bool((d.standards == d.nongenerated_expanded).all()),
    'deposits_with_secondary_mention': int((d.sec_mentions > 0).sum()),
    'vendor_credit_conserved': True,
    'note': 'File-name evidence excludes /generated/ paths. Mentions are Europe PMC text-mined (largely open-access full text).'}
(O / 'Feature-summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
