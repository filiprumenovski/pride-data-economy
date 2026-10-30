"""Hand-check sheet for the file-layer classification (layers.py).

Draws 100 datasets at random (fixed seed), stratified 50/50 by release period (2013-2019,
2020-2025) and with at least 25 DIA datasets, and lists every non-generated file with the
layer it was assigned. A reviewer marks each dataset's layer flags as right or wrong; the
agreement rate is the classification accuracy we report.
"""
import numpy as np, pandas as pd
from common import R, O, SEED
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect
from layer_rules import RULES  # run layers.py first

OUT = O / 'Audit'; OUT.mkdir(exist_ok=True)
L = pd.read_parquet(O / 'Layers/Dataset-layers.parquet')
rng = np.random.default_rng(SEED)
dia = L[L.dia_expanded & (L.release_year >= 2020)].sample(25, random_state=SEED)
rest = L[~L.accession.isin(dia.accession)]
early = rest[rest.release_year <= 2019].sample(50, random_state=SEED + 1)
late = rest[(rest.release_year >= 2020)].sample(25, random_state=SEED + 2)
pick = pd.concat([dia.assign(stratum='DIA 2020-2025'), early.assign(stratum='2013-2019'), late.assign(stratum='2020-2025 other')])
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
c = connect(); acc = "','".join(pick.accession)
files = c.execute(f"""
with f as (select project_accession a, file_name, lower(file_name) n,
  lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
  coalesce(file_category_value, '') cat, file_size_bytes
  from pride_files where project_accession in ('{acc}') and coalesce(primary_public_location_value, '') not like '%/generated/%')
select a accession, file_name, cat pride_category, round(file_size_bytes / 1e6, 1) size_mb, {case} assigned_layer from f order by 1, 5, 2""").df()
files['reviewer_layer_if_wrong'] = ''
LAY = ['vendor_raw', 'open_spectra', 'standard_results', 'tool_results', 'unidentified_tables', 'sequence_db', 'spectral_library', 'sample_metadata']
ds = pick[['accession', 'stratum', 'release_year', 'dia_expanded', 'vendor_single'] + ['has_' + l for l in LAY]].copy()
ds.insert(1, 'pride_url', 'https://www.ebi.ac.uk/pride/archive/projects/' + ds.accession)
for l in LAY:
    ds['ok_' + l] = ''   # reviewer: Y if the has_ flag is right, N if wrong
ds['reviewer_notes'] = ''
with pd.ExcelWriter(OUT / 'Layer-audit-sheet.xlsx') as xw:
    pd.DataFrame({'instructions': [
        'One row per dataset in "datasets"; every non-generated file and its assigned layer in "files".',
        'For each layer, put Y in ok_<layer> if has_<layer> (True/False) is correct for that dataset, N if not.',
        'Layers: vendor_raw = instrument raw files; open_spectra = mzML/mzXML/MGF/MS2/DTA/PKL; standard_results = mzIdentML/mzTab;',
        'tool_results = known outputs of search/quant software (MaxQuant txt, PD .msf/.pdResult, Mascot .dat, DIA-NN report, Spectronaut .sne, FragPipe tsv, pepXML, Skyline...);',
        'unidentified_tables = other csv/tsv/xlsx/txt; sequence_db = FASTA; spectral_library = .speclib/.blib/.msp/.kit...; sample_metadata = SDRF.',
        'Optional: in "files", write the correct layer in reviewer_layer_if_wrong for misclassified files.']}).to_excel(xw, sheet_name='README', index=False)
    ds.to_excel(xw, sheet_name='datasets', index=False)
    files.to_excel(xw, sheet_name='files', index=False)
ds.to_csv(OUT / 'Layer-audit-datasets.csv', index=False); files.to_csv(OUT / 'Layer-audit-files.csv', index=False)
print(len(ds), 'datasets;', len(files), 'files;', ds.stratum.value_counts().to_dict())
