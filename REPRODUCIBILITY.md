# Reproducibility coverage

## Offline entry points

`reproduce_all.py` runs the released numerical analyses, recomputes figure statistics, renders all 19 numbered figures at 7.2 inches and 600 dpi, and builds the 21 analytical and provenance table CSVs plus the curated supporting-comparisons table (S22). It requires the pinned Python environment, Liberation Sans, and the matching `Source-metadata.zip` and `Extraction-caches.zip` bundles. It writes to a new directory outside the repository. No PRIDE, publication or institutional API is queried during this run.

```sh
python -m pip install -r requirements.txt
python scripts/reproduce_all.py \
  --source /path/to/extracted-source \
  --cache-dir /path/to/extracted-caches \
  --output /path/to/new-reconstruction
```

On Linux, install `fonts-liberation`; elsewhere use `--font-dir /path/to/liberation-ttfs` if necessary. Font substitution is rejected for gallery rendering. `Stages.json` records each completed stage and its exit status. The final stage compares all 79 required release tables with the checked manuscript snapshot. Missing tables or numerical differences fail the run. Integer counts are exact; floating-point comparisons use relative tolerance 10⁻⁸ and absolute tolerance 10⁻⁹. `Numerical-comparison.json` records each result.

## Current validation

The local source rebuild completed all 22 analysis and rendering stages, then matched all 78 release tables. The 24-page manuscript PDF was rendered and visually inspected. The styled workbook was rebuilt from the current CSVs: all 7,634 data cells across 21 tables matched on reopening, and changed sheets were visually checked. `validation/` contains the stage logs, numerical comparisons, historical differences and artifact checks. A complete Linux run also completed all 22 stages. It exposed one platform-dependent density curve, traced to equal-sized FASTA candidates. The repaired database-date selection and dependent clock/figure outputs were rerun separately; the local eight-table comparison passed without changing the reported supplementary-table values. The release’s `Cloud-validation.json` records the independent runner checks and their exact commits.

## Source and inference checks

The source bundle contains 20 checksummed inputs. It includes 3,363,331 PRIDE file records, scientific project annotations, instrument and species metadata, publication records and links, and deidentified geography evidence. Cohort selection independently recovers the same 36,249 releases in 2013–2025. Source annotations reproduce all 22 software flags. Species flags and submission/release dates are checked separately.

Country reconstruction uses PRIDE contacts and historical MassIVE reference seeds. It reproduces all cohort country assignments, submitter-reported countries, inference flags and contributor-group partitions. Full affiliations, emails and original contact hashes are excluded. Retained equality identifiers and country suffixes support the inference without disclosing the original matching strings. The source exporter has been checked against all 20 input files by logical contents.

The run-count suite recomputes the raw-file classification, all run-count sensitivities and 71 early/later contrasts. The timing suite rebuilds publication and release dates, bulk-release checks, accession-order diagnostics, and timing sensitivities. DuckDB date casts explicitly use UTC. The acquisition-date parser, instrument harmonisation, storage analysis, software geography, documentation analyses, clocks and figure statistics run from the public source bundle and extraction observations.

## Historical and fresh evidence

The 7 April 2026 PRIDE metadata snapshot and 30 September 2026 literature snapshot remain fixed. The original October 2 archive-directory and top-level-header response caches were not recovered. A repeat retrieval on 11–12 October produced 46,147 ZIP records, 35,960,852 archive-member rows and 59,227 table-header records. Their SHA-256 values and exact target coverage are checked before reclassification.

Relative to the original extraction, nine additional archives have readable lists. Archive-level dataset recognition is unchanged. Header evidence changes for five datasets, with a final recognition change for one (PXD060265); annual recognition percentages remain equal at the reported one-decimal precision. The complete rebuild uses the repeat retrieval for cache-dependent results. Historical derived inputs and their dedicated validation entry points are retained so this change can be audited.

Some earlier bootstrap draws depended on SQL row order. Current breadth and timing resampling fixes sampling order. Point estimates remain comparable; intervals and rarefied richness can differ from the historical Monte Carlo draws. Equal-sized dated FASTA candidates use an explicit filename/location tie-break: the earlier size-only sort could select different dates on macOS and Linux. Contributor-group diversity is descriptive: invalid intervals from resampling contributor identities as their own categories are omitted from current outputs.

## Reuse analysis

The mention model now uses 2,000 contributor-group bootstrap draws, plus the full-cohort fit. Its design still has 20,684 rows and 33 columns. All 2,001 fits returned finite coefficients without fit warnings, with maximum absolute penalised score below 1e-5. The first 200 draws reproduce the preceding interval endpoints, and point estimates are unchanged. The Thermo Fisher-only interval now includes one in the full run and both split halves. The standard-file interval also includes one; DIA and human-sample intervals remain above one. The shared random stream advances through the additional model draws, so later A4/A5/A6 bootstrap intervals also changed while their point estimates and populations remained fixed. `validation/late/Bootstrap-upgrade.json` records these changes; `Bootstrap-odds-ratios.csv` preserves all 2,000 draws. Mentions are an incomplete attention proxy and do not establish reanalysis.

## Explicit boundaries

- Supplementary Table S21 is assembled from the archived findings and review ledger. Its 107 entries are historical assessments, not 107 independently rerun validation tests. The pipeline retains this provenance table; it does not repeat human or model-assisted adjudication.
- Software versions and ZIP-contained table-header flags begin at archived extraction observations; their original response bytes are not present. Downstream calculations can be rerun, but this is not a claim that every original network extraction is recovered.
- Scientific annotations, affiliations and filename/header rules are incomplete evidence, not ground truth. Restoring raw source records does not independently validate their scientific interpretation.
- Rebuilding analysis CSVs and figures does not by itself validate final PDF pagination or workbook styling. Those artifacts need their separate visual and cell-content checks.
- A Zenodo DOI is not yet confirmed. The BioRender overview remains a connector preview without a verified publication export/license. Author review and submission have not occurred.

## Editorial reporting table

Supplementary Table S22 was added on 17 October 2026 when secondary numerical comparisons were moved out of the Results. It preserves the previously reported values, populations, periods and uncertainty, with source references. Its checked-in source is `work/stack-and-reuse/reporting_comparisons.json`; the table builder includes it in offline reconstructions. This curated reporting table is not a new statistical analysis. The earlier 78-table and 21-sheet validation counts above describe the preceding package; the expected comparison list now contains 79 tables.

## Submission review update

The October 17 update changes the mention-model bootstrap and its dependent Figure 5, Table S20 and supporting comparisons in S22. Table S21 now distinguishes historical model-assisted checks from independent accuracy validation. Recognition shares are explicitly rule-based detections with unmeasured sensitivity and false-positive rates. The complete 22-stage local and Linux runs described above preceded this update; they are historical validation, not a claim that the entire source pipeline was repeated for these edits. The current model is checked separately by a complete repeat run, and the modified presentation artifacts are checked against its outputs. Author approval is still pending.
