# The data economy of public proteomics

Analysis code and frozen derived inputs for **The data economy of public proteomics: thirteen years of PRIDE**.

Filip Rumenovski · Department of Chemistry, Wayne State University · filip.rumenovski@wayne.edu

This study describes 36,249 public PRIDE deposits released during 2013–2025. The PRIDE snapshot is dated 7 April 2026; literature enrichment is dated 30 September 2026. Deposits, publications and recurring contributor groups are distinct units. The manuscript is a pre-submission research draft, not a peer-reviewed publication.

## Contents

- `data/early/`: seven Parquet tables containing public analytical features, operational file layers, experimental scale, breadth, model credits and species credits; frozen comparison tables and a field inventory.
- `data/late/`: dataset-level flags, publication links, mention pairs, model covariates, and comparison tables. These are CSV inputs; they need no private database.
- `scripts/`: offline reproduction entry points. See the commands below.
- `work/stack-and-reuse/`: original analysis, extraction and figure source for provenance. Some historical scripts require upstream inputs outside this release; they are not the offline entry points.
- `outputs/Stack-and-reuse/`: aggregate evidence, frozen figure blocks, final figure previews and manuscript materials.
- `validation/`: archived numerical validation reports.

## Reproduce the checked analyses

Use Python 3.12 or newer. Exact validation used Python 3.14 and the pinned versions in `requirements.txt`. Run from the repository root:

```sh
python -m pip install -r requirements.txt
python scripts/verify_release.py
python scripts/release_early.py validate --data-dir data/early --output-dir results/early
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/release_late.py validate --data data/late --output results/late
python scripts/release_late_support.py --data data/late --output results/late-support
python work/stack-and-reuse/compose_journal.py --no-render
```

Use fresh output directories for each run. The first analysis recomputes 26 checks including annual breadth points, model credits, software-country contrasts and the 1,000-replicate headline run-count bootstrap. The late analysis recomputes mention eligibility, rates, and all 200 adjusted-model bootstrap replicates. The support analysis checks six downstream tables, including software export capability and journal concentration intervals. Figure composition uses frozen blocks at 7.2 inches wide and 600 dpi.

These checks establish specific numerical reproducibility from frozen derived inputs. They do not reconstruct every classification decision or every manuscript sensitivity analysis from raw files. See `REPRODUCIBILITY.md` for the coverage and limits.

## Data interpretation

Contributor IDs are replacement labels for recurring contact-derived groups, not named laboratories or verified independent research teams. Their private identity crosswalk is excluded. No detection of a method or file type means “not detected”; it does not establish absence. File recognition does not establish validity, usability or actual reuse. Standard-result flags exclude repository-generated paths. Later accession mentions are an incomplete attention proxy and do not prove reanalysis.

## Licensing and citation

Code: MIT (`LICENSE`). Author-created derived tables, manuscript and figures: CC BY 4.0 (`DATA-LICENSE`). Third-party source records retain their original terms; see `NOTICE.md`. Citation metadata is in `CITATION.cff`. A Zenodo DOI will be added once the archive exists.
