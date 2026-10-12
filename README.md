# The data economy of public proteomics

Analysis code and frozen derived inputs for **The data economy of public proteomics: thirteen years of PRIDE**.

Filip Rumenovski · Department of Chemistry, Wayne State University · filip.rumenovski@wayne.edu

This study describes 36,249 public PRIDE deposits released during 2013–2025. The PRIDE snapshot is dated 7 April 2026; literature enrichment is dated 30 September 2026. Deposits, publications and recurring contributor groups are distinct units. The manuscript is a pre-submission research draft, not a peer-reviewed publication.

[Read the manuscript PDF](outputs/Stack-and-reuse/Paper-2col.pdf) · [Manuscript source](outputs/Stack-and-reuse/Paper.md) · [Figure legends](outputs/Stack-and-reuse/Journal-figures/Legends.md)

## Contents

- `data/early/`: seven Parquet tables containing public analytical features, operational file layers, experimental scale, breadth, model credits and species credits; frozen comparison tables and a field inventory.
- `data/late/`: dataset-level flags, publication links, mention pairs, model covariates, and comparison tables. These are CSV inputs; they need no private database.
- `data/extraction/`: frozen target lists, cache checksums and the expected fresh-retrieval comparison. The large cache files are in the separate `Extraction-caches.zip` release asset.
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
python scripts/release_clocks.py --data data/late --output results/clocks
python work/stack-and-reuse/compose_journal.py --no-render
```

Use fresh output directories for each run. The first analysis recomputes 26 checks including annual breadth points, model credits, software-country contrasts and the 1,000-replicate headline run-count bootstrap. The late analysis recomputes mention eligibility, rates, and all 200 adjusted-model bootstrap replicates. The support analysis checks six downstream tables, including software export capability and journal concentration intervals. The clock analysis recomputes eleven date differences and nineteen median summaries with 1,000 bootstrap draws each. Figure composition uses frozen blocks at 7.2 inches wide and 600 dpi.

These checks establish specific numerical reproducibility from frozen derived inputs. They do not reconstruct every classification decision or every manuscript sensitivity analysis from raw files. See `REPRODUCIBILITY.md` for the coverage and limits.

## Data interpretation

Contributor IDs are replacement labels for recurring contact-derived groups, not named laboratories or verified independent research teams. Their private identity crosswalk is excluded. These labels are pseudonymous and remain linkable through public project accessions; they are not guaranteed anonymous. No detection of a method or file type means “not detected”; it does not establish absence. File recognition does not establish validity, usability or actual reuse. Standard-result flags exclude repository-generated paths. Later accession mentions are an incomplete attention proxy and do not prove reanalysis.

## Licensing and citation

Code: MIT (`LICENSE`). Author-created derived tables, manuscript and figures: CC BY 4.0 (`DATA-LICENSE`). Third-party source records retain their original terms; see `NOTICE.md`. Citation metadata is in `CITATION.cff`. A Zenodo DOI will be added once the archive exists.

## Rebuild the manuscript PDF

After figure composition, install Chrome or Chromium and point `CHROME` at its executable:

```sh
CHROME=/path/to/chrome python work/stack-and-reuse/build_paper_2col.py
```

The checked PDF uses the existing 600-dpi figures without downsampling. Font availability can change pagination; inspect a fresh render before submission.

## Reclassify the recovered extraction caches

Extract `Extraction-caches.zip` from the matching release, then run:

```sh
python scripts/release_caches.py --cache-dir /path/to/extracted-caches --output /path/to/new-results
```

This offline run verifies checksums and complete URL coverage, reclassifies 35,960,852 archive members and 59,227 header records, and checks the differences from the frozen manuscript inputs. Allow several GB of disk space for intermediate files. The fresh retrieval changes final recognition for one dataset; all annual recognition percentages remain equal at the reported one-decimal precision. See `validation/Extraction-reclassification.json`. The release remains a draft until the archive handoff and submission checks are complete.
