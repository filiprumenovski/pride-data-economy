# Reproducibility coverage

## Frozen evidence

The public package begins with derived analytical inputs for 36,249 deposits. The seven early Parquets include reduced cohort, deposit-feature and dataset-layer tables, plus scale, breadth, model and species inputs. They preserve required numerical and categorical fields but exclude free-text protocols, original contact hashes and the contributor identity crosswalk. These public tables are deliberately not byte-identical copies of the richer private working tables.

## Executed checks

- Early analyses: 26 checks against frozen numerical outputs, including point estimates and the 1,000-replicate median raw-run ratio bootstrap. The ratio is 16/30, with a 95% percentile interval of 0.4848485–0.6296296.
- Late analyses: mention pair derivations and six outcomes, 168 rate rows, reused-pool composition, and the complete A1–A7 sensitivity result object. The adjusted design has 20,684 rows and 33 columns; the main fit and 200 bootstrap fits were finite and issued no warnings in the recorded run. Interval tails have limited Monte Carlo precision at 200 draws. The Thermo Fisher-only lower endpoint crossed one when the bootstrap sample was split in half; that near-null association should not be treated as robust evidence.
- Support analyses: six frozen summary tables covering archive/header recognition, software capability, documentation persistence and journal concentration, including the implemented bootstrap intervals.
- Figures: composition of all 19 figures from frozen blocks. This is independent of rebuilding every underlying analysis or reclassifying archive members.

New deterministic breadth resampling is saved separately by the early script. It does not replace the historical frozen intervals or rarefied-richness estimates in the manuscript. Contributor-group diversity is presented without the invalid historical confidence intervals.

## Fresh retrieval and offline reclassification

The external drive is now available. All 1,874 recovered source Parquets matched the source archive payload checksums and had readable footers. That broader PRIDE/GEO/MassIVE collection is retained locally; it is not redistributed wholesale in this release. The original October 2 extraction caches were not located.

On October 11–12, the frozen extraction targets were retrieved anew: 46,147 ZIP listings, 35,960,852 archive-member rows and 59,227 header records. All three files have checked SHA-256 digests in `data/extraction/Cache-manifest.json`; the compressed files are supplied separately as `Extraction-caches.zip`. No transient network failures remained at completion, but unsupported formats, inaccessible files and parsing failures remain explicitly recorded.

`release_caches.py` verifies file hashes, exact unique target coverage and member counts before applying the original classifiers offline. Compared with the frozen inputs, archive recognition outcomes are unchanged; one dataset has changed inner-content flags and six have changed unread-archive flags. Header evidence changes for five datasets, affecting final recognition for one (PXD060265). Every annual recognition percentage agrees at the reported one-decimal precision. The frozen manuscript results are retained; the fresh records are a dated retrieval sensitivity, not a byte-identical historical cache. Two archive targets and 145 header targets had observed file sizes differing from the frozen inventory.

## Additional clock validation

`release_clocks.py` recomputes eleven date differences for 36,249 datasets, checks the exclusion rules, and reproduces all nine main clock medians and ten software-specific medians, including their 1,000-draw cluster-bootstrap intervals. This starts from public extracted dates; it does not repeat date extraction from filenames.

## Remaining boundaries

- Geography inference from earlier cross-repository evidence is not reconstructed by the public offline scripts.
- Detection of instrument, software, acquisition-date and other upstream labels still begins from frozen derived inputs for most analyses.
- Historical breadth confidence intervals and rarefied richness are not reproduced exactly by the public deterministic resampling; contributor diversity has no inferential intervals.
- The full run-count sensitivity suite, timing sensitivity suite and every supplementary analysis are not yet covered by the public offline entry points.
- Reassembling 19 figures from frozen blocks verifies rendering, not all underlying statistics. The current 21-table workbook is preserved unchanged from the styled author copy.
- The Zenodo deposit/DOI is not yet confirmed. The BioRender graphical abstract needs a verified publication export and license.

The package therefore supports the specific executed checks above. It must not be described as a complete raw-record-to-paper reconstruction.
