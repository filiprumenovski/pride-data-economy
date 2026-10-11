# Reproducibility coverage

## Frozen evidence

The public package begins with derived analytical inputs for 36,249 deposits. The seven early Parquets include reduced cohort, deposit-feature and dataset-layer tables, plus scale, breadth, model and species inputs. They preserve required numerical and categorical fields but exclude free-text protocols, original contact hashes and the contributor identity crosswalk. These public tables are deliberately not byte-identical copies of the richer private working tables.

## Executed checks

- Early analyses: 26 checks against frozen numerical outputs, including point estimates and the 1,000-replicate median raw-run ratio bootstrap. The ratio is 16/30, with a 95% percentile interval of 0.4848485–0.6296296.
- Late analyses: mention pair derivations and six outcomes, 168 rate rows, reused-pool composition, and the complete A1–A7 sensitivity result object. The adjusted design has 20,684 rows and 33 columns; the main fit and 200 bootstrap fits were finite and issued no warnings in the recorded run.
- Support analyses: six frozen summary tables covering archive/header recognition, software capability, documentation persistence and journal concentration, including the implemented bootstrap intervals.
- Figures: composition of all 19 figures from frozen blocks. This is independent of rebuilding every underlying analysis or reclassifying archive members.

New deterministic breadth resampling is saved separately by the early script. It does not replace the historical frozen intervals or rarefied-richness estimates in the manuscript. Contributor-group diversity is presented without the invalid historical confidence intervals.

## Deferred upstream cache transfer

The author confirmed on 11 October 2026 that these extraction caches are on a local NVMe drive that is currently disconnected, and deferred their transfer:

- `Archive-contents/Members.parquet`
- `Archive-contents/Listings.jsonl`
- `Header-sniff/Headers.jsonl`

They are not included in this release and have not been inspected or checksummed here. The public package retains derived recognition flags and aggregate evidence. A later release can add the frozen caches after transfer, integrity checking and scope/rights review. Their absence does not prevent the specific derived-input checks above; it does prevent claiming full reconstruction from raw extraction records.

Other boundaries: geography inference from earlier cross-repository evidence is not reconstructed; original extraction requires the source databases/caches; not every historical sensitivity is rerun by the offline entry points. Third-party data retrieval can change over time and is not equivalent to the frozen snapshot.
