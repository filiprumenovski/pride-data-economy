# Inputs for the offline reconstruction

Download `Source-metadata.zip` and `Extraction-caches.zip` from the matching release. Extract them into separate directories. `Source-manifest.json` records the path, size, SHA-256 and extraction query or provenance for each source input.

The source bundle contains reduced PRIDE file inventories, instrument/species records, scientific project annotations, linked-publication metadata, deidentified geography evidence and a bootstrap-order map. Raw experimental spectra are not included. Public project prose may contain contributor-provided names or contact details; dedicated personal contact fields and original contact hashes are excluded.

Geography evidence includes 77,061 PRIDE contact records and 15,745 historical MassIVE reference-seed records, reduced to equality identifiers, country fields and matched country-name suffixes. The assigned country is not supplied as an input. The inference runner reconstructs all cohort countries, reported submitter countries, inference flags and contributor-group membership.

`Group-order.csv` preserves cluster membership, missingness and ordering relative to the original singleton fallbacks. It contains replacement labels, not the original identity crosswalk. Group labels remain linkable through public accessions and should not be described as anonymous.

Two observations remain at an extraction boundary: recovered software-version strings and flags from headers of tables inside ZIP files. Their released per-file/per-dataset observations support the downstream analyses; the original network-response caches for those two extractions were not recovered. The larger archive-directory and top-level-header caches were retrieved again on 11–12 October 2026 and are separately archived.

`export_source.py` documents and repeats the reduction from the author's original database snapshot. Re-exported Parquet/GZIP bytes may vary by serialization metadata; the recorded export check compares logical table contents. Readers should use the pinned release bundles for exact input hashes.
