# Source attribution and third-party rights

The study relies on public metadata and identifiers from PRIDE/ProteomeXchange, Europe PMC, OpenAlex and iCite, and on the scientific work of their contributors. Cite the original resource publications listed in the manuscript as well as this analysis. Public project accessions and publication identifiers are retained to make the analysis traceable.

- PRIDE: https://www.ebi.ac.uk/pride/ — archive and deposit metadata, file inventories and source links.
- MassIVE: https://massive.ucsd.edu/ — historical country-reference seeds, reduced to deidentified matching evidence.
- Europe PMC: https://europepmc.org/ — publication metadata and accession-mention evidence.
- OpenAlex: https://openalex.org/ — linked publication topics and institutional/geographic evidence.
- iCite: https://icite.od.nih.gov/ — publication citation metrics.

The MIT and CC BY licenses cover the author's original contributions, not a relicensing of third-party records or deposited experimental data. No bulk experimental spectra, full-text articles or raw contributor contact records are included. Upstream acquisition/normalization decisions are preserved in the analysis source and manuscript methods. The source-metadata bundle includes deidentified PRIDE contact evidence and MassIVE reference seeds needed to reconstruct geography. Equality identifiers and matched country suffixes are retained; full contact and affiliation strings are excluded. Public scientific project descriptions and protocols are retained as source metadata and may contain contributor-provided names or contact details.

Historical contributor-diversity bootstrap intervals in aggregate provenance tables are invalid and must not be interpreted. The current manuscript omits them and the public early comparison tables leave those cells empty; contributor-diversity point estimates remain descriptive.

## Bundled reference assets

- `work/stack-and-reuse/reference/pride-submission-validation.yaml`: unmodified PRIDE submission rules under Apache-2.0. The source commit, checksum and full license are retained in the adjacent reference README and `pride-resources-LICENSE`.
- `work/stack-and-reuse/ne_50m_admin_0_countries.geojson`: Natural Earth map data. Natural Earth places its map datasets in the public domain: https://www.naturalearthdata.com/about/terms-of-use/.

- `work/vendor-patterns/psi-ms.obo`: HUPO-PSI mass-spectrometry controlled vocabulary, version 4.2.2, under CC BY 4.0. Source: https://github.com/HUPO-PSI/psi-ms-CV. License retained as `reference/psi-ms-LICENSE`.
- `reference/iso3166-1.json`: historical pycountry ISO country-name database, unchanged from the geography compiler cache. Source: https://github.com/pycountry/pycountry. License retained as `reference/pycountry-LICENSE` (LGPL-2.1).

## Fresh extraction caches and graphical abstract

`Extraction-caches.zip` contains public archive-member names/sizes and table column headings, with source URLs and retrieval metadata. It contains no spectra or table data rows. These source-derived records retain the rights of their original contributors; they are not covered by a blanket grant of CC BY from this author. Filenames and column labels can contain contributor-provided identifiers. The wider PRIDE/GEO/MassIVE source collection on the author’s storage is not part of this release.

The BioRender graphical abstract (`Journal-figures/Graphical-abstract/PRIDE-overview.jpeg`, also embedded in the draft PDF) is excluded from this repository’s CC BY grant. It is a connector preview at 2752 pixels wide (about 382 dpi at 7.2 inches), not a verified publication export. Its attribution and provenance are retained; publication licensing must be confirmed before submission. The 19 numbered figure masters are 7.2 inches wide at 600 dpi.
