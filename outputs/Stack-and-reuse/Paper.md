# The data economy of public proteomics: thirteen years of PRIDE

## Abstract

PRIDE's public archive makes changes in proteomics research practice visible across an expanding community. We examined all 36,249 datasets released between 2013 and 2025 (2.71 petabytes), reading file lists, 45,985 results archives and table headers. Paper links covered 75% of datasets and countries could be resolved for 78%. Effective numbers of contributor groups and publication outlets grew 5.0-fold and 5.4-fold, respectively, between the periods 2013 to 2015 and 2023 to 2025. Human and mouse accounted for increasing shares of deposits as species diversity fell by a third. Large deposits became more numerous, while the median raw-run count declined, reaching 16 in 2023 to 2025 against 30 in 2013 to 2016. Annual medians for dataset size stayed within 19 to 25 GB. We recognised identification results in 80% of 2025 deposits by reading tool-specific outputs, while 16% included standard result files. In the subset with dated raw files, the median acquisition-to-release interval was 2.2 years; most elapsed before submission. Qualifying non-primary papers more often named human datasets and those acquired by data-independent acquisition (DIA). We did not detect an association of these mentions with standard result files; mentions alone cannot confirm reanalysis. PRIDE's accumulated deposits provide a foundation for improving discovery of the analytical content already shared by its contributors.

## Introduction

Researchers package proteomics outputs for public release and publication. Later work can build on them. We use data economy to mean this circulation of data after an experiment. Costs and prices are outside our scope. PRIDE, the largest repository in the ProteomeXchange consortium, preserves mass-spectrometry proteomics data from a growing community, makes them publicly accessible and links deposits to the published literature [PRIDE2025] [PX2026].

Evidence in PRIDE's deposited files connects the growth recorded in repository reports to changing experiments and analytical workflows. Raw files record experimental runs and instrument labels. We identify the analysis software from logs and parameter files. With dates and publication links, these records trace the history of a study and make changes in research practice visible across the archive.

An archive-wide view of PRIDE complements studies of public-data impact through citations [PerezRiverol2019] and examinations of selected deposits, including six reanalysed projects [Vadadokhau2026] and an audit of 500 datasets [Leite2026]. Deposited files also show how developments described elsewhere, including large cohorts [Geyer2017] [Messner2023] and common result formats [PSI2023], appear in the public record.

We examine all 36,249 PRIDE datasets released between 2013 and 2025. We follow the expanding contributor community through deposit scale and recorded workflows, then read result files to identify analytical content beyond standard formats. Literature links extend this account to publication and later attention. The analysis describes deposits rather than instrument ownership or adoption dates. We report selection denominators and, where computed, 95% intervals from resampling contributor groups. Later accession mentions provide an incomplete proxy for reuse and do not establish reanalysis.

## Results

### PRIDE reached a broader contributor community

PRIDE expanded across contributor groups and the journals publishing their work. Comparing 2023 to 2025 with 2013 to 2015, the effective number of groups increased 5.0-fold and that of journals 5.4-fold (95% interval 4.7 to 6.3) (Fig. 1a,b; Supplementary Tables S2 and S3). By publication year, the number of distinct journals carrying linked papers rose from 78 in 2013 to 627 in 2025 (Supplementary Table S18). Effective numbers account for both the number of categories and their relative shares (Hill numbers; Box 1).

Broader participation accompanied more even representation of instrument models. Their effective number increased 1.8-fold (1.6 to 2.0), although rarefied richness did not increase (61 against 66 distinct models per 1,500 datasets). Species diversity moved in the opposite direction, reaching 0.67 times its earlier value (0.58 to 0.81). For each of these four dimensions, the direction of change agreed between orders 1 and 2 (Supplementary Table S3).

At order 2, Supplementary Table S3 shows no detectable change in the effective numbers of analysis tools, research subfields, countries and instrument vendors. The annual series revealed different timing: instrument models increased mainly from 2020, journals until 2022, and tools only in 2024 and 2025 (Fig. 1a).

Human and mouse accounted for most of the concentration in species. Human samples were named in 41% of datasets released in 2013 and 52% of those released in 2025; mouse was named in 11% and 19%. Datasets naming neither species made up 46% of deposits during 2013 to 2016, compared with 33% during 2023 to 2025 (Supplementary Fig. 9). Among the 354 groups contributing in both periods, this share also declined (from 40% to 32%), so contributor turnover explained only part of the change.

New contributors continued to supply a substantial share of PRIDE's deposits. Groups with no earlier PRIDE release since 2013 contributed between 42% and 48% of annual deposits during 2020 to 2025 (44% in 2025; Fig. 1c,d). Two-thirds of all groups deposited a single dataset, accounting for a quarter of all datasets (Table S5; Supplementary Fig. 1).

Datasets from groups with a resolved country increased 3.3-fold when comparing 2023 to 2025 against 2016 to 2018. Half of the additional datasets came from Europe (Fig. 1e to g; Supplementary Table S6). Datasets from China increased 5.2-fold over these periods, compared with 2.8-fold for Germany. China's share of country-resolved datasets rose from 3% in 2013 to 15% in 2025. These counts concern PRIDE only and exclude groups that deposit only in iProX, jPOST or MassIVE.

### Archive growth encompassed experiments at different scales

PRIDE's growth [PRIDE2025] accommodated experiments at different scales. The median number of raw runs per dataset was 16 for 2023 to 2025 (n = 15,248), compared with 30 for 2013 to 2016 (n = 2,498), a ratio of 0.53 (95% interval 0.48 to 0.63). The whole distribution shifted towards fewer runs (Fig. 2a).

Large deposits nevertheless became more numerous. For datasets with at least 100 runs, the annual mean count was 124 during 2013 to 2016 and 459 during 2023 to 2025, while their share fell from 19.8% to 9.0% (Fig. 2b). Such deposits are relevant to the scale proposed for cohort proteomics [Geyer2017] [Messner2023]. Run counts differ from sample counts because isobaric labelling places several samples in one run: the mean number of channels per labelled dataset increased from 6 in 2014 to 13 in 2025 as TMTpro replaced iTRAQ (Supplementary Fig. 8).

The decrease in median runs persisted across five sensitivity analyses: combining datasets from the same paper; counting exact runs only; collapsing fraction-numbered runs; reweighting to the 2013 to 2016 composition of vendor, acquisition mode, labelling, species and submission type; and stratifying by country (Fig. 2c; Supplementary Table S8). TMT-labelled, SCIEX and Bruker datasets showed no detectable decrease in median raw runs (early Bruker n = 56).

Within the 3,456 groups contributing in multiple years, the within-group slope of log2 raw runs was about 38% of the pooled slope on the same 23,486 datasets. For 2018 to 2025 releases, the matched geometric-mean ratio of run counts in first-time versus returning groups was 0.64 (0.61 to 0.68), stratifying by release year, vendor, acquisition mode, human-sample status and continent.

Individual runs occupied more storage. Median raw-run size increased from 0.35 GB in 2013 to 1.07 GB in 2025. Median raw-run sizes were 0.44 GB on the LTQ Orbitrap and 3.95 GB on the Orbitrap Astral (Fig. 2d; Supplementary Table S7). Among single-model datasets with exact run counts, a decomposition of mean log10 bytes per run compared 2014 to 2016 with 2023 to 2025. The changing model mix contributed 2.30-fold and the within-model component 0.85-fold (0.80 to 1.10). With fewer but larger runs, median dataset size remained between 19 and 25 GB in every year, while the mean and upper tail grew.

Storage was concentrated in a small fraction of deposits. The largest 1% held 24% of the 2.71 PB and the largest 10% held 61% (Fig. 2e). Vendor raw files accounted for 84% of all bytes (Fig. 2f; Supplementary Fig. 6).

### Deposits recorded technological change alongside established workflows

PRIDE’s instrument labels document changing workflows alongside persistent earlier technologies. Across years, Thermo Fisher labels held 79 to 89% of model-labelled dataset credit; the repository reports a similar share of stored raw files [PRIDE2025]. The LTQ Orbitrap accounted for 40% of fractional model-labelled dataset credit in 2013, the first year of the series. Successive leading models had progressively lower peak shares of this credit, with 38% for the Q Exactive and 19% for the Q Exactive HF (Fig. 3a; Supplementary Fig. 2).

On harmonised model names, the effective number of models rose from 4.7 (95% interval 4.0 to 5.7) in 2013 to 12.1 (11.5 to 12.8) in 2025 (Fig. 3b). Every model named in at least 50 datasets was still represented among 2025 deposits. This record includes reanalyses and data released after a delay from the experiment. Datasets labelled with timsTOF models represented 13% of the 2025 total, compared with 1.8% in 2021. The family had been described in 2021 as becoming widely used [Meier2021].

New models appeared in deposits from Asia about one year later than expected from its dataset count. Africa's longer delay was close to that expected from its small number of datasets (Supplementary Fig. 3). The null model used for these comparisons ignores clustering within contributor groups.

DIA-NN grew rapidly among deposited analysis tools. It was named in 2 datasets released in 2020 and 722 in 2025, ranking third behind MaxQuant (28% of datasets) and Proteome Discoverer (20%) (Fig. 3c).

Tool use also varied by country. In 18 of 20 pairs of a tool and its developers’ country, groups in the developers' country named the tool more often than groups in all countries combined. The median ratio across pairs was 1.8-fold (permutation P = 0.0002; Fig. 3d). For 11 pairs the 95% interval excluded a ratio of 1. Examples include Spectronaut in Switzerland (3.4-fold), Comet in the United States (2.9-fold) and MaxQuant in Germany (1.5-fold). The developers' own laboratories were not excluded. A similar pattern had been remarked for open-source pipelines in western Europe [Orsburn2021].

Analysis files preserved the software versions used. Parameter files and logs showed that the MaxQuant 1.6 release line, first recorded in PRIDE in 2017, had a share of 47% among MaxQuant datasets released in 2025 for which a version could be recovered. Among 2025 submissions with a recoverable version, the share was 40%; versions were recoverable for 361 of 576 datasets naming MaxQuant. DIA-NN 1.8 accounted for 61% (48 to 74%) of DIA-NN datasets submitted in 2025 with a recoverable version. Versions were recoverable for 110 of 331 submissions naming DIA-NN (Supplementary Fig. 4 and Table S12).

Measurement age tracked the release lines in use. We compared MaxQuant datasets submitted in 2023 to 2025 by the age of their newest raw file at submission. A 1.x release was used in 35% of datasets when that file was less than one year old and in 80% when it was more than three years old.

### Tool-specific files revealed analytical content beyond standard formats

The deposited files contained far more identifiable results than standard formats alone revealed. Among datasets released in 2025 (n = 6,231), reading table headers brought detection to 80%, compared with 16% containing standard result files. Detection had reached 59% from file names and 68% after inspecting zip contents (Fig. 4c; Supplementary Table S13). Even at full depth, the detectable share fell from 92% in 2015 to 80% in 2025 (Fig. 4d). These searches provide lower bounds on content. Detection alone cannot show whether a file is valid or usable.

Proteomics Standards Initiative formats provide common representations of results; mzIdentML has been implemented in over 35 software tools [PSI2023]. We counted submitter-deposited mzIdentML or mzTab as standard result files, excluding repository-generated files (Box 1); these appeared in 16 to 20% of datasets each year since 2015 (Supplementary Table S13).

Standard-file rates differed across the named analysis tools. Among single-tool deposits, rates ranged from 0.2% for DIA-NN to 76% for MS-GF+. Spectronaut and MaxQuant also had low rates; Mascot deposits more often included standard files (Fig. 4a). Supplementary Table S9 gives the tool-specific rates.

Datasets from tools writing standard files by default had a deposition rate of 55% (95% interval 44 to 66%). The rate was 15% (14 to 18%) for datasets from tools offering them as an option and 3% (2 to 4%) for datasets whose analysis tools could not write these files (Fig. 4b; Supplementary Table S10). Only MS-GF+ and MetaMorpheus comprised the default group (144 datasets), and MetaMorpheus deposits included standard results in 23% of cases despite the software writing them by default.

Software composition also changed over time. Between the 2019 to 2021 and 2023 to 2025 periods, the archive-wide standard-file rate changed by −1.4 percentage points (−3.0 to 0.3). A decomposition by software family assigned −3.0 points (−3.4 to −2.5) to changing composition, mainly the growth of DIA-NN and Spectronaut and the decline of Mascot. Change within families contributed +1.6 points (0.1 to 3.2).

Among datasets with no DIA detected, the standard-file rate remained at 19 to 21% each year since 2020. DIA grew from 2.8% of datasets in 2015 to 21.1% in 2025, but about 3% of DIA datasets released in 2023 to 2025 included standard result files. Walzer et al. described heterogeneous spectral-library formats as one obstacle to DIA reanalysis [Walzer2022].

Scan-level files and intermediate search outputs occupied 66% of uncompressed results-archive bytes, compared with 3.9% for final result tables (Supplementary Fig. 6). Open spectrum formats appeared in 34% of datasets in 2013 to 2016 and 21% in 2023 to 2025. The halving of complete submissions, which must include peak lists, accounted for 5 to 10 of the 13 percentage points (Supplementary Fig. 7).

Documentation was uncommon but often recurred within contributor groups. When a group's previous dataset included a README, its next did so in 41% of cases, against an expected 1.5% for datasets released in the same years. For sample and data relationship format (SDRF) files, the corresponding rates were 44% and 2.5% (Fig. 4e; Supplementary Table S14). These recurrence rates are unadjusted for pairs submitted on the same day or for the same paper; some SDRF files were added to several datasets of a group after release. Across the archive, SDRF presence rose from 0.7% in 2020 to 2022 releases to 2.7% in 2023 to 2025, consistent with reports of limited use [VanDenBossche2026].

### Release dates placed public deposits in the history of an experiment

Dates preserved in deposited files connected experiments and analyses to public release (Fig. 5a; Supplementary Table S15). Raw-file names supplied an acquisition date for 11,868 of the 35,364 datasets with raw files (34%). In this subset, the median time from the latest acquisition to public release was 2.2 years (95% interval 2.1 to 2.3).

Deposited files also dated analysis resources. The date in the main sequence database file name (n = 2,558) preceded public release by a median of 2.5 years (2.4 to 2.6). The median age of the software release used (n = 10,466), measured from its first appearance in PRIDE, was at least 3.4 years (3.3 to 3.5).

Of the 1,006 datasets with both dates, 71% had a database file dated earlier than the latest acquisition. The software release had appeared in PRIDE before the latest acquisition in 79% of 3,986 datasets.

Most of the acquisition-to-release interval elapsed before submission (Supplementary Table S15). Submission preceded the linked paper by 141 days at the median; the median delay between the paper and public release was 25 days.

Public releases also clustered in time. From October 2018 onwards, at least 50 datasets were made public on each of 62 days. Datasets released on these days became public a median of 95 days after their paper, compared with 11 days on other days (Supplementary Fig. 13 and Table S16). Who set each release date is unrecorded.

### PRIDE-linked publications reached beyond specialist journals

PRIDE-linked papers reached an increasingly broad range of journals. Titles containing "proteom" published 51% of these papers in 2013 and 10% in 2025, even as their count of linked papers rose from 142 to 259 (Fig. 5b,c; Supplementary Tables S17 and S18). The ten leading journals' combined share fell from 66% to 31% over the same period. In 2025, *Nature Communications* was the largest single journal (Supplementary Table S17).

The research fields of linked papers also shifted. Chemistry (almost entirely spectroscopy, the subfield to which mass-spectrometry methods papers are assigned) fell from 26% of linked datasets in 2013 to 6% in 2025. Medicine rose from 18% to 29% (Supplementary Fig. 12). Coverage limits these comparisons: paper links were available for 94% of datasets released in 2013 but 61 to 67% of those released since 2021.

### Accession mentions provided evidence of continued attention to public data

Europe PMC text mining identified papers that named PRIDE accessions. A qualifying paper was not marked as the primary reference, named datasets from at least two contributor groups and mentioned no more than 100 datasets. Such papers named 2.3% of the 20,684 datasets released in 2013 to 2022 during the release calendar year or the next three calendar years (Fig. 5d; Supplementary Table S19). Publication and release order within a year was unresolved. This proxy for reuse is incomplete, especially for accessions confined to supplementary tables or literature outside open-access full text. A mention does not establish reanalysis.

Qualifying papers named DIA and human datasets more frequently. Adjusted odds were higher for DIA datasets (odds ratio 1.71; 95% interval 1.22 to 2.35) and for human datasets (1.39; 1.17 to 1.71). The Thermo Fisher-only point estimate was also above one, with an interval close to the null (Fig. 5e; Supplementary Table S20). The model accounts for the other dataset characteristics and study-level covariates specified in Methods; its DIA cohort includes releases only up to 2022.

We did not detect a relationship between standard result files and qualifying mentions (adjusted odds ratio 0.89; 0.66 to 1.17). The crude rate ratio was 1.04 (0.72 to 1.39) when the window covered the release calendar year and the next five calendar years, using datasets released by 2020. A restricted cohort of human datasets assigned solely to Thermo Fisher, with no DIA detected and naming MaxQuant or Proteome Discoverer, gave an imprecise crude rate ratio of 0.74 (0.26 to 1.38) (Supplementary Table S20). A modest advantage therefore cannot be excluded.

SDRF files require a different interpretation. Datasets with SDRF files had higher odds of a qualifying mention (4.08; 2.29 to 7.26), but 91 of the 205 were released before the format existed. Those files were added later, often by reanalysis projects, so the direction of the association cannot be determined. Methods and machine-learning papers were among those naming the largest numbers of datasets, alongside resource papers.

## Discussion

Thirteen years of PRIDE deposits document the work of an expanding proteomics community. New contributors continued to account for a substantial share of annual releases, and linked papers appeared across a much wider range of journals. This collective record, built by the repository and its contributors, makes changes in research practice visible at a scale difficult to achieve through individual studies.

File-level analysis gives PRIDE's growth an experimental context. PRIDE accommodates increasing numbers of large deposits alongside many datasets with fewer raw runs. The mix of instrument models helps explain why a typical run occupies more storage. Older models remained represented in deposits, which include delayed releases and reanalyses. Deposited software versions also preserve evidence of earlier workflows. Participation can broaden even as biological coverage becomes more concentrated on human and mouse.

Identification results could be recognised in 80% of 2025 deposits when tool-specific outputs were examined, although 16% contained standard result files. Recognising these outputs more systematically could make contributors' results easier to discover. File presence does not establish validity or readiness for reuse, so extraction would need validation and clear provenance.

Software-specific deposition rates suggest where improvements could be tested collaboratively with tool developers. Standard files accompanied 55% of datasets whose tools write them by default; for tools lacking this capability, the rate was 3%. The MetaMorpheus result shows that a default alone is insufficient, and tool choice may reflect depositor intentions. Sample annotations supplied during analysis could populate SDRF files and reduce reconstruction work for later users [Dai2021] [VanDenBossche2026]. Capturing instrument, software version and database date from deposited files would help readers assess provenance, a concern addressed in reporting recommendations for single-cell proteomics [Gatto2023].

Release timing helps explain which stages of a study PRIDE records. Most of the measured acquisition-to-release interval elapsed before submission, so it cannot be attributed wholly to repository handling. Older software often accompanied older measurements. A public release date therefore needs to be read alongside the available acquisition and analysis dates.

Qualifying non-primary accession mentions provide an incomplete view of attention to deposits. These papers more often named DIA and human datasets. Standard result files showed no detectable association with those mentions, but the intervals do not exclude a modest benefit. Mentions alone do not establish reanalysis. Calendar-year mention windows cannot establish within-year chronology. The SDRF association also requires caution because many files were added after release, sometimes during reanalysis. Evaluations of future improvements should record when changes occur and distinguish subsequent use from metadata added during that use.

The scope and coverage of the evidence constrain these conclusions. PRIDE represents only part of ProteomeXchange, which matters especially for geographical comparisons. Hashed principal-investigator contacts approximate contributor groups: changed addresses can create apparent new groups, while core facilities serve several laboratories. First-time status is left-censored at 2013. File and text searches provide lower bounds on the presence of software, acquisition modes, labelling and identifiable results. Paper links cover 75% of datasets overall but 61 to 67% of recent releases; dated raw-file names cover 34% of deposits containing raw files. Text mining also misses later mentions, particularly those confined to supplementary material. The reconstructed snapshot cannot recover the history of metadata edits or withdrawn datasets and excludes recent submissions that remain private. Our comparisons are descriptive; the Supplementary Tables give denominators and intervals where computed, and Table S21 records the sensitivity analyses.

Making the analytical results already shared by PRIDE's contributors easier to find and documenting their provenance would extend the value of that work.

## Box 1. Data, definitions and limits

**Cohort.** The 7 April 2026 PRIDE snapshot contains 36,249 public datasets with release dates between 2013 and 2025 (Supplementary Table S1). Release year is the year of public release in PRIDE; where submission year is used instead, this is stated.

**Units.** Results are reported per dataset unless stated otherwise. Papers denote the OpenAlex works matched to dataset publication records. Contributor groups are defined by hashed contacts of principal investigators; a first-time group is a group with no earlier PRIDE dataset released since 2013. When a dataset names several instrument models, tools or species, each receives an equal fractional share.

**Effective number.** Category shares pᵢ among datasets define the Hill numbers. We use 1/Σpᵢ² at order 2. At order 1 the expression is exp(−Σpᵢ ln pᵢ). Order 0 measures richness, rarefied here to equal numbers of datasets.

**File evidence.** Top-level file lists (3.03 million files); members of 45,985 results archives, read from zip central directories by HTTP range requests; column headers of deposited tables, matched to output signatures of analysis tools; software versions from parameter files, logs and workflow files; dates from raw and database file names. Each is a lower bound on file content.

**Standard result files.** mzIdentML or mzTab files deposited by submitters, including those identified by expanded name matching; files generated by the repository are excluded.

**Uncertainty.** Bootstrap samples of contributor groups supplied 95% percentile intervals (seed 20261001; 200 to 2,000 replicates depending on the analysis). Contributor-group effective numbers are descriptive, without bootstrap intervals. Mention windows use calendar years and do not resolve within-year chronology.

**What the data cannot show.** Instrument ownership or purchase dates; laboratory size; the reasons for packaging decisions; who set a release date; reuse that leaves no text-mined mention. Country is unresolved for 22% of datasets.

## Methods

**Data sources.** The analysis used a frozen snapshot (`data/Data-manifest.json`) of the PRIDE archive tables (projects, files, instruments, species and publications), ProteomeXchange registry links, OpenAlex records of linked papers (works, authorships and institutions) and Europe PMC records (journal titles and accession mentions).

**Reading of files.** Members of zip archives were listed from their central directories, including zip64 archives, using byte-range requests. Column headers were read from the first bytes of deposited text and xlsx tables and matched to signatures of the outputs of analysis tools. Software release strings were parsed from parameter files and logs. Archive members were classified with the same name-based rules as top-level files. Header signatures underwent model-assisted review and revision; this did not provide an independent estimate of classification accuracy. Recognition of a file type does not verify its contents or usability.

**Runs.** Raw runs are vendor raw acquisitions identified from file names and extensions. Bundled formats were counted as single runs, so counts are lower bounds; exact counts are reported separately.

**Instrument models.** Instrument labels were harmonised to models; labels that name a series, a vendor or no instrument were excluded from model-level analyses.

**Statistical uncertainty.** We obtained 95% percentile intervals by bootstrapping contributor groups. For non-group category axes, ratios of effective numbers between periods used groups resampled independently within each period. Effective numbers of contributor groups and their period ratios are descriptive point estimates without bootstrap intervals: resampling groups while retaining their category identities biases that concentration statistic.

**Comparisons.** Shift-share decompositions separate changes within strata from changes in composition. Rates of qualifying accession mentions were directly standardised to the release-year × human-status composition of the cohort. The expected first appearance of instrument models by continent was obtained by drawing model labels independently of continent within each release year.

**Accession-mention outcome and model.** We modelled a qualifying non-primary accession mention in the release calendar year or the next three calendar years for datasets released in 2013–2022. Qualifying papers named deposits from at least two contributor groups and no more than 100 datasets. Same-year chronology could not be resolved. A ridge-penalised binomial logistic model (penalty 10⁻⁴) included standard result files, detected DIA, SDRF, human samples and a Thermo Fisher-only assignment from instrument labels. Adjustment included categorical release year, grouped software family and research field, log(1 + citations), a missing-citation indicator and log₁₀(dataset bytes). The model included 20,684 datasets, 483 outcomes and 8,363 contributor groups. Percentile intervals used 200 resamples of whole contributor groups. Supplementary Methods specify transformations, missingness and reference levels. Sensitivity analyses reported unadjusted rate ratios in a restricted cohort with no DIA detected and over the release calendar year plus five subsequent calendar years.

**Sensitivity analyses.** Each principal result was re-examined using alternative definitions, units of analysis and subsets. The outcome of each check (supported, weakened or refuted) and the resulting wording are recorded in Supplementary Table S21.

**Figures.** Figures were generated from code using frozen analytical outputs; captions state the relevant selection criteria, denominators and estimands. Raster exports have a width of 7.20 inches at 600 dpi.


**AI assistance.** OpenAI Codex assisted with code development and review, figure preparation, literature checking and manuscript drafting and editing. The author’s final verification and approval of the analysis and text remain pending before submission.

## Supplementary methods

**Mention ascertainment.** Accession–PMID pairs marked as primary references were excluded. The main outcome required a mentioning paper to name no more than 100 datasets spanning at least two contributor groups. These paper-level counts were computed before filtering on publication-year minus release-year lags of 0–3 inclusive. The main outcome does not additionally exclude primary papers of another deposit from the same group. Same-year mentions can precede the exact release date. A five-subsequent-calendar-years sensitivity used lags 0–5 among 2013–2020 releases (11,741 datasets). The restricted cohort comprised 4,399 human datasets assigned solely to Thermo Fisher from instrument labels, with no DIA detected, naming MaxQuant or Proteome Discoverer. Neither restriction establishes raw-file presence or confirmed DDA.

**Model construction.** The binomial GLM used the statsmodels logit link and `fit_regularized(alpha=1e-4, L1_wt=0)`, with an intercept. Missing primary-publication citation counts (5,556 datasets) were replaced with zero before the natural-log transformation log(1 + citations), with a separate missing-citation indicator. This indicator does not prove absence of a publication. Citation counts were the maximum iCite counts across linked publication records. Listed bytes were floored at 10⁶ bytes before log₁₀ transformation; missing values would be imputed at that floor. No size values were missing in the model cohort, and 25 fell below the floor. Missing research field was represented as “None” (4,460 datasets); the eight most frequent resulting categories were retained and others pooled as “Other”. The seven most frequent software-family categories were retained and the remainder pooled. Multiple named major tools formed a distinct group. Categorical reference levels were 2013, Mascot, and Agricultural and Biological Sciences. Undetected human and SDRF flags were false; missing group identities would receive accession-specific identifiers, although no such fallback occurred in this cohort. All 20,684 rows entered the 33-column design including the intercept.

**Intervals and estimands.** Logistic-model point estimates used all cohort rows. Two hundred cluster-bootstrap fits resampled whole contributor groups with replacement using the same fixed design; exponentiated coefficients were summarised by their 2.5th and 97.5th percentiles. The seed was 20261001; the source script shares its random-number stream with preceding analyses, so the execution order affects the exact bootstrap sequence. An offline rerun exactly reproduced the frozen results: all 201 fits had finite coefficients, no fit warnings and maximum absolute penalised scores below 1.0 × 10⁻⁵. Interval tails have limited Monte Carlo precision with 200 draws; the Thermo Fisher-only lower endpoint crossed one in a split-half diagnostic. Figure 5d instead shows rates directly standardised over release-year × human-status strata to the full cohort distribution, with 1,000 group-bootstrap draws. Empty strata within a subset or replicate were omitted and weights renormalised over represented strata. Figure 5e shows adjusted odds ratios. Sensitivity rate ratios are unadjusted feature-positive versus feature-negative contrasts with 1,000 group-bootstrap draws.

**Descriptive diversity.** Contributor-group effective numbers are reported without intervals because the original bootstrap retained duplicated group identities as the same category, changing category concentration. Other category axes retain the frozen group-bootstrap intervals. Rarefaction uses the frozen pooled figure-data values; dataset/category ordering is now deterministic for subsequent reruns. The instrument richness comparison rounds to 61 versus 66.

## Supplementary information

**Supplementary Figures.** 1, Concentration of datasets and bytes across contributor groups. 2, History of every common instrument model. 3, First appearance of new instrument models by continent. 4, Software release lines in use. 5, Age of software releases and sequence databases. 6, Composition of the archive's bytes. 7, File formats. 8, Isobaric labelling and post-translational modifications. 9, Organisms. 10, Countries. 11, International co-authorship. 12, Research fields of linked papers. 13, Timing of public release. 14, Submission timestamps.

**Supplementary Tables S1 to S21.** Provided in Supplementary-Tables.xlsx; the contents sheet gives the definition and source of each table.


The Supplementary Methods document mention ascertainment, model construction and inferential limits.

## Data and code availability

Code, derived inputs and figure evidence are at [github.com/filiprumenovski/pride-data-economy](https://github.com/filiprumenovski/pride-data-economy). Offline checks reproduce the headline run-count bootstrap, adjusted mention analysis and selected summaries. Contributor labels are pseudonymised; the identity crosswalk is excluded. Code uses MIT; author-created data, text and figures use CC BY 4.0, retaining third-party rights. The Zenodo archive is pending. Three extraction caches await transfer from local storage, and upstream geography reconstruction is incomplete in this package. Full reconstruction from raw records is not claimed; the repository documents coverage and limits.

## Author declarations

**Author.** Filip Rumenovski, Department of Chemistry, Wayne State University.

**Correspondence.** [filip.rumenovski@wayne.edu](mailto:filip.rumenovski@wayne.edu)

**Funding.** This research received no specific funding.

**Competing interests.** The author declares no competing interests.

**Author contributions.** Pending confirmation of the contribution statement and final approval.

## References

[PRIDE2025] Perez-Riverol, Y. et al. The PRIDE database at 20 years: 2025 update. *Nucleic Acids Res.* **53**, D543–D553 (2025). doi:10.1093/nar/gkae1011

[PX2026] Deutsch, E. W. et al. The ProteomeXchange consortium in 2026: making proteomics data FAIR. *Nucleic Acids Res.* **54**, D459–D469 (2026). doi:10.1093/nar/gkaf1146

[PSI2023] Deutsch, E. W. et al. Proteomics Standards Initiative at Twenty Years: Current Activities and Future Work. *J. Proteome Res.* **22**, 287–301 (2023). doi:10.1021/acs.jproteome.2c00637

[Geyer2017] Geyer, P. E., Holdt, L. M., Teupser, D. & Mann, M. Revisiting biomarker discovery by plasma proteomics. *Mol. Syst. Biol.* **13**, 942 (2017). doi:10.15252/msb.20156297

[Messner2023] Messner, C. B. et al. Mass spectrometry-based high-throughput proteomics and its role in biomedical studies and systems biology. *Proteomics* **23**, e2200013 (2023). doi:10.1002/pmic.202200013

[PerezRiverol2019] Perez-Riverol, Y. et al. Quantifying the impact of public omics data. *Nat. Commun.* **10**, 3512 (2019). doi:10.1038/s41467-019-11461-w

[Vadadokhau2026] Vadadokhau, U. et al. Preventing Proteomics Data Tombs Through Collective Responsibility and Community Engagement. *Sci. Data* **13**, 287 (2026). doi:10.1038/s41597-026-06614-8

[Leite2026] Leite, G. G. F., Carregari, V. C., Tashima, A. K. & Salomão, R. How Reusable Are Public Proteomics Datasets in Practice? *Proteomics* Early View, e70175 (2026). doi:10.1002/pmic.70175

[Meier2021] Meier, F., Park, M. A. & Mann, M. Trapped Ion Mobility Spectrometry and Parallel Accumulation–Serial Fragmentation in Proteomics. *Mol. Cell. Proteomics* **20**, 100138 (2021). doi:10.1016/j.mcpro.2021.100138

[Orsburn2021] Orsburn, B. C. Proteome Discoverer—A Community Enhanced Data Processing Suite for Protein Informatics. *Proteomes* **9**, 15 (2021). doi:10.3390/proteomes9010015

[Walzer2022] Walzer, M. et al. Implementing the reuse of public DIA proteomics datasets: from the PRIDE database to Expression Atlas. *Sci. Data* **9**, 335 (2022). doi:10.1038/s41597-022-01380-9

[VanDenBossche2026] Van Den Bossche, T. et al. Unlocking the Next Decade of Proteomics with Standardized, Structured Metadata. *J. Proteome Res.* **25**, 556–561 (2026). doi:10.1021/acs.jproteome.5c00958

[Dai2021] Dai, C. et al. A proteomics sample metadata representation for multiomics integration and big data analysis. *Nat. Commun.* **12**, 5854 (2021). doi:10.1038/s41467-021-26111-3

[Gatto2023] Gatto, L. et al. Initial recommendations for performing, benchmarking and reporting single-cell proteomics experiments. *Nat. Methods* **20**, 375–386 (2023). doi:10.1038/s41592-023-01785-3
