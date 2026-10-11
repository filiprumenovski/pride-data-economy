"""Supplementary Tables for the paper, assembled from frozen outputs (no new estimates except Table S1 counts).

Writes outputs/Stack-and-reuse/Supplementary-Tables/Supplementary-Tables.xlsx (one sheet per table plus a
contents sheet with each table's title, definition and source files) and one CSV per table.
"""
import json, re, sys
import pandas as pd
from common import O, features

OUT = O / 'Supplementary-Tables'; OUT.mkdir(exist_ok=True)
GD = O / 'Gallery/data'; A = O / 'Arsenal'
T = []   # (id, title, note, sources, frame)


def add(tid, title, note, sources, df):
    T.append((tid, title, note, sources, df.reset_index(drop=True)))


# S1 cohort and denominators
F = features()
H = pd.read_csv(O / 'Habits/Dataset-form.csv.gz', usecols=['accession', 'recognisable'])
F = F.merge(H, on='accession', how='left')
F['first_year'] = F.group_first_year.eq(F.release_year)
g = F.groupby('release_year')
s1 = pd.DataFrame({
    'datasets': g.size(), 'contributor_groups': g.pi_group.nunique(),
    'country_resolved_pct': 100 * g.country.apply(lambda x: x.notna().mean()),
    'linked_paper_pct': 100 * g.work_id.apply(lambda x: x.notna().mean()),
    'from_first_year_group_pct': 100 * g.first_year.mean(),
    'dia_pct': 100 * g.dia_expanded.mean(),
    'thermo_single_vendor_pct': 100 * g.vendor_single.apply(lambda x: x.eq('Thermo Fisher').mean()),
    'bruker_single_vendor_pct': 100 * g.vendor_single.apply(lambda x: x.eq('Bruker').mean()),
    'standard_result_files_pct': 100 * g.standards.mean(),
    'identifiable_results_pct': 100 * g.recognisable.apply(lambda x: x.fillna(False).astype(bool).mean()),
    'human_pct': 100 * g.human.apply(lambda x: x.fillna(False).astype(bool).mean()),
}).round(1).reset_index()
add('S1', 'Cohort and denominators by release year',
    'All 36,249 PRIDE datasets released 2013–2025. Standard result files = corrected non-generated expanded mzIdentML/mzTab definition. '
    'Identifiable results = file names, archive members and table headers (lower bound). Contributor group = hashed PI contact.',
    'Deposit-features.parquet; Habits/Dataset-form.csv.gz', s1)

# S2 effective numbers by axis
b = pd.read_csv(GD / '42-breadth.csv')
# Contributor identity is also the diversity category: ordinary cluster resampling
# repeats category labels and depresses effective numbers. Report this axis descriptively.
b.loc[b.axis.eq('Contributor groups'), ['lo', 'hi']] = float('nan')
add('S2', 'Effective number of categories by axis and release year (Fig. 1a)',
    'Inverse Simpson index of category shares among each year’s datasets; fractional credit for multiple categories; 95% intervals by bootstrap over contributor groups (300 replicates). Contributor-group diversity is descriptive; its intervals are omitted.',
    'Gallery/data/42-breadth.csv', b.round(3))
pooled_breadth = pd.read_csv(GD / '42-breadth-pooled.csv').round(3)
pooled_breadth.loc[pooled_breadth.axis.eq('Contributor groups'), ['q2_lo', 'q2_hi', 'q1_lo', 'q1_hi']] = float('nan')
add('S2b', 'Breadth of the archive: 2023–25 vs 2013–15 at three diversity orders (Fig. 1b)',
    'Hill numbers: order 2 = inverse Simpson, order 1 = exp(Shannon), order 0 = richness rarefied to the smaller period’s dataset count (mean of 50 '
    'subsamples). Ratios with 95% intervals from bootstrap over contributor groups within each period (orders 1 and 2; 300 replicates). '
    'Instrument models harmonised as in the fleet analysis; vendors are the harmonised models’ vendors; contributor groups credit each dataset to its group. Contributor-group diversity is descriptive; its intervals are omitted.',
    'Gallery/data/42-breadth-pooled.csv', pooled_breadth)

# S3 concentration of supply
c1 = pd.read_csv(GD / '34-contributor-concentration.csv'); c2 = pd.read_csv(GD / '34-contributor-concentration-classes.csv')
add('S3', 'Concentration of datasets and bytes across contributor groups (Supplementary Fig. 1)',
    'Shares supplied by the largest groups, by period; Gini coefficients; 95% bootstrap over groups.',
    'Gallery/data/34-contributor-concentration*.csv', c1.round(4))
add('S3b', 'Datasets, groups and bytes by group size class (Supplementary Fig. 1)', '', 'Gallery/data/34-contributor-concentration-classes.csv', c2.round(4))

# S4 countries
w = pd.read_csv(GD / '20-world-growth.csv')[['country', 'continent', 'd_2016_2018', 'd_2023_2025', 'added', 'multiple', 'multiple_lo', 'multiple_hi']]
col = pd.read_csv(GD / '35-collaboration-by-country.csv')[['country', 'papers', 'est', 'lo', 'hi']].rename(
    columns={'papers': 'linked_papers', 'est': 'international_share', 'lo': 'intl_lo', 'hi': 'intl_hi'})
age = pd.read_csv(GD / '39-data-age-by-country.csv')[['country', 'll_datasets', 'll_median_days', 'll_ci_lo', 'll_ci_hi', 'coverage']].rename(
    columns={'ll_datasets': 'data_age_datasets', 'll_median_days': 'data_age_median_days', 'll_ci_lo': 'age_lo', 'll_ci_hi': 'age_hi',
             'coverage': 'dated_raw_share'})
age['country'] = age.country.replace({'UK': 'United Kingdom', 'USA': 'United States'})
s4 = w.merge(col, on='country', how='left').merge(age, on='country', how='left').sort_values('d_2023_2025', ascending=False)
add('S4', 'Countries: growth, international co-authorship and data age (Fig. 1e–g; Supplementary Figs 10–11)',
    'Datasets by contributing group’s country (country-resolved datasets only). Growth multiple 2016–18 → 2023–25 with 95% group-bootstrap interval. '
    'International share = linked papers with author institutions in ≥2 countries. Data age = days from latest acquisition date in raw-file names to submission, '
    'Thermo-only non-DIA datasets submitted 2016–2025.', 'Gallery/data/20-world-growth.csv; 35-collaboration-by-country.csv; 39-data-age-by-country.csv',
    s4.round(4))

# S5 instrument models
lc = pd.read_csv(A / 'H1-H3-scale-instruments/H3-model-life-cycles.csv')[['model', 'vendor', 'datasets', 'groups_total', 'first_year_any', 'first_year_ge5',
                                                                           'peak_year_share', 'peak_share', 'share_2025', 'left_censored']]
rs = pd.read_csv(A / 'H1-H3-scale-instruments/H1-per-run-size-by-model.csv')[['model', 'datasets', 'groups', 'median_runs', 'median_per_run_GB', 'p25', 'p75', 'dia_share']].rename(
    columns={'datasets': 'run_size_datasets', 'groups': 'run_size_groups'})
s5 = lc.merge(rs, on='model', how='left').sort_values('datasets', ascending=False)
add('S5', 'Instrument models: datasets, life cycle and run size (Figs 2d, 3a; Supplementary Figs 2–3)',
    'Models named in ≥50 datasets. peak_share and share_2025 are unweighted proportions of model-labelled datasets naming the model; '
    'a dataset naming several models counts once for each. peak_year_share is the year of the largest unweighted proportion. '
    'These differ from the fractional-credit shares in Fig. 3a and Supplementary Fig. 2. first_year_any and first_year_ge5 use ≥1 and ≥5 datasets, respectively; '
    'Supplementary Fig. 2 orders models by the first year with ≥3 datasets. datasets/groups_total count all deposits/groups naming the model. '
    'run_size_datasets/run_size_groups count the single-model subset with exact run counts and positive resolved bytes per run; summaries require ≥50 eligible datasets. '
    'Within this subset, median_runs is the median exact run count; median_per_run_GB and p25/p75 are the median and 25th/75th percentiles of dataset-level '
    'resolved bytes divided by exact run count, in GB (10^9 bytes); dia_share is the fraction with DIA detected. '
    'Blank run-size fields mean the model did not meet the subset threshold. A deposited label is not an ownership record.', 'Arsenal/H1-H3-scale-instruments/H3-model-life-cycles.csv; H1-per-run-size-by-model.csv', s5.round(4))

# S6 runs-per-dataset robustness
add('S6', 'Raw runs per dataset, 2013–16 vs 2023–25, under alternative definitions (Fig. 2c)',
    'Median ratio and change in the share with ≥100 runs, with 95% group-bootstrap intervals.',
    'Attacks/runs/Contrasts-2013-16-vs-2023-25.csv', pd.read_csv(O / 'Attacks/runs/Contrasts-2013-16-vs-2023-25.csv').round(4))

# S7 software
st = pd.read_csv(O / 'Tool-outputs/Standard-rate-by-tool-single-tool-datasets.csv')
ex = pd.read_csv(O / 'Tool-outputs/Standard-rate-by-export-capability.csv')
add('S7', 'Standard result files by analysis tool and export capability (Fig. 4a–b)',
    'Single-tool datasets. Export capability from the software registry (Software-registry.md): default = tool writes mzIdentML/mzTab by default; built-in = option; none = cannot.',
    'Tool-outputs/*.csv; Software-registry.md', st.round(2))
add('S7b', 'Standard result files by export capability, with 95% group-bootstrap intervals (Fig. 4b)', '', 'Tool-outputs/Standard-rate-by-export-capability.csv', ex.round(2))
add('S8', 'Age of the software release at public release, by tool (Supplementary Fig. 5)',
    'Years from the first PRIDE submission recording the exact release to public release of the dataset (a lower bound); medians with group-bootstrap intervals; '
    'ceilings = oldest release of the tool in PRIDE.', 'Gallery/data/33-tool-age-at-release.csv', pd.read_csv(GD / '33-tool-age-at-release.csv').round(3))
add('S9', 'Software release lines in 2025 submissions (Supplementary Fig. 4)', 'Share of datasets with a recovered version.',
    'Gallery/data/31-version-streams-2025-shares.csv', pd.read_csv(GD / '31-version-streams-2025-shares.csv').round(3))

# S10 results by reading depth
ba = pd.read_csv(O / 'Archive-contents/Before-after-by-year.csv')
hh = pd.read_csv(O / 'Header-sniff/Recognisable-by-year-three-ways.csv')[['release_year', 'recognisable_with_headers', 'adhoc_after_headers']]
add('S10', 'What travels with the data, by release year and reading depth (Fig. 4c–d)',
    'Percent of datasets. before = top-level file names; after = adding archive members; with_headers = adding table headers.',
    'Archive-contents/Before-after-by-year.csv; Header-sniff/Recognisable-by-year-three-ways.csv', ba.merge(hh, on='release_year'))
add('S11', 'Documentation persistence within contributor groups (Fig. 4e)',
    'Next dataset of the same group carrying the file, by whether the previous dataset did, vs the expected rate for datasets of the same years.',
    'Arsenal/H5-H6-H8-strategies-software-docs/H8-documentation-persistence-within-groups.csv',
    pd.read_csv(A / 'H5-H6-H8-strategies-software-docs/H8-documentation-persistence-within-groups.csv'))

# S12 clocks and release timing
ck = json.load(open(A / 'H11-clocks/Summary.json'))['clocks']
s12 = pd.DataFrame([{'interval': k, **{kk: (str(v) if isinstance(v, list) else v) for kk, v in d.items()}} for k, d in ck.items()])
add('S12', 'Clocks of a deposit: medians with 95% group-bootstrap intervals (Fig. 5a)',
    'acq = latest acquisition date in raw file names; db = date in the main FASTA name; sw = first PRIDE submission recording the software release (lower bound); '
    'sub = submission; rel = public release; paper = OpenAlex publication date.', 'Arsenal/H11-clocks/Summary.json', s12)
add('S13', 'Release relative to the linked paper, bulk release days vs other days (Supplementary Fig. 13)', '',
    'Gallery/data/41-paper-and-release-summary.csv', pd.read_csv(GD / '41-paper-and-release-summary.csv'))

# S14 journals
add('S14', 'Journals of dataset-linked papers (Fig. 5b–c)', 'Top 15 journals by linked datasets and share of each year’s linked papers.',
    'Gallery/data/37-journals-top15.csv; 37-journals-concentration.csv', pd.read_csv(GD / '37-journals-top15.csv'))
add('S14b', 'Concentration of dataset-linked papers across journals by year (Fig. 5c)', '', 'Gallery/data/37-journals-concentration.csv',
    pd.read_csv(GD / '37-journals-concentration.csv').round(4))

# S15 reuse
rr = pd.read_csv(O / 'Reuse-mentions/Reuse-rates.csv')
rr = rr.replace('Thermo raw', 'Thermo Fisher only')
add('S15', 'Non-primary cross-group accession mentions: standardised rates by feature (Fig. 5d)',
    'D1 = non-primary mention in release year or next three calendar years by a paper naming 2 or more contributor groups and at most 100 datasets. Same-year chronology is unresolved. '
    'D2 = publication-year lag at least 2, excluding primary papers linked to another deposit of the same group; D3 = both. D1_all removes the 100-dataset ceiling. '
    'Rates directly standardised over release-year by human-status strata, with 1000 group-bootstrap draws; empty strata omitted and weights renormalised. Thermo Fisher only denotes instrument-vendor assignment.',
    'Reuse-mentions/Reuse-rates.csv', rr)
R = json.load(open(O / 'Attacks/reuse-mentions/Results.json'))
rob = [{'test': 'Adjusted odds ratio (A3)', 'feature': 'Thermo Fisher only' if k == 'thermo' else k, 'estimate': v[0], 'lo': v[1], 'hi': v[2]} for k, v in R['A3_adjusted_odds_ratios'].items()]
for test in ['A2_excluding_method_like', 'A5_window_5y_release_to_2020']:
    for k, v in R[test].items():
        if isinstance(v, list):
            rob.append({'test': f'Crude rate ratio, {test}', 'feature': 'Thermo Fisher only' if k == 'thermo' else k, 'estimate': v[0], 'lo': v[1], 'hi': v[2]})
rob.append({'test': 'Crude rate ratio, A4 human; Thermo Fisher only; no DIA detected; MaxQuant/PD', 'feature': 'standard', 'estimate': R['A4_human_thermo_dda_mq_pd']['standard'][0],
            'lo': R['A4_human_thermo_dda_mq_pd']['standard'][1], 'hi': R['A4_human_thermo_dda_mq_pd']['standard'][2]})
add('S16', 'Reuse: adjusted odds ratios and robustness checks (Fig. 5e)',
    'Ridge binomial logit GLM with alpha 1e-4 and L1_wt 0; 200 contributor-group bootstrap draws. '
    'Predictors: standard files, detected DIA, SDRF, human and Thermo Fisher-only vendor assignment. Adjustment: release year, pooled field/software categories, log10 bytes floored at 1e6, log(1+citations) and missing-citation indicator. '
    'See Supplementary Methods for missingness and category reference levels. A4/A5 report crude rate ratios with 1000 group-bootstrap draws. A5 uses release year plus five calendar years among 2013-2020 releases. SDRF direction is not interpretable because files may be added after release.',
    'Attacks/reuse-mentions/Results.json', pd.DataFrame(rob))

# S17 findings ledger status
rows = []
for line in (O / 'Findings-ledger.md').read_text().splitlines():
    m = re.match(r'^\| ([A-Z]\d+(?:\.\d+)?[a-z]?) \| (.*?) \|', line)
    if m:
        cells = [x.strip() for x in line.strip('|').split('|')]
        rows.append({'id': cells[0], 'finding': cells[1], 'status': cells[-2] if len(cells) > 3 else ''})
# verdicts of the robustness attacks override ledger statuses
txt = (O / 'Findings-ledger.md').read_text()
sec = txt[txt.index('## Robustness attacks (verdicts)'):]
sec = sec[:sec.index('\n## ', 5)]
over = {}
for line in sec.splitlines():
    m = re.match(r'^\| ([^|]+?) \| ([^|]+?) \| (.*?) \|\s*$', line)
    if m and not m.group(1).startswith(('Finding', '---')):
        for tok in re.findall(r'[A-H]\d+(?:\.\d+)?', m.group(1)):
            over[tok] = (m.group(2).strip(), m.group(3).strip())
R20 = pd.DataFrame(rows)
R20['attack_verdict'] = R20.id.map(lambda i: over.get(i, ('', ''))[0])
R20['corrected_wording'] = R20.id.map(lambda i: over.get(i, ('', ''))[1])
# Reporting corrections supersede exploratory ledger wording; retain frozen source ledger.
reporting = {
 'G5': 'Qualifying non-primary cross-group mentions in release year or the next three calendar years; within-year chronology unresolved.',
 'G6': 'No detectable adjusted association between standard result files and qualifying mentions; a modest advantage remains possible.',
 'G7': 'Detected DIA, human samples and Thermo Fisher-only vendor assignments associated with qualifying mentions; mentions do not establish reanalysis.',
 'G8': 'SDRF associated with qualifying mentions; temporal direction unresolved because metadata may be added after release.',
 'H1.4': 'Mean log10 bytes per run decomposed in single-model, exact-count datasets: 2014-2016 versus 2023-2025. Model mix contributes 2.30-fold and within-model change 0.85-fold (0.80-1.10).',
}
for key, wording in reporting.items():
    R20.loc[R20.id.eq(key), ['finding', 'corrected_wording']] = wording
for col in ['finding', 'corrected_wording']:
    R20[col] = R20[col].str.replace('within DDA', 'among datasets with no DIA detected', regex=False).str.replace('Within DDA', 'Among datasets with no DIA detected', regex=False).str.replace('non-DIA', 'no DIA detected', regex=False)
    R20[col] = R20[col].str.replace('About 38% of the decline is within continuing groups.', 'Within-group log2 run-count slope is 38% of the pooled slope among multi-year groups; this is not a decomposition of the median decline.', regex=False)
    R20[col] = R20[col].str.replace("First-year groups' datasets are 0.64× the size of comparable returning groups'.", 'The matched geometric-mean run-count ratio is 0.64 for first-time versus returning groups among 2018-2025 releases.', regex=False)
add('S17', 'Findings ledger: every result with its robustness status',
    'Status: V verified; W weakened; S single analysis, not yet attacked; X dropped. attack_verdict and corrected_wording come from the robustness attacks and override status where present. Full numbers, denominators and caveats in Findings-ledger.md.',
    'Findings-ledger.md', R20)

T = [(f'S{i}', ti, n, s_, df) for i, (_, ti, n, s_, df) in enumerate(T, 1)]   # sequential numbering
for old in OUT.glob('Table-*.csv'):
    old.unlink()
contents = pd.DataFrame([{'table': t, 'title': ti, 'definition': n, 'source': s, 'rows': len(df)} for t, ti, n, s, df in T])
if '--csv-only' not in sys.argv:
    with pd.ExcelWriter(OUT / 'Supplementary-Tables.xlsx') as xw:
        contents.to_excel(xw, sheet_name='Contents', index=False)
        for t, ti, n, s, df in T:
            df.to_excel(xw, sheet_name=t, index=False)
for t, ti, n, s, df in T:
    df.to_csv(OUT / f'Table-{t}.csv', index=False)
contents.to_csv(OUT / 'Contents.csv', index=False)
print(contents[['table', 'title', 'rows']].to_string())
