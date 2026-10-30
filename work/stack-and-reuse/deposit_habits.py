"""Is the way results are deposited a habit of the lab, and who drives the change over time?

Results form of each dataset (first match wins; archives opened, see archive_summary.py; table
headers read, see header_sniff.py, when available):
  standard                    mzIdentML/mzTab anywhere (top level or inside an archive)
  tool output, loose          known tool output among the top-level files
  tool output, zipped only    known tool output only inside archives
  tool output, renamed        no known file name, but a table header (top level, or inside a tables-only
                              zip: arsenal_zip_headers.py) identifies a tool output
  ad hoc tables only          tables of no identifiable origin, or results archives we could not read
  no results files            raw data / spectra only

1. Persistence. For consecutive datasets of the same contributor group (hashed `pi_group`,
   ordered by submission date), the observed share whose form equals the previous dataset's form
   is compared with the share expected if the form were drawn from all OTHER groups' datasets in
   the same stratum (release year x software family x DIA). Ratio > 1 = group-specific practice.
   Stratified by days between deposits, same vs different paper, and same vs switched software.
   Group-bootstrap 95% intervals (1,000 resamples of groups). Descriptive.
2. Turnover. Recognisable-results rate for a group's first-year datasets vs returning groups, and
   by entry cohort; exact shift-share between 2015-17 and 2023-25 (composition term centred).
"""
import json
import numpy as np, pandas as pd
from common import R, O, SEED, features
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect
from layer_rules import RULES

OUT = O / 'Habits'; OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(SEED)

D = pd.read_csv(O / 'Archive-contents/Dataset-before-after.csv.gz').set_index('accession')
case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
loose = set(connect().execute(f"""
with f as (select project_accession a, lower(file_name) n, lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
  coalesce(file_category_value, '') cat from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
  and project_accession in (select accession from study))
select distinct a from f where ({case}) = 'tool_results' and not regexp_matches(n, '\\.(zip|tgz|tar|rar|7z)$|\\.tar\\.gz$')""").df().a)
D['loose_tool'] = D.index.isin(loose)
H = O / 'Header-sniff/Dataset-headers.csv.gz'
D['header_tool'] = pd.read_csv(H).set_index('accession').any_identified.reindex(D.index).fillna(False).astype(bool) if H.exists() else False
HZ = O / 'Arsenal/H10-zip-headers/Dataset-zip-headers.csv'   # headers of tables inside tables-only zips
if HZ.exists():
    D['header_tool'] |= pd.read_csv(HZ).set_index('accession').zip_header_tool.reindex(D.index).fillna(False).astype(bool)
D['form'] = np.select([D.standard_after, D.loose_tool, D.recognisable_after, D.adhoc_only_after & D.header_tool, D.adhoc_only_after],
                      ['standard', 'tool output, loose', 'tool output, zipped only', 'tool output, renamed', 'ad hoc tables only'],
                      'no results files')
D['recognisable'] = D.form.isin(['standard', 'tool output, loose', 'tool output, zipped only', 'tool output, renamed'])
F = features().set_index('accession')
D = D.join(F[['pi_group', 'group_first_year', 'family', 'submission_date', 'work_id']])
D[['release_year', 'form', 'recognisable', 'pi_group', 'family', 'dia_expanded']].to_csv(OUT / 'Dataset-form.csv.gz')
forms = pd.crosstab(D.release_year, D.form, normalize='index').mul(100).round(1)
forms['datasets'] = D.groupby('release_year').size(); forms.to_csv(OUT / 'Form-by-year.csv')

# ---- 1. persistence ------------------------------------------------------------------------
G = D[D.pi_group.notna()].reset_index().sort_values(['pi_group', 'submission_date', 'accession'])
G['stratum'] = G.release_year.astype(str) + '|' + G.family + '|' + G.dia_expanded.astype(str)
cnt = G.groupby(['stratum', 'form']).size(); gcnt = G.groupby(['pi_group', 'stratum', 'form']).size()
tot = G.groupby('stratum').size(); gtot = G.groupby(['pi_group', 'stratum']).size()
for c in ['form', 'family', 'dia_expanded', 'submission_date', 'work_id', 'accession']:
    G['prev_' + c] = G.groupby('pi_group')[c].shift()
P = G[G.prev_form.notna()].copy()
num = np.array([cnt.get((s, f), 0) - gcnt.get((g, s, f), 0) for g, s, f in zip(P.pi_group, P.stratum, P.prev_form)], float)
den = np.array([tot[s] - gtot[(g, s)] for g, s in zip(P.pi_group, P.stratum)], float)
P['expected'] = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
P['observed'] = (P.form == P.prev_form).astype(float)
P = P[P.expected.notna()]
P['gap_days'] = (pd.to_datetime(P.submission_date) - pd.to_datetime(P.prev_submission_date)).dt.days
P['same_paper'] = P.work_id.notna() & (P.work_id == P.prev_work_id)
P['software'] = np.where(P.family == P.prev_family, 'same software family', 'switched software family')
P['mode'] = np.where(P.dia_expanded == P.prev_dia_expanded, 'same DIA/DDA', 'switched DIA/DDA')
P['gap'] = pd.cut(P.gap_days, [-1, 30, 180, 365, 730, 1460, 100000],
                  labels=['<=30 d', '31-180 d', '181-365 d', '1-2 y', '2-4 y', '>4 y']).astype(str)
idx = P.groupby('pi_group').indices; groups = np.array(list(idx))


def ratio_ci(x):
    gi = x.groupby('pi_group').indices; gs = np.array(list(gi)); o, e = x.observed.to_numpy(), x.expected.to_numpy()
    so = np.array([o[gi[g]].sum() for g in gs]); se = np.array([e[gi[g]].sum() for g in gs])
    bs = [so[s].sum() / se[s].sum() for s in (rng.integers(0, len(gs), len(gs)) for _ in range(1000))]
    return pd.Series({'pairs': len(x), 'groups': len(gs), 'observed_pct': 100 * o.mean(), 'expected_pct': 100 * e.mean(),
                      'ratio': o.sum() / e.sum(), 'ratio_lo': np.percentile(bs, 2.5), 'ratio_hi': np.percentile(bs, 97.5)})


rows = [ratio_ci(P).rename('all pairs')]
for col in ['gap', 'software', 'mode', 'prev_form']:
    for k, x in P.groupby(col):
        rows.append(ratio_ci(x).rename(f'{col}: {k}'))
strict = P[(P.gap_days >= 365) & ~P.same_paper & (P.software == 'switched software family')]
rows.append(ratio_ci(P[~P.same_paper]).rename('different paper'))
rows.append(ratio_ci(strict).rename('>=1 y apart, different paper, switched software'))
pers = pd.DataFrame(rows).round(3); pers.to_csv(OUT / 'Persistence.csv')
trans = pd.crosstab(P.prev_form, P.form, normalize='index').mul(100).round(1); trans.to_csv(OUT / 'Transition-matrix.csv')

# ---- 2. turnover ---------------------------------------------------------------------------
D['who'] = np.where(D.pi_group.isna(), 'no group id', np.where(D.group_first_year == D.release_year, 'first-year group', 'returning group'))
turn = D.groupby(['release_year', 'who']).agg(datasets=('recognisable', 'size'), recognisable_pct=('recognisable', 'mean'))
turn['recognisable_pct'] *= 100; turn.round(1).to_csv(OUT / 'Turnover-by-year.csv')
cuts, labs = [2012, 2016, 2019, 2022, 2025], ['2013-16', '2017-19', '2020-22', '2023-25']
D['period'] = pd.cut(D.release_year, cuts, labels=labs); D['cohort'] = pd.cut(D.group_first_year, cuts, labels=['entered ' + l for l in labs])
coh = D.groupby(['cohort', 'period', 'dia_expanded'], observed=True).agg(datasets=('recognisable', 'size'), recognisable_pct=('recognisable', 'mean'))
coh['recognisable_pct'] *= 100; coh.round(1).to_csv(OUT / 'Cohort-by-period.csv')
fam = D[~D.dia_expanded].groupby(['family', 'period'], observed=True).agg(datasets=('recognisable', 'size'), recognisable_pct=('recognisable', 'mean'))
fam['recognisable_pct'] *= 100; fam.round(1).to_csv(OUT / 'Non-DIA-by-software-period.csv')
a, b = D[D.release_year.between(2015, 2017) & D.pi_group.notna()], D[D.release_year.between(2023, 2025) & D.pi_group.notna()]
sa, sb = a.who.value_counts(normalize=True), b.who.value_counts(normalize=True)
ra, rb = a.groupby('who').recognisable.mean(), b.groupby('who').recognisable.mean(); cbar = D.recognisable.mean()
shift = {'total_pp': 100 * (b.recognisable.mean() - a.recognisable.mean()),
         'within_groups_pp': 100 * ((sa + sb) / 2 * (rb - ra)).sum(),
         'first_year_vs_returning_mix_pp': 100 * ((sb - sa) * ((ra + rb) / 2 - cbar)).sum(),
         'first_year_share_2015_17': sa.get('first-year group'), 'first_year_share_2023_25': sb.get('first-year group')}
summary = {'datasets': int(len(D)), 'with_group_id': int(D.pi_group.notna().sum()), 'pairs': int(len(P)),
           'header_results_used': bool(H.exists()), 'persistence_all': pers.loc['all pairs'].to_dict(),
           'persistence_strict': pers.loc['>=1 y apart, different paper, switched software'].to_dict(), 'turnover_shift_share': shift}
(OUT / 'Summary.json').write_text(json.dumps(summary, indent=2, default=float))
pd.set_option('display.width', 250)
print(json.dumps(summary, indent=1, default=float)); print(pers.to_string()); print(forms.to_string())
