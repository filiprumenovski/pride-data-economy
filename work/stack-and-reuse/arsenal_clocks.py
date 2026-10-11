"""The clocks inside a deposit (descriptive).

Every dataset carries several dates:
  acquisition   latest acquisition date written in raw file names (Arsenal/H2-time; 34% of datasets)
  database      date written in the name of the main FASTA file, high confidence (Arsenal/H4; ~2.6k datasets)
  software      first submission in this archive that used the same software release (Arsenal/H9;
                a lower bound on how long the release had been available)
  submission, paper (OpenAlex publication date of the linked paper), public release (PRIDE).
For each clock: age at submission and at public release, and the order of events within datasets
that carry several clocks. Medians with group-bootstrap 95% intervals (pi_group), seed 20261001.
Outputs: Arsenal/H11-clocks/.
"""
import json
import numpy as np, pandas as pd
from common import R, O, SEED
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Arsenal/H11-clocks'; OUT.mkdir(parents=True, exist_ok=True)
A = O / 'Arsenal'
rng = np.random.default_rng(SEED)

F = pd.read_parquet(O / 'Deposit-features.parquet', columns=['accession', 'submission_date', 'release_date', 'work_id', 'release_year', 'pi_group'])
F['submission_date'] = pd.to_datetime(F.submission_date); F['release_date'] = pd.to_datetime(F.release_date)
W = connect().execute('select work_id, publication_date from lit.openalex_works').df()
W['paper_date'] = pd.to_datetime(W.publication_date, errors='coerce')
F = F.merge(W[['work_id', 'paper_date']], on='work_id', how='left')
T = pd.read_csv(A / 'H2-time/Dataset-time.csv.gz', usecols=['accession', 'latest_date', 'analysed'])
T = T[T.analysed].assign(acq_date=lambda x: pd.to_datetime(x.latest_date))[['accession', 'acq_date']]
B = pd.read_csv(A / 'H4-H7-databases-storage/H4-dataset-database-date.csv.gz')
B = B[B.date_conf == 'high'].assign(db_date=lambda x: pd.to_datetime(x.name_date))[['accession', 'db_date']]
V = pd.read_csv(A / 'H9-versions/Dataset-tool-version.csv.gz', usecols=['accession', 'software', 'release', 'first_seen'])
V['sw_date'] = pd.to_datetime(V.first_seen)
V1 = V.sort_values('sw_date').groupby('accession').head(1)[['accession', 'software', 'sw_date']]   # oldest native tool release
D = F.merge(T, on='accession', how='left').merge(B, on='accession', how='left').merge(V1, on='accession', how='left')
D['grp'] = D.pi_group.fillna('solo:' + D.accession)
days = lambda a, b: (D[a] - D[b]).dt.days
D['acq_to_sub'] = days('submission_date', 'acq_date'); D['acq_to_rel'] = days('release_date', 'acq_date')
D['db_to_sub'] = days('submission_date', 'db_date'); D['db_to_rel'] = days('release_date', 'db_date')
D['sw_to_sub'] = days('submission_date', 'sw_date'); D['sw_to_rel'] = days('release_date', 'sw_date')
D['sub_to_paper'] = days('paper_date', 'submission_date'); D['paper_to_rel'] = days('release_date', 'paper_date')
D['sub_to_rel'] = days('release_date', 'submission_date')
D['db_minus_acq'] = days('db_date', 'acq_date'); D['sw_minus_acq'] = days('sw_date', 'acq_date')
# drop impossible orderings (a clock dated after submission by >31 days; release before submission)
for c in ['acq_to_sub', 'db_to_sub', 'sw_to_sub']:
    D.loc[D[c] < -31, c.replace('_sub', '_rel')] = np.nan; D.loc[D[c] < -31, c] = np.nan
D.loc[D.sub_to_rel < 0, ['sub_to_rel', 'paper_to_rel']] = np.nan
D.to_csv(OUT / 'Dataset-clocks.csv.gz', index=False)


def med_ci(col, df=None, reps=1000):
    x = (D if df is None else df)[[col, 'grp']].dropna()
    gi = x.groupby('grp').indices; gs = list(gi); v = x[col].to_numpy()
    bs = [np.median(v[np.concatenate([gi[gs[k]] for k in rng.integers(0, len(gs), len(gs))])]) for _ in range(reps)]
    return {'datasets': int(len(x)), 'groups': int(len(gs)), 'median_days': float(np.median(v)), 'median_years': round(float(np.median(v)) / 365.25, 2),
            'ci_years': [round(np.percentile(bs, 2.5) / 365.25, 2), round(np.percentile(bs, 97.5) / 365.25, 2)],
            'share_gt_1y': round(float((v > 365.25).mean()), 3), 'share_gt_2y': round(float((v > 730.5).mean()), 3)}


rows = {c: med_ci(c) for c in ['acq_to_sub', 'acq_to_rel', 'db_to_sub', 'db_to_rel', 'sw_to_sub', 'sw_to_rel',
                               'sub_to_paper', 'paper_to_rel', 'sub_to_rel']}
by_tool = {t: {c: med_ci(c, D[D.software == t]) for c in ['sw_to_sub', 'sw_to_rel']} for t in ['MaxQuant', 'DIA-NN', 'FragPipe', 'Skyline', 'Mascot']}
# order of events in datasets with several clocks
both_db = D.dropna(subset=['acq_to_sub', 'db_minus_acq']); both_sw = D.dropna(subset=['acq_to_sub', 'sw_minus_acq'])
coh = {'acquisition_and_database': {'datasets': int(len(both_db)),
                                    'median_database_minus_acquisition_days': float(both_db.db_minus_acq.median()),
                                    'share_database_dated_before_latest_acquisition': round(float((both_db.db_minus_acq <= 0).mean()), 3),
                                    'median_acq_to_sub_days_same_datasets': float(both_db.acq_to_sub.median())},
       'acquisition_and_software': {'datasets': int(len(both_sw)),
                                    'median_software_first_seen_minus_acquisition_days': float(both_sw.sw_minus_acq.median()),
                                    'share_software_release_first_seen_before_latest_acquisition': round(float((both_sw.sw_minus_acq <= 0).mean()), 3),
                                    'spearman_data_age_vs_software_age': round(float(both_sw[['acq_to_sub', 'sw_to_sub']].corr('spearman').iloc[0, 1]), 3)}}
three = D.dropna(subset=['acq_to_rel', 'db_to_rel', 'sw_to_rel'])
summary = {'clocks': rows, 'software_by_tool': by_tool, 'order_of_events': coh, 'datasets_with_all_three_clocks': int(len(three)),
           'denominators': {'cohort': int(len(D)), 'with_paper_date': int(D.paper_date.notna().sum()), 'with_acquisition_date': int(D.acq_date.notna().sum()),
                            'with_database_date': int(D.db_date.notna().sum()), 'with_software_release': int(D.sw_date.notna().sum())}}
(OUT / 'Summary.json').write_text(json.dumps(summary, indent=2, default=str))
print(json.dumps(summary, indent=1, default=str))
