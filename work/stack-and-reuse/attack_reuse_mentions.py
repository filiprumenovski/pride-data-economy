"""Attack the mention-based reuse contrasts (Reuse-mentions/).

Claims under attack (D1, released 2013-2022, window release year .. +3):
  standard result files are not associated with more cross-group reuse; DIA ~2x; SDRF ~8x; human ~1.8x.
Attacks:
  A1 paper concentration  leave-one-paper-out over the 25 papers with most links; number of distinct papers per contrast
  A2 benchmark/method     drop datasets whose title/description reads as a benchmark, method or standard sample
  A3 adjusted association logistic model with year, field, species, original-paper citations, size, vendor, software family;
                          odds ratios with group-bootstrap intervals (descriptive adjustment, not causal)
  A4 like-for-like        human, Thermo, DDA, MaxQuant or Proteome Discoverer only: standard vs not
  A5 window               release 2013-2020, mentions to +5 years
  A6 SDRF                 concentration of SDRF datasets and their reusing papers across groups
  A7 SDRF timing          datasets released 2013-2018 predate the SDRF format, so an SDRF there was added afterwards
Outputs: outputs/Stack-and-reuse/Attacks/reuse-mentions/.
"""
import json, re, sys
import numpy as np, pandas as pd
import statsmodels.api as sm
from common import O, R, SEED
sys.path.insert(0, str(R / 'scripts'))
from query import connect

IN = O / 'Reuse-mentions'; OUT = O / 'Attacks/reuse-mentions'; OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(SEED)
C = pd.read_csv(IN / 'Dataset-reuse.csv.gz', index_col=0)
M = pd.read_csv(IN / 'Mention-pairs.csv.gz')
F = pd.read_parquet(O / 'Deposit-features.parquet', columns=['accession', 'primary_field', 'primary_citations', 'listed_bytes', 'f_raw_files']).set_index('accession')
C = C.join(F)
CONTRASTS = ['standard', 'dia', 'sdrf', 'human', 'thermo']
res = {}


def ratio(df, col, out='D1'):
    a, b = df[df[col]][out].mean(), df[~df[col]][out].mean()
    return a / b if b > 0 else np.nan


def ratio_ci(df, col, out='D1', reps=1000):
    gi = df.groupby('grp').indices; gs = np.array(list(gi))
    bs = []
    for _ in range(reps):
        sub = df.iloc[np.concatenate([gi[g] for g in gs[rng.integers(0, len(gs), len(gs))]])]
        bs.append(ratio(sub, col, out))
    return [round(ratio(df, col, out), 2), round(np.nanpercentile(bs, 2.5), 2), round(np.nanpercentile(bs, 97.5), 2)]


# A1 paper concentration
win = M[M.lag.between(0, 3) & (M.kg >= 2) & (M.k <= 100) & M.accession.isin(C.index)]
top = win.groupby('pmid').accession.nunique().sort_values(ascending=False).head(25).index
a1 = {}
for col in CONTRASTS:
    ids = C.index[C[col]]
    a1[col] = {'papers_linking_feature_datasets': int(win[win.accession.isin(ids)].pmid.nunique()),
               'ratio': round(ratio(C, col), 2)}
    loo = []
    for p in top:
        keep = set(win[win.pmid != p].accession)
        loo.append(ratio(C.assign(D1=C.index.isin(keep)), col))
    a1[col]['leave_one_paper_out_range'] = [round(min(loo), 2), round(max(loo), 2)]
res['A1_paper_concentration'] = a1

# A2 benchmark / method datasets
P = connect().execute("select project_accession accession, project_title t, project_description d from pride_projects").df().set_index('accession')
pat = re.compile(r'benchmark|standard (?:mixture|sample|protein)|spike[- ]?in|lfqbench|ups1|ups2|hela (?:digest|standard)|'
                 r'method(?:s|ology)? (?:development|comparison|evaluation)|comparison of|evaluation of|workflow|technical|'
                 r'software|algorithm|tool|quality control|\bqc\b|instrument performance', re.I)
txt = (P.t.fillna('') + ' ' + P.d.fillna('')).reindex(C.index).fillna('')
C['method_like'] = txt.str.contains(pat)
C2 = C[~C.method_like]
res['A2_excluding_method_like'] = {'method_like_share': round(float(C.method_like.mean()), 3),
                                   'method_like_reuse_ratio': round(ratio(C, 'method_like'), 2),
                                   **{col: ratio_ci(C2, col) for col in CONTRASTS}}

# A3 adjusted logistic model, group bootstrap of odds ratios
X = pd.DataFrame({'standard': C.standard, 'dia': C.dia, 'sdrf': C.sdrf, 'human': C.human, 'thermo': C.thermo,
                  'log_citations': np.log1p(C.primary_citations.fillna(0)), 'no_paper': C.primary_citations.isna(),
                  'log_bytes': np.log10(C.listed_bytes.clip(lower=1e6).fillna(1e6))}).astype(float)
fam = C.family.where(C.family.isin(C.family.value_counts().head(7).index), 'Other')
field = C.primary_field.fillna('None'); field = field.where(field.isin(field.value_counts().head(8).index), 'Other')
X = pd.concat([X, pd.get_dummies(C.release_year.astype(str), prefix='y', drop_first=True, dtype=float),
               pd.get_dummies(fam, prefix='fam', drop_first=True, dtype=float),
               pd.get_dummies(field, prefix='field', drop_first=True, dtype=float)], axis=1)
X = sm.add_constant(X)
y = C.D1.astype(float).to_numpy(); Xa = X.to_numpy()


def fit_or(idx):
    m = sm.GLM(y[idx], Xa[idx], family=sm.families.Binomial()).fit_regularized(alpha=1e-4, L1_wt=0)
    return np.exp(np.asarray(m.params)[[X.columns.get_loc(k) for k in ['standard', 'dia', 'sdrf', 'human', 'thermo', 'log_citations']]])


est = fit_or(np.arange(len(C)))
gi = C.reset_index().groupby('grp').indices; gs = np.array(list(gi))
bs = np.array([fit_or(np.concatenate([gi[g] for g in gs[rng.integers(0, len(gs), len(gs))]])) for _ in range(200)])
res['A3_adjusted_odds_ratios'] = {k: [round(e, 2), round(lo, 2), round(hi, 2)] for k, e, lo, hi in
                                  zip(['standard', 'dia', 'sdrf', 'human', 'thermo', 'log_citations'], est,
                                      np.percentile(bs, 2.5, axis=0), np.percentile(bs, 97.5, axis=0))}

# A4 like-for-like
L4 = C[C.human & C.thermo & ~C.dia & C.family.isin(['MaxQuant', 'Proteome Discoverer'])]
res['A4_human_thermo_dda_mq_pd'] = {'n': int(len(L4)), 'standard': ratio_ci(L4, 'standard'),
                                    'rates_pct': [round(100 * L4[L4.standard].D1.mean(), 2), round(100 * L4[~L4.standard].D1.mean(), 2)]}

# A5 longer window
w5 = M[M.lag.between(0, 5) & (M.kg >= 2) & (M.k <= 100)]
C5 = C[C.release_year <= 2020].assign(D1=lambda x: x.index.isin(w5.accession))
res['A5_window_5y_release_to_2020'] = {'n': int(len(C5)), 'reused_pct': round(100 * C5.D1.mean(), 2),
                                       **{col: ratio_ci(C5, col) for col in CONTRASTS}}

# A6 SDRF concentration
S = C[C.sdrf]
res['A6_sdrf'] = {'datasets': int(len(S)), 'groups': int(S.grp.nunique()),
                  'top_group_share': round(float(S.grp.value_counts(normalize=True).iloc[0]), 3),
                  'reused': int(S.D1.sum()), 'reusing_papers': int(win[win.accession.isin(S.index)].pmid.nunique()),
                  'top_paper_links': win[win.accession.isin(S.index)].pmid.value_counts().head(5).to_dict(),
                  'sdrf_by_year': S.release_year.value_counts().sort_index().to_dict(),
                  'excluding_top_group': ratio_ci(C[C.grp != S.grp.value_counts().index[0]], 'sdrf')}
# A7 SDRF timing: file dates equal project dates in pride_files, so release era is the test
era = pd.cut(C.release_year, [2012, 2018, 2022], labels=['2013-2018', '2019-2022']).astype(str)
res['A7_sdrf_by_era'] = {f'{e}|sdrf={s}': {'n': int(len(g)), 'reused_pct': round(100 * g.D1.mean(), 1)}
                         for (e, s), g in C.groupby([era, C.sdrf])}
(OUT / 'Results.json').write_text(json.dumps(res, indent=2, default=str))
print(json.dumps(res, indent=1, default=str))
