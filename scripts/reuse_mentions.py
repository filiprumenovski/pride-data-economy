"""Which PRIDE datasets do later papers actually reuse? (descriptive; local Europe PMC mentions only)

Source: europepmc_mentions (text-mined accessions, largely open-access full text; snapshot 2026-04-07).
A mention is a (paper, dataset) pair from a paper that is not the dataset's recorded primary reference.
Most such pairs are the dataset's own paper under another identifier (median lag 0 years; ledger G3), so
reuse is defined conservatively, with three definitions used as sensitivities:
  D1 cross-group   the paper names datasets from >= 2 contributor groups (pi_group; solo datasets are their own group)
  D2 later         mention >= 2 years after release, and the paper is not the primary paper of a dataset of the same group
  D3 both          D1 and D2
Papers naming > 100 datasets are excluded as mining artefacts (one record links 558; sensitivity keeps them).
Exposure: mentions in the release year to release + 3 years, datasets released 2013-2022 (complete window
before the 2025 publication year). Rates are compared within release year and species (human vs other)
by direct standardisation to the cohort's year x human mix, with group-bootstrap 95% intervals (pi_group).
Outputs: outputs/Stack-and-reuse/Reuse-mentions/.
"""
import json, sys
import numpy as np, pandas as pd
from common import O, R, SEED, features
sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Reuse-mentions'; OUT.mkdir(parents=True, exist_ok=True)
REPS = 1000
rng = np.random.default_rng(SEED)
c = connect()

F = features().set_index('accession')
L = pd.read_parquet(O / 'Layers/Dataset-layers.parquet').set_index('accession')
H = pd.read_csv(O / 'Habits/Dataset-form.csv.gz').set_index('accession')
d = F[['release_year', 'country', 'dia_expanded', 'vendor_single', 'standards', 'family', 'pi_group', 'human']] \
    .join(L[['has_sample_metadata']]).join(H[['recognisable']])
d['grp'] = d.pi_group.fillna('solo:' + d.index.to_series())
d['standard'] = d.standards.astype(bool)
d['thermo'] = d.vendor_single.eq('Thermo Fisher')
d['bruker'] = d.vendor_single.eq('Bruker')
d['sdrf'] = d.has_sample_metadata.fillna(False).astype(bool)
d['deep'] = d.recognisable.fillna(False).astype(bool)
d['dia'] = d.dia_expanded.astype(bool)
d['dia_tool_only'] = d.family.isin(['DIA-NN', 'Spectronaut'])
d['human'] = d.human.fillna(False).astype(bool)

M = c.execute("""select project_accession accession, mentioning_pmid pmid, any_value(mentioning_pub_year) y, any_value(mentioning_title) title
                 from europepmc_mentions where not is_primary_reference and mentioning_pmid > 0 group by 1, 2""").df()
M = M[M.accession.isin(d.index)].join(d[['release_year', 'grp']], on='accession')
M['lag'] = M.y - M.release_year
P = c.execute("select project_accession accession, try_cast(publication_pubmed_id as bigint) pmid from pride_publications").df().dropna()
own = P.join(d[['grp']], on='accession').dropna()[['pmid', 'grp']].drop_duplicates()
M['same_group_paper'] = M.set_index(['pmid', 'grp']).index.isin(own.set_index(['pmid', 'grp']).index)
size = M.groupby('pmid').agg(k=('accession', 'nunique'), kg=('grp', 'nunique'))
M = M.join(size, on='pmid')


def links(M, max_k=100):
    m = M[M.k <= max_k] if max_k else M
    m = m[m.lag.between(0, 3)]
    return {'D1': m[m.kg >= 2], 'D2': m[(m.lag >= 2) & ~m.same_group_paper], 'D3': m[(m.kg >= 2) & (m.lag >= 2) & ~m.same_group_paper]}


C = d[d.release_year.between(2013, 2022)].copy()
for name, m in links(M).items():
    C[name] = C.index.isin(m.accession)
for name, m in links(M, None).items():
    C[name + '_all'] = C.index.isin(m.accession)
C['stratum'] = C.release_year.astype(str) + '|' + C.human.astype(str)
C.drop(columns=['recognisable', 'has_sample_metadata']).to_csv(OUT / 'Dataset-reuse.csv.gz')
M.drop(columns=['title']).to_csv(OUT / 'Mention-pairs.csv.gz', index=False)
W = C.stratum.value_counts(normalize=True)


def boot(df, col):
    """Rate directly standardised to the cohort year x human mix; group bootstrap as multinomial group weights."""
    g = pd.factorize(df.grp)[0]; s = pd.Categorical(df.stratum, categories=W.index).codes
    N = np.zeros((g.max() + 1, len(W))); Y = np.zeros_like(N)
    np.add.at(N, (g, s), 1); np.add.at(Y, (g, s), df[col].to_numpy(float))
    w = rng.multinomial(len(N), np.full(len(N), 1 / len(N)), size=REPS)        # REPS x groups
    wt = W.to_numpy()

    def rate(n, y):
        ok = n > 0
        return (np.where(ok, y / np.where(ok, n, 1), 0) * wt * ok).sum(-1) / (wt * ok).sum(-1)
    est = rate(N.sum(0), Y.sum(0)); bs = rate(w @ N, w @ Y)
    return est, np.percentile(bs, 2.5), np.percentile(bs, 97.5)


FEATS = {'Standard files': 'standard', 'Thermo raw': 'thermo', 'Bruker': 'bruker', 'DIA': 'dia', 'DIA-NN/Spectronaut': 'dia_tool_only',
         'SDRF': 'sdrf', 'Results identifiable from files': 'deep', 'Human': 'human'}
rows = []
for outcome in ['D1', 'D2', 'D3', 'D1_all']:
    rows.append({'outcome': outcome, 'feature': 'All', 'level': 'all', 'n': len(C), 'raw_pct': 100 * C[outcome].mean(),
                 **dict(zip(['std_pct', 'lo', 'hi'], [100 * x for x in boot(C, outcome)]))})
    for lab, col in FEATS.items():
        for val in [True, False]:
            S = C[C[col] == val]
            if len(S) < 50:
                continue
            est = boot(S, outcome)
            rows.append({'outcome': outcome, 'feature': lab, 'level': str(val), 'n': len(S), 'raw_pct': 100 * S[outcome].mean(),
                         **dict(zip(['std_pct', 'lo', 'hi'], [100 * x for x in est]))})
    for by, k in [('family', 8), ('country', 12), ('vendor_single', 5)]:
        for val in C[by].value_counts().head(k).index:
            S = C[C[by] == val]
            rows.append({'outcome': outcome, 'feature': by, 'level': val, 'n': len(S), 'raw_pct': 100 * S[outcome].mean(),
                         **dict(zip(['std_pct', 'lo', 'hi'], [100 * x for x in boot(S, outcome)]))})
T = pd.DataFrame(rows).round(2)
T.to_csv(OUT / 'Reuse-rates.csv', index=False)

# composition of the reused pool vs the cohort, within the same release years
comp = []
for outcome in ['D1', 'D3']:
    Rp = C[C[outcome]]
    for lab, col in FEATS.items():
        comp.append({'outcome': outcome, 'feature': lab, 'cohort_pct': 100 * C[col].mean(), 'reused_pct': 100 * Rp[col].mean(),
                     'reused_n': len(Rp)})
pd.DataFrame(comp).round(1).to_csv(OUT / 'Reused-pool-composition.csv', index=False)

# the reusing papers: largest cross-group papers (titles are public metadata)
top = M[(M.kg >= 2)].groupby('pmid').agg(datasets=('accession', 'nunique'), groups=('grp', 'nunique'), year=('y', 'first'),
                                          title=('title', 'first')).sort_values('datasets', ascending=False)
top.head(60).to_csv(OUT / 'Largest-cross-group-papers.csv')
summary = {
    'cohort': 'PRIDE datasets released 2013-2022; mentions in release year to +3 years', 'datasets': int(len(C)), 'groups': int(C.grp.nunique()),
    'mention_pairs_secondary': int(len(M)), 'papers_secondary': int(M.pmid.nunique()),
    'lag_share_0_or_less': round(float((M.lag <= 0).mean()), 3),
    'cross_group_papers': int((size.kg >= 2).sum()), 'papers_over_100_datasets': int((size.k > 100).sum()),
    'reused_datasets': {o: int(C[o].sum()) for o in ['D1', 'D2', 'D3', 'D1_all']},
    'notes': ['Mentions are text-mined, largely from open-access full text; reuse cited only in supplementary tables is missed.',
              'A mention is not proof of reanalysis; D1-D3 approximate third-party use.',
              'Rates are standardised to the cohort release-year x human mix.'],
}
(OUT / 'Summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=1))
print(T[T.outcome.isin(['D1', 'D3'])].to_string())
print(pd.read_csv(OUT / 'Reused-pool-composition.csv').to_string())
print(top.head(25).to_string())
