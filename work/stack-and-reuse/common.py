"""Shared definitions for the stack-and-reuse experiments (new, October 2026)."""
from pathlib import Path
import sys
import numpy as np, pandas as pd, duckdb
from scipy.optimize import minimize
from scipy.special import expit

R = Path(__file__).resolve().parents[2]
O = R / 'outputs/Stack-and-reuse'
O.mkdir(parents=True, exist_ok=True)
SEED = 20261001

# Search/quantification software detected in the structured software field or the
# data-processing protocol. A match means "named"; no match means "not detected".
SOFTWARE = {
        'mascot': r'mascot', 'sequest': r'sequest',
    # "peaks" alone is ordinary prose ("peaks were picked"); require a product context.
    'peaks': r'peaks ?(studio|x\b|xpro|x\+|db|online|software|ab\b|au\b|\d)|bioinformatics solutions|"name":"peaks|\(peaks\)',
    'scaffold': r'scaffold',
    'progenesis': r'progenesis', 'msgf': r'ms-?gf', 'peptideshaker': r'peptideshaker|searchgui',
    'comet': r'\bcomet\b', 'xtandem': r'x!? ?tandem|xtandem', 'byonic': r'byonic',
    'proteinpilot': r'proteinpilot|protein pilot|paragon', 'plgs': r'plgs|proteinlynx',
    'skyline': r'skyline', 'msfragger': r'msfragger|fragpipe', 'pfind': r'pfind|plink',
    'openms': r'openms|openswath', 'percolator': r'percolator', 'metamorpheus': r'metamorpheus',
    'maxquant': r'maxquant', 'pd': r'proteome discoverer', 'diann': r'dia-?nn',
    'spectronaut': r'spectronaut',
}
# Plain software grouping (no new categories): the six analysis platforms most often named in
# 2023-2025 deposits. SEQUEST and Percolator are not separate entries because about 80% of
# deposits naming them also name Proteome Discoverer, inside which they usually run.
MAIN_TOOLS = {'Mascot': 'mascot', 'Proteome Discoverer': 'pd', 'MaxQuant': 'maxquant',
              'FragPipe': 'msfragger', 'Spectronaut': 'spectronaut', 'DIA-NN': 'diann'}
FAMILY_ORDER = list(MAIN_TOOLS) + ['Several of these', 'Other tools only', 'No tool named']


def add_families(d):
    """`family` = the one main platform a deposit names; 'Several of these' if it names more than
    one; 'Other tools only' if it names only other software; 'No tool named' otherwise."""
    F = d[['sw_' + v for v in MAIN_TOOLS.values()]].to_numpy()
    n = F.sum(axis=1); anyother = d[['sw_' + k for k in SOFTWARE]].any(axis=1).to_numpy()
    single = np.array(list(MAIN_TOOLS))[F.argmax(axis=1)]
    d['family'] = np.where(n == 1, single, np.where(n > 1, 'Several of these',
                                                    np.where(anyother, 'Other tools only', 'No tool named')))
    d['dia_native_tool'] = d.sw_diann | d.sw_spectronaut   # "DIA-NN or Spectronaut named"
    return d


def software_flags(c=None):
    """Return accession x sw_* boolean frame from the tracked cohort text fields."""
    c = c or duckdb.connect()
    s = c.execute(f"select accession, lower(coalesce(software,'')) || ' ' || lower(coalesce(processing,'')) t "
                  f"from read_parquet('{R / 'data/Study-cohort.parquet'}')").df()
    for k, p in SOFTWARE.items():
        s['sw_' + k] = s.t.str.contains(p, regex=True)
    return s.drop(columns='t')


def features():
    d = pd.read_parquet(O / 'Deposit-features.parquet')
    return add_families(d.merge(software_flags(), on='accession', validate='one_to_one'))


# Fixed-penalty logistic regression identical to the frozen geography models.
ALPHA = .001


def fit(X, y):
    def objective(b):
        z = X @ b; loss = np.mean(np.logaddexp(0, z) - y * z) + ALPHA * np.dot(b[1:], b[1:]) / 2
        grad = X.T @ (expit(z) - y) / len(y); grad[1:] += ALPHA * b[1:]
        return loss, grad
    r = minimize(objective, np.zeros(X.shape[1]), jac=True, method='L-BFGS-B',
                 options={'maxiter': 1000, 'gtol': 1e-7, 'ftol': 1e-10})
    assert r.success, r.message
    return r.x


def matrices(tr, te, cols, rare=('instruments', 'primary_subfield', 'primary_topic')):
    X = [np.ones((len(tr), 1))]; T = [np.ones((len(te), 1))]
    for col in cols:
        counts = tr[col].value_counts(); levels = set(counts[counts >= (25 if col in rare else 1)].index)
        x = tr[col].where(tr[col].isin(levels), 'Other'); y = te[col].where(te[col].isin(levels), 'Other')
        for lev in sorted(set(x) | {'Other'})[1:]:
            X.append((x == lev).to_numpy(float)[:, None]); T.append((y == lev).to_numpy(float)[:, None])
    return np.hstack(X), np.hstack(T)


def cluster_boot(values, groups, reps=2000, seed=SEED):
    """Ratio-of-sums percentile interval resampling whole clusters."""
    z = pd.DataFrame({'g': groups, 'v': values, 'n': 1}).groupby('g')[['v', 'n']].sum().to_numpy()
    rng = np.random.default_rng(seed)
    b = [(w := z[rng.integers(0, len(z), len(z))])[:, 0].sum() / w[:, 1].sum() for _ in range(reps)]
    return np.quantile(b, .025), np.quantile(b, .975)
