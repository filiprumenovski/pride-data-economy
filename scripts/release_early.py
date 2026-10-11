"""Export shareable early-analysis inputs; recompute selected results fully offline.

Private handoff: python release_early.py export --source-root /path/to/handoff --data-dir DIR
Public release:  python release_early.py validate --data-dir DIR --output-dir RESULTS

Validation reads only the exported Parquets/CSV references, never private databases or
network resources. It covers breadth point estimates, scale, instrument credits,
software-country rates and the headline run-count bootstrap. It is not a rerun of
upstream geography inference, annotation detection or every sensitivity analysis.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd

SEED = 20261001


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def export(source, dest):
    sys.path.insert(0, str(source / 'work/stack-and-reuse'))
    from common import features
    from gallery4_breadth import long_frame
    out = source / 'outputs/Stack-and-reuse'
    d = features()
    assert len(d) == 36249 and d.accession.is_unique
    keys = d.pi_group.fillna('acc:' + d.accession)
    ordered = sorted(keys.unique())
    cross = pd.DataFrame({'original_group': ordered,
                          'contributor_group': [f'G{i:06d}' for i in range(1, len(ordered) + 1)]})
    private = out / 'Release-preparation/private'; private.mkdir(parents=True, exist_ok=True)
    cross.to_csv(private / 'Contributor-group-crosswalk.csv', index=False)
    mapping = cross.set_index('original_group').contributor_group
    groups = pd.Series(keys.map(mapping).to_numpy(), index=d.accession)
    known = pd.Series(d.pi_group.notna().to_numpy(), index=d.accession)
    dest.mkdir(parents=True, exist_ok=True)

    def public_frame(x):
        x = x.copy()
        if 'pi_group' in x:
            x['contributor_group_known'] = x.accession.map(known)
            x['pi_group'] = x.accession.map(groups)
            assert x.pi_group.notna().all()
        return x

    # Explicitly remove identity-bearing anchors, and retain no free-text protocol.
    safe = d.drop(columns=['anchor_institution'])
    safe = public_frame(safe)
    safe.to_parquet(dest / 'Deposit-features.parquet', index=False)
    study = pd.read_parquet(source / 'data/Study-cohort.parquet')
    study = study.drop(columns=['anchor_institution', 'publication_title', 'title',
                                'description', 'keywords', 'quant', 'experiments', 'software', 'processing'])
    public_frame(study).to_parquet(dest / 'Study-cohort.parquet', index=False)
    public_frame(pd.read_parquet(out / 'Layers/Dataset-layers.parquet')).to_parquet(dest / 'Dataset-layers.parquet', index=False)
    scale = pd.read_csv(out / 'Arsenal/H1-H3-scale-instruments/H1-dataset-scale.csv.gz')
    # Preserve feature input order for first-occurrence factorization in the bootstrap.
    scale = d[['accession']].merge(scale, on='accession', validate='one_to_one')
    public_frame(scale).to_parquet(dest / 'Dataset-scale.parquet', index=False)
    long = long_frame()
    long['grp'] = long.accession.map(groups)
    long.loc[long.axis == 'Contributor groups', 'cat'] = long.loc[long.axis == 'Contributor groups', 'grp']
    long = long.sort_values(['axis', 'accession', 'cat', 'grp']).reset_index(drop=True)
    long.to_parquet(dest / 'Breadth-long.parquet', index=False)
    for axis, name in [('Instrument models', 'Model-credits'), ('Species', 'Species-credits')]:
        long.loc[long.axis == axis].drop(columns='axis').to_parquet(dest / f'{name}.parquet', index=False)
    refs = {
        'breadth-annual.csv': 'Gallery/data/42-breadth.csv',
        'breadth-pooled.csv': 'Gallery/data/42-breadth-pooled.csv',
        'scale-annual.csv': 'Arsenal/H1-H3-scale-instruments/H1-scale-by-year.csv',
        'run-contrasts.csv': 'Attacks/runs/Contrasts-2013-16-vs-2023-25.csv',
        'models-annual.csv': 'Arsenal/H1-H3-scale-instruments/H3-model-by-year.csv',
        'software-country.csv': 'Software-geography/Tool-country-lift.csv',
    }
    refdir = dest / 'frozen-reference'; refdir.mkdir(exist_ok=True)
    for name, rel in refs.items():
        shutil.copyfile(out / rel, refdir / name)
        if name.startswith('breadth-'):
            ref = pd.read_csv(refdir / name)
            interval_columns = [c for c in ref if c in ['lo', 'hi', 'q1_lo', 'q1_hi', 'q2_lo', 'q2_hi']]
            ref.loc[ref.axis == 'Contributor groups', interval_columns] = np.nan
            ref.to_csv(refdir / name, index=False)
    records = []
    for p in sorted(dest.glob('*.parquet')):
        x = pd.read_parquet(p)
        assert not any(c in x for c in ['anchor_institution', 'email', 'description', 'processing', 'title'])
        records.append({'path': p.name, 'rows': len(x), 'bytes': p.stat().st_size,
                        'sha256': hashlib.sha256(p.read_bytes()).hexdigest(),
                        'schema': {c: str(t) for c, t in x.dtypes.items()}})
    write_json(dest / 'Data-manifest.json', {'snapshot': '2026-04-07', 'literature_retrieval': '2026-09-30',
        'files': records, 'reference_correction': 'Unsupported historical contributor-diversity uncertainty intervals are blanked in public reference copies; point estimates are unchanged.', 'grouping': 'Public G identifiers are a one-to-one relabeling of original groups; missing identities use accession-specific groups. contributor_group_known preserves original missingness.',
        'limits': 'Frozen derived metadata, not experimental mass-spectrometry files. Upstream annotations/geography cannot be reconstructed from these redacted inputs alone.'})
    print(json.dumps({'parquet_files': len(records), 'groups': len(cross), 'rows': len(d)}))


def hill(x, q=2):
    p = np.asarray(x, float); p = p / p.sum()
    return float(1 / np.square(p).sum() if q == 2 else np.exp(-np.sum(p[p > 0] * np.log(p[p > 0]))))


def weighted_median(x, w):
    k = w > 0; x, w = x[k], w[k]
    order = np.argsort(x, kind='stable'); x, w = x[order], w[order]
    cumulative = np.cumsum(w); midpoint = cumulative[-1] / 2
    i = np.searchsorted(cumulative, midpoint)
    return (x[i] + x[i + 1]) / 2 if np.isclose(cumulative[i], midpoint) and i + 1 < len(x) else x[i]


def run_bootstrap(scale, reps):
    z = scale[(scale.runs_lb > 0) & (scale.bytes_raw_layer > 0) &
              (scale.release_year.between(2013, 2016) | scale.release_year.between(2023, 2025))]
    late = z.release_year.to_numpy() >= 2023; v = z.runs_lb.to_numpy(float)
    group = pd.factorize(z.pi_group)[0]; ng = group.max() + 1
    rng = np.random.default_rng(SEED)
    vals = []
    for _ in range(reps):
        wt = np.bincount(rng.integers(0, ng, ng), minlength=ng)[group].astype(float)
        vals.append(weighted_median(v[late], wt[late]) / weighted_median(v[~late], wt[~late]))
    return {'n_early': int((~late).sum()), 'n_late': int(late.sum()),
            'median_early': float(np.median(v[~late])), 'median_late': float(np.median(v[late])),
            'median_ratio': float(np.median(v[late]) / np.median(v[~late])),
            'median_ratio_lo': float(np.quantile(vals, .025)), 'median_ratio_hi': float(np.quantile(vals, .975))}


def breadth_resampling(long, reps=300):
    """New deterministic resampling, separate from frozen manuscript Monte Carlo values."""
    from scipy.sparse import coo_matrix
    rng = np.random.default_rng(SEED)
    rows = []
    for axis in ['Instrument models', 'Journals', 'Species']:
        a = long[long.axis == axis].sort_values(['accession', 'cat', 'grp'])
        cats = pd.Index(sorted(a.cat.unique()))
        samples = {}; richness = {}; sizes = {}
        for label, window in [('early', (2013, 2015)), ('late', (2023, 2025))]:
            y = a[a.release_year.between(*window)]
            groups = pd.factorize(y.grp, sort=True)[0]
            matrix = coo_matrix((y.w, (groups, cats.get_indexer(y.cat))),
                                shape=(groups.max() + 1, len(cats))).tocsr()
            ng = matrix.shape[0]
            values = []
            for _ in range(reps):
                weights = rng.multinomial(ng, np.full(ng, 1 / ng))
                values.append(hill(np.asarray(weights @ matrix).ravel()))
            samples[label] = np.asarray(values)
            sizes[label] = y.accession.nunique()
        n = min(sizes.values())
        for label, window in [('early', (2013, 2015)), ('late', (2023, 2025))]:
            y = a[a.release_year.between(*window)]
            accessions = y.accession.unique()
            richness[label] = float(np.mean([y[y.accession.isin(rng.choice(accessions, n, replace=False))].cat.nunique() for _ in range(50)]))
        ratios = samples['late'] / samples['early']
        rows.append({'axis': axis, 'q2_lo': float(np.quantile(ratios, .025)),
                     'q2_hi': float(np.quantile(ratios, .975)), 'replicates': reps,
                     'q0_early_rarefied': richness['early'], 'q0_late_rarefied': richness['late'],
                     'rarefaction_replicates': 50, 'rarefied_to': n, 'seed': SEED})
    return pd.DataFrame(rows)


def validate(dest, result, reps):
    result.mkdir(parents=True, exist_ok=True)
    d = pd.read_parquet(dest / 'Deposit-features.parquet')
    scale = pd.read_parquet(dest / 'Dataset-scale.parquet')
    long = pd.read_parquet(dest / 'Breadth-long.parquet')
    checks = []

    def check(name, actual, expected, tolerance=1e-10):
        actual, expected = np.asarray(actual, float), np.asarray(expected, float)
        passed = bool(np.allclose(actual, expected, rtol=tolerance, atol=tolerance, equal_nan=True))
        checks.append({'name': name, 'passed': passed, 'values_compared': int(actual.size),
                       'max_absolute_error': float(np.nanmax(np.abs(actual - expected)))})
        assert passed, name

    check('cohort size and resolved geography', [len(d), d.country.notna().sum()], [36249, 28225])
    assert d.accession.is_unique and d.pi_group.str.fullmatch(r'G\d{6}').all()
    check('breadth credit sums within dataset and axis', long.groupby(['accession', 'axis']).w.sum(), 1)
    assert d.standards.equals(d.nongenerated_expanded)
    annual = []
    pooled = []
    for axis, a in long.groupby('axis'):
        for year, y in a[a.release_year >= 2016].groupby('release_year'):
            annual.append({'axis': axis, 'year': year, 'datasets': y.accession.nunique(), 'categories': y.cat.nunique(), 'eff': hill(y.groupby('cat').w.sum())})
        row = {'axis': axis}
        for label, window in [('early', (2013, 2015)), ('late', (2023, 2025))]:
            y = a[a.release_year.between(*window)]
            row['datasets_' + label] = y.accession.nunique()
            for q in (1, 2): row[f'q{q}_{label}'] = hill(y.groupby('cat').w.sum(), q)
        for q in (1, 2): row[f'q{q}_ratio'] = row[f'q{q}_late'] / row[f'q{q}_early']
        pooled.append(row)
    annual = pd.DataFrame(annual).sort_values(['axis', 'year']); pooled = pd.DataFrame(pooled).sort_values('axis')
    ref = pd.read_csv(dest / 'frozen-reference/breadth-annual.csv').sort_values(['axis', 'year'])
    for c in ['datasets', 'categories', 'eff']: check('annual breadth ' + c, annual[c], ref[c])
    ref = pd.read_csv(dest / 'frozen-reference/breadth-pooled.csv').sort_values('axis')
    for c in pooled.columns.drop('axis'): check('pooled breadth ' + c, pooled[c], ref[c])
    annual.to_csv(result / 'Breadth-annual-points.csv', index=False); pooled.to_csv(result / 'Breadth-pooled-points.csv', index=False)
    s = scale[scale.runs_lb >= 1]
    annual_scale = s.groupby('release_year').agg(datasets_with_raw=('accession', 'size'), median_runs_lb=('runs_lb', 'median'),
        share_ge100_runs=('runs_lb', lambda v: (v >= 100).mean())).reset_index()
    ref = pd.read_csv(dest / 'frozen-reference/scale-annual.csv').sort_values('release_year')
    for c in ['datasets_with_raw', 'median_runs_lb', 'share_ge100_runs']: check('annual scale ' + c, annual_scale[c], ref[c])
    annual_scale.to_csv(result / 'Scale-annual-points.csv', index=False)
    boot = run_bootstrap(scale, reps)
    ref = pd.read_csv(dest / 'frozen-reference/run-contrasts.csv').set_index('label').loc['all datasets with runs (claim population)']
    for c in ['n_early', 'n_late', 'median_early', 'median_late', 'median_ratio']: check('run contrast ' + c, boot[c], ref[c])
    if reps == 1000:
        for c in ['median_ratio_lo', 'median_ratio_hi']: check('run bootstrap ' + c, boot[c], ref[c])
    write_json(result / 'Run-count-bootstrap.json', {'replicates': reps, 'seed': SEED, **boot})
    models = long[long.axis == 'Instrument models'].groupby(['release_year', 'cat']).agg(datasets=('accession', 'nunique'), credit=('w', 'sum')).reset_index().rename(columns={'cat': 'model'})
    ref = pd.read_csv(dest / 'frozen-reference/models-annual.csv')
    joined = models.merge(ref, on=['release_year', 'model'], validate='one_to_one', suffixes=('', '_ref'))
    check('instrument model datasets', joined.datasets, joined.datasets_ref); check('instrument fractional credits', joined.credit, joined.credit_ref)
    models.to_csv(result / 'Model-annual-credits.csv', index=False)
    # Original equal-group weights are within each country, including groups observed in multiple countries.
    geo = d[d.release_year.between(2019, 2025) & d.country.notna() & d.contributor_group_known].copy()
    keep = geo.groupby('country').pi_group.nunique(); geo = geo[geo.country.isin(keep[keep >= 60].index)].copy()
    geo['weight'] = 1 / geo.groupby(['country', 'pi_group']).accession.transform('size')
    ref = pd.read_csv(dest / 'frozen-reference/software-country.csv'); rates = []
    for tool, expected in ref.groupby('tool', sort=False):
        flag = geo['sw_' + tool].astype(float)
        global_share = (flag * geo.weight).sum() / geo.weight.sum()
        for r in expected.itertuples():
            z = geo[geo.country == r.country]
            share = (flag[z.index] * z.weight).groupby(z.pi_group).sum().mean()
            rates.append({'tool': tool, 'country': r.country, 'lift': share / global_share})
    rates = pd.DataFrame(rates); joined = rates.merge(ref, on=['tool', 'country'], suffixes=('', '_ref'), validate='one_to_one')
    check('software-country equal-group lifts', joined.lift, joined.lift_ref)
    rates.to_csv(result / 'Software-country-points.csv', index=False)
    breadth_resampling(long).to_csv(result / 'Breadth-new-deterministic-resampling.csv', index=False)
    write_json(result / 'Validation.json', {'status': 'passed', 'checks': checks,
        'scope': 'Offline point-estimate recomputation plus headline run-count cluster bootstrap from public derived inputs.',
        'not_recomputed': ['raw-file classification', 'geography inference', 'all sensitivity analyses', 'frozen breadth bootstrap intervals and rarefied richness'],
        'contributor_diversity': 'Descriptive point estimates only; invalid historical cluster intervals are not recomputed or endorsed.'})
    print(json.dumps({'status': 'passed', 'checks': len(checks), 'bootstrap_replicates': reps}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['export', 'validate'])
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--source-root', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--reps', type=int, default=1000)
    args = parser.parse_args()
    if args.command == 'export':
        if not args.source_root: parser.error('export requires --source-root')
        export(args.source_root.resolve(), args.data_dir)
    else:
        if not args.output_dir: parser.error('validate requires --output-dir')
        if args.reps < 1: parser.error('--reps must be positive')
        validate(args.data_dir, args.output_dir, args.reps)
