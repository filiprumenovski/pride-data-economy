"""Export safe frozen late-results inputs and rerun mention analyses offline.

export reads the private workspace; validate reads only --data and the two original
analysis scripts beside this file. Outputs never overwrite frozen study results.
No source text, file names, contact hashes, emails or author identities are exported.
"""
from pathlib import Path
import argparse
import contextlib
import hashlib
import json
import re
import sys
import time
import warnings

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / 'outputs/Stack-and-reuse'
PREP = BASE / 'Release-preparation'
SEED = 20261001


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def export(out):
    from common import features
    sys.path.insert(0, str(ROOT / 'scripts'))
    from query import connect
    out.mkdir(parents=True, exist_ok=False)
    f = features()
    keys = f.pi_group.fillna('acc:' + f.accession)
    mapping = {key: f'G{i:06d}' for i, key in enumerate(sorted(keys.unique()), 1)}
    crosswalk = PREP / 'private/Contributor-group-crosswalk.csv'
    if crosswalk.exists():
        saved = pd.read_csv(crosswalk).set_index('original_group').contributor_group.to_dict()
        assert mapping == saved, 'Global contributor mapping differs'
    f['group_id'] = keys.map(mapping)
    f['group_known'] = f.pi_group.notna()
    legacy_keys = f.pi_group.fillna('solo:' + f.accession)
    legacy_order = {key: i for i, key in enumerate(sorted(legacy_keys.unique()))}
    f['bootstrap_group_order'] = legacy_keys.map(legacy_order)
    g = f.set_index('accession').group_id
    # Keep known/missing distinction and original row order; sorted group order is preserved.
    fields = ['accession', 'release_year', 'country', 'dia_expanded', 'vendor_single',
              'standards', 'family', 'human', 'group_id', 'group_known', 'bootstrap_group_order', 'primary_field',
              'primary_citations', 'listed_bytes', 'submission_date', 'release_date', 'work_id']
    d = f[fields].set_index('accession')
    layers = pd.read_parquet(BASE / 'Layers/Dataset-layers.parquet').set_index('accession')
    habits = pd.read_csv(BASE / 'Habits/Dataset-form.csv.gz').set_index('accession')
    d = d.join(layers[['has_sample_metadata']]).join(habits[['recognisable']])
    con = connect()
    # Match the original A2 regex exactly, but publish only the Boolean result.
    source = (HERE / 'attack_reuse_mentions.py').read_text()
    regex_code = source[source.index('pat = re.compile('):source.index('\ntxt = ')]
    ns = {'re': re}; exec(regex_code, ns)
    projects = con.execute('select project_accession accession, project_title t, project_description d from pride_projects').df().set_index('accession')
    text = (projects.t.fillna('') + ' ' + projects.d.fillna('')).reindex(d.index).fillna('')
    d['method_like'] = text.str.contains(ns['pat'])
    d.to_csv(out / 'Dataset-features.csv', index=True)
    mentions = pd.read_csv(BASE / 'Reuse-mentions/Mention-pairs.csv.gz')
    mentions[['accession', 'pmid', 'y']].to_csv(out / 'Nonprimary-mention-pairs.csv', index=False)
    primary = con.execute('select project_accession accession, try_cast(publication_pubmed_id as bigint) pmid from pride_publications').df().dropna()
    primary = primary[primary.accession.isin(d.index)]
    primary.pmid = primary.pmid.astype('int64')
    primary.to_csv(out / 'Primary-publication-links.csv', index=False)
    expected = pd.read_csv(BASE / 'Reuse-mentions/Dataset-reuse.csv.gz')
    assert expected.accession.tolist() == d[d.release_year.between(2013, 2022)].index.tolist()
    expected[['accession', 'D1', 'D2', 'D3', 'D1_all', 'D2_all', 'D3_all']].to_csv(out / 'Expected-outcomes.csv', index=False)
    expected_pairs = mentions.drop(columns=['grp'])
    expected_pairs.to_csv(out / 'Expected-pair-derivations.csv', index=False)
    tables = {
        'Dataset-tool-outcomes.csv': ('Tool-outputs/Dataset-tool-outcomes.csv.gz',
            ['accession', 'tool', 'export', 'since', 'native_visible', 'standard', 'tool_archive', 'outcome', 'tools_named', 'release_year', 'dia_expanded']),
        'Dataset-documentation.csv': ('Arsenal/H5-H6-H8-strategies-software-docs/Dataset-strategy-software-docs.csv.gz',
            ['accession', 'release_year', 'doc_readme', 'doc_sdrf', 'doc_design_author', 'doc_params']),
        'Dataset-clocks.csv': ('Arsenal/H11-clocks/Dataset-clocks.csv.gz',
            ['accession', 'submission_date', 'release_date', 'work_id', 'release_year', 'paper_date', 'acq_date', 'db_date', 'software', 'sw_date',
             'acq_to_sub', 'acq_to_rel', 'db_to_sub', 'db_to_rel', 'sw_to_sub', 'sw_to_rel', 'sub_to_paper', 'paper_to_rel', 'sub_to_rel', 'db_minus_acq', 'sw_minus_acq']),
        'Dataset-archive-flags.csv': ('Archive-contents/Dataset-before-after.csv.gz', None),
        'Dataset-header-flags.csv': ('Header-sniff/Dataset-headers.csv.gz', None),
    }
    for name, (rel, cols) in tables.items():
        x = pd.read_csv(BASE / rel, usecols=cols)
        if name in ('Dataset-tool-outcomes.csv', 'Dataset-documentation.csv', 'Dataset-clocks.csv'):
            x['group_id'] = x.accession.map(g)
            x['group_known'] = x.accession.map(f.set_index('accession').group_known)
        if name == 'Dataset-documentation.csv':
            x['submission_date'] = x.accession.map(d.submission_date)
        x.to_csv(out / name, index=False)
    # Figure 5b-c needs paper grain and the paper's modal contributor group, not author data.
    works = con.execute('''select w.work_id, w.publication_year, e.journal
        from lit.openalex_works w left join lit.europepmc_works e on w.pmid=e.pmid''').df()
    assert works.work_id.is_unique
    dep = f[f.work_id.notna()].copy()
    oldgrp = dep.pi_group.fillna('acc:' + dep.accession)
    dep['grp'] = oldgrp
    paper = dep.groupby('work_id').agg(n_deposits=('accession', 'size'), grp=('grp', lambda s: s.mode().iloc[0]))
    paper['group_id'] = paper.pop('grp').map(mapping)
    paper = paper.join(works.set_index('work_id'))
    paper.to_csv(out / 'Linked-papers.csv')
    dep[['accession', 'work_id']].to_csv(out / 'Dataset-paper-links.csv', index=False)
    # Frozen expected aggregates are test fixtures, not inputs to the recomputation.
    reference = out / 'expected'; reference.mkdir()
    for rel in ['Reuse-mentions/Reuse-rates.csv', 'Reuse-mentions/Reused-pool-composition.csv',
                'Reuse-mentions/Summary.json', 'Attacks/reuse-mentions/Results.json',
                'Tool-outputs/Standard-rate-by-tool-single-tool-datasets.csv', 'Tool-outputs/Standard-rate-by-export-capability.csv',
                'Arsenal/H5-H6-H8-strategies-software-docs/H8-documentation-persistence-within-groups.csv',
                'Archive-contents/Before-after-by-year.csv', 'Header-sniff/Recognisable-by-year-three-ways.csv',
                'Gallery/data/37-journals-concentration.csv', 'Gallery/data/37-journals-top15.csv']:
        (reference / Path(rel).name).write_bytes((BASE / rel).read_bytes())
    files = []
    for p in sorted(out.rglob('*')):
        if p.is_file():
            content = p.read_text()
            assert not re.search(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}', content), p
            item = {'path': str(p.relative_to(out)), 'sha256': digest(p), 'bytes': p.stat().st_size}
            if p.suffix == '.csv':
                x = pd.read_csv(p); item.update(rows=len(x), columns=x.columns.tolist())
                assert not any(c in x for c in ['pi_group', 'grp', 'title', 'description', 'email', 'file_name', 'processing', 'software_protocol'])
            files.append(item)
    manifest = {'scope': 'Derived numerical/categorical inputs; extraction from raw files is not reproduced by this package.',
        'group_mapping': 'Full feature cohort pi_group.fillna(acc: + accession), sorted, numbered G000001; missingness retained as group_known. Private crosswalk excluded.',
        'bootstrap_group_order': 'Integer ordering of original groups with solo:accession fallback; retains frozen RNG group sampling order without publishing identities.',
        'method_like': 'Computed privately from exact A2 regex on frozen project title+description; only Boolean released.',
        'nonprimary_pairs': 'Frozen non-primary accession-PMID pairs preserve the original any_value publication-year selection; titles excluded.',
        'source_hashes': {p.name: digest(p) for p in [HERE / 'reuse_mentions.py', HERE / 'attack_reuse_mentions.py', Path(__file__)]},
        'files': files}
    (out / 'Input-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'export': str(out), 'files': len(files), 'datasets': len(d), 'groups': len(mapping)}, indent=2))


def bootstrap_stability(draws):
    """Split-half endpoint diagnostic; does not replace frozen 200-draw intervals."""
    columns = ['standard', 'dia', 'sdrf', 'human', 'thermo', 'log_citations']
    bounds = np.percentile(draws, [2.5, 97.5], axis=0)
    halves = [np.percentile(x, [2.5, 97.5], axis=0) for x in np.array_split(draws, 2)]
    return {'replicates': len(draws), 'note': 'Split-half comparison describes Monte Carlo sensitivity only, not confidence-interval validation; 200 draws leave about five draws per percentile tail.',
        'terms': {key: {'full_interval': bounds[:, i].tolist(), 'first_half_interval': halves[0][:, i].tolist(),
                        'second_half_interval': halves[1][:, i].tolist(),
                        'maximum_endpoint_difference': float(np.max(np.abs(halves[0][:, i] - halves[1][:, i])))}
                  for i, key in enumerate(columns)}}


def replace_once(text, before, after):
    assert text.count(before) == 1, f'Original source changed: {before[:80]}'
    return text.replace(before, after, 1)


def validate(data, out):
    """Run original analysis bodies with only their I/O replaced; retain RNG call order."""
    import statsmodels.api as sm
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads((data / 'Input-manifest.json').read_text())
    for item in manifest['files']:
        assert digest(data / item['path']) == item['sha256'], item['path']
    for name in ['reuse_mentions.py', 'attack_reuse_mentions.py']:
        assert digest(HERE / name) == manifest['source_hashes'][name], f'Source changed: {name}'
    d = pd.read_csv(data / 'Dataset-features.csv', index_col='accession')
    d['pi_group'] = d.group_id.where(d.group_known)
    d['grp'] = d.bootstrap_group_order
    for target, source in [('standard', 'standards'), ('dia', 'dia_expanded'), ('deep', 'recognisable')]:
        d[target] = d[source].fillna(False).astype(bool)
    d['thermo'] = d.vendor_single.eq('Thermo Fisher'); d['bruker'] = d.vendor_single.eq('Bruker')
    d['sdrf'] = d.has_sample_metadata.fillna(False).astype(bool)
    d['dia_tool_only'] = d.family.isin(['DIA-NN', 'Spectronaut']); d['human'] = d.human.fillna(False).astype(bool)
    m = pd.read_csv(data / 'Nonprimary-mention-pairs.csv')
    m = m.join(d[['release_year', 'grp']], on='accession'); m['lag'] = m.y - m.release_year
    primary = pd.read_csv(data / 'Primary-publication-links.csv')
    own = primary.join(d[['grp']], on='accession').dropna()[['pmid', 'grp']].drop_duplicates()
    m['same_group_paper'] = m.set_index(['pmid', 'grp']).index.isin(own.set_index(['pmid', 'grp']).index)
    size = m.groupby('pmid').agg(k=('accession', 'nunique'), kg=('grp', 'nunique')); m = m.join(size, on='pmid')
    expected_pairs = pd.read_csv(data / 'Expected-pair-derivations.csv')
    pd.testing.assert_frame_equal(m.drop(columns='grp')[expected_pairs.columns], expected_pairs, check_dtype=False)
    cols = ['release_year', 'country', 'dia_expanded', 'vendor_single', 'standards', 'family', 'pi_group', 'human',
            'has_sample_metadata', 'recognisable', 'grp', 'standard', 'thermo', 'bruker', 'sdrf', 'deep', 'dia', 'dia_tool_only']
    stage = out / 'Reuse-mentions'; stage.mkdir()
    source = (HERE / 'reuse_mentions.py').read_text(); body = source[source.index('def links('):]
    body = replace_once(body, ",\n                                          title=('title', 'first')", '')
    body = replace_once(body, "M.drop(columns=['title']).to_csv", "M.to_csv")
    ns = {'np': np, 'pd': pd, 'json': json, 'd': d[cols], 'M': m, 'size': size,
          'OUT': stage, 'REPS': 1000, 'rng': np.random.default_rng(SEED)}
    start = time.monotonic()
    with (out / 'Run.log').open('w') as log, contextlib.redirect_stdout(log):
        exec(compile(body, 'reuse_mentions.py [offline I/O]', 'exec'), ns)
    outcomes = pd.read_csv(stage / 'Dataset-reuse.csv.gz').set_index('accession')
    expected = pd.read_csv(data / 'Expected-outcomes.csv').set_index('accession')
    pd.testing.assert_frame_equal(outcomes[expected.columns], expected)
    checks = {'pair_derivations_exact': True, 'outcomes_exact': True}
    for name in ['Reuse-rates.csv', 'Reused-pool-composition.csv']:
        pd.testing.assert_frame_equal(pd.read_csv(stage / name), pd.read_csv(data / 'expected' / name), check_dtype=False)
        checks[name] = 'exact'
    print('Mention-rate recomputation matched frozen tables.', flush=True)
    source = (HERE / 'attack_reuse_mentions.py').read_text()
    body = source[source.index('CONTRASTS = '):]
    a2_start = body.index('P = connect().execute('); a2_end = body.index('\nC2 = ', a2_start)
    body = body[:a2_start] + "C['method_like'] = safe_features.method_like.reindex(C.index)\n" + body[a2_end:]
    fit_start = body.index('def fit_or(idx):'); fit_end = body.index('\n\nest = ', fit_start)
    body = body[:fit_start] + '''def fit_or(idx):
    model = sm.GLM(y[idx], Xa[idx], family=sm.families.Binomial())
    begun = time.monotonic()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        fit = model.fit_regularized(alpha=1e-4, L1_wt=0)
    params = np.asarray(fit.params)
    gradient = -model.score(params) / len(idx) + 1e-4 * params
    diagnostics.append({'fit': len(diagnostics), 'rows': len(idx), 'events': int(y[idx].sum()),
        'finite_params': bool(np.isfinite(params).all()), 'score_inf_norm': float(np.max(np.abs(gradient))),
        'warnings': [str(w.message) for w in caught], 'seconds': time.monotonic() - begun,
        'converged_attribute': bool(fit.converged) if hasattr(fit, 'converged') else None})
    assert np.isfinite(params).all(), 'Nonfinite fitted coefficient'
    if len(diagnostics) % 25 == 0:
        progress(len(diagnostics))
    return np.exp(params[[X.columns.get_loc(k) for k in ['standard', 'dia', 'sdrf', 'human', 'thermo', 'log_citations']]])
''' + body[fit_end:]
    attack = out / 'Attacks/reuse-mentions'; attack.mkdir(parents=True)
    covars = d[['primary_field', 'primary_citations', 'listed_bytes']]
    C = pd.read_csv(stage / 'Dataset-reuse.csv.gz', index_col=0).join(covars)
    diagnostics = []
    # Progress writes to original stdout, outside the analysis log redirection.
    stdout = sys.stdout
    ns = {'np': np, 'pd': pd, 'json': json, 're': re, 'sm': sm, 'C': C, 'M': m,
          'safe_features': d, 'OUT': attack, 'rng': np.random.default_rng(SEED), 'time': time,
          'warnings': warnings, 'diagnostics': diagnostics,
          'progress': lambda n: print(f'Model fits completed: {n}/201', file=stdout, flush=True)}
    with (out / 'Run.log').open('a') as log, contextlib.redirect_stdout(log):
        exec(compile(body, 'attack_reuse_mentions.py [offline I/O + diagnostics]', 'exec'), ns)
    observed = json.loads((attack / 'Results.json').read_text())
    frozen = json.loads((data / 'expected/Results.json').read_text())
    checks['attack_results_exact'] = observed == frozen
    checks['model_design_shape'] = list(ns['X'].shape)
    checks['model_columns'] = ns['X'].columns.tolist()
    checks['fits'] = len(diagnostics)
    checks['finite_fits'] = sum(x['finite_params'] for x in diagnostics)
    checks['fit_warnings'] = sum(bool(x['warnings']) for x in diagnostics)
    checks['maximum_penalized_score_inf_norm'] = max(x['score_inf_norm'] for x in diagnostics)
    checks['seconds'] = time.monotonic() - start
    (out / 'Bootstrap-stability.json').write_text(json.dumps(bootstrap_stability(ns['bs']), indent=2) + '\n')
    (out / 'Model-diagnostics.json').write_text(json.dumps(diagnostics, indent=2) + '\n')
    (out / 'Validation.json').write_text(json.dumps(checks, indent=2) + '\n')
    pd.DataFrame(ns['bs'], columns=['standard', 'dia', 'sdrf', 'human', 'thermo', 'log_citations']).to_csv(out / 'Bootstrap-odds-ratios.csv', index=False)
    ns['X'].to_csv(out / 'Model-design.csv', index_label='accession')
    print(json.dumps(checks, indent=2))
    assert checks['attack_results_exact'], 'Frozen attack results differ; inspect diagnostic outputs'
    assert len(diagnostics) == 201


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['export', 'validate'])
    parser.add_argument('--data', type=Path, default=PREP / 'public-data/late')
    parser.add_argument('--output', type=Path, default=PREP / 'validation/late')
    args = parser.parse_args()
    if args.action == 'export':
        export(args.data)
    else:
        validate(args.data, args.output)
