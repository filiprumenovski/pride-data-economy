"""Recompute late-result summaries from the public, identity-free derived inputs."""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd

SEED = 20261001


def validate(data, out):
    out.mkdir(parents=True, exist_ok=False)
    checks = {}

    def compare(name, actual):
        expected = pd.read_csv(data / 'expected' / name)
        actual = actual[expected.columns]
        actual.to_csv(out / name, index=False)
        pd.testing.assert_frame_equal(actual.reset_index(drop=True), expected, check_dtype=False, rtol=1e-10, atol=1e-10)
        checks[name] = {'matched': True, 'rows': len(actual)}

    f = pd.read_csv(data / 'Dataset-features.csv', index_col='accession', low_memory=False)
    d = pd.read_csv(data / 'Dataset-archive-flags.csv')
    expected = pd.read_csv(data / 'expected/Before-after-by-year.csv')
    columns = expected.columns.drop(['release_year', 'datasets']).tolist()
    actual = d.groupby('release_year')[columns].mean().mul(100).round(1)
    actual['datasets'] = d.groupby('release_year').size()
    compare('Before-after-by-year.csv', actual.reset_index())
    h = pd.read_csv(data / 'Dataset-header-flags.csv')
    merged = d.merge(h[['accession', 'recognisable_with_headers', 'adhoc_after_headers']], on='accession', validate='one_to_one')
    cols = pd.read_csv(data / 'expected/Recognisable-by-year-three-ways.csv').columns.drop('release_year').tolist()
    compare('Recognisable-by-year-three-ways.csv', merged.groupby('release_year')[cols].mean().mul(100).round(1).reset_index())
    d = pd.read_csv(data / 'Dataset-tool-outcomes.csv')
    s = d[(d.tools_named == 1) & (d.release_year >= 2019)].copy()
    s = s[s.since.isna() | (s.release_year > s.since)]
    s['cluster'] = s.accession.map(f.bootstrap_group_order)
    per = s.groupby(['export', 'tool']).agg(datasets=('accession', 'size'), standard_pct=('standard', lambda x: 100 * x.mean()), native_visible_pct=('native_visible', lambda x: 100 * x.mean())).round(1)
    compare('Standard-rate-by-tool-single-tool-datasets.csv', per.reset_index())
    cap = []
    for cls, g in s.groupby('export'):
        z = pd.DataFrame({'g': g.cluster.to_numpy(), 'v': g.standard.to_numpy(float), 'n': 1}).groupby('g')[['v', 'n']].sum().to_numpy()
        rng = np.random.default_rng(SEED)
        bs = [(w := z[rng.integers(0, len(z), len(z))])[:, 0].sum() / w[:, 1].sum() for _ in range(1000)]
        lo, hi = np.quantile(bs, [.025, .975])
        cap.append({'export_class': cls, 'datasets': len(g), 'groups': g.cluster.nunique(), 'tools': ', '.join(sorted(g.tool.unique())), 'standard_pct': 100 * g.standard.mean(), 'ci_low': 100 * lo, 'ci_high': 100 * hi, 'native_visible_pct': 100 * g.native_visible.mean(), 'archive_only_pct': 100 * (g.outcome == 'archive').mean()})
    compare('Standard-rate-by-export-capability.csv', pd.DataFrame(cap))
    d = pd.read_csv(data / 'Dataset-documentation.csv')
    g = d[d.group_known].sort_values(['group_id', 'submission_date', 'accession']).copy()
    records = []
    for k in ['doc_readme', 'doc_sdrf', 'doc_design_author', 'doc_params']:
        g['prev'] = g.groupby('group_id')[k].shift()
        x = g[g.prev.notna()]
        base = x.groupby('release_year')[k].mean()
        for pv in [True, False]:
            s = x[x.prev == pv]
            records.append({'documentation': k, 'previous_dataset_had_it': pv, 'pairs': len(s), 'groups': s.group_id.nunique(), 'pct_now': round(100 * s[k].mean(), 1), 'pct_expected_same_years': round(100 * s.release_year.map(base).mean(), 1)})
    compare('H8-documentation-persistence-within-groups.csv', pd.DataFrame(records))
    p = pd.read_csv(data / 'Linked-papers.csv', index_col='work_id')
    p = p[p.publication_year.between(2013, 2025)].copy()
    p['jn'] = p.journal.str.lower().str.strip()
    p = p[p.jn.notna()].copy()
    p['prot'] = p.jn.str.contains('proteom').astype(float)
    p.publication_year = p.publication_year.astype(int)
    years = np.arange(2013, 2026)
    gc, gu = pd.factorize(p.group_id)
    jc, ju = pd.factorize(p.jn)
    yi = p.publication_year.to_numpy() - 2013

    def top10(weights):
        ww = weights[gc]
        result = []
        for k in range(len(years)):
            mask = yi == k
            counts = np.bincount(jc[mask], weights=ww[mask], minlength=len(ju))
            result.append(np.sort(counts)[::-1][:10].sum() / counts.sum())
        return result

    rng = np.random.default_rng(SEED)
    estimate = top10(np.ones(len(gu)))
    bs = np.array([top10(np.bincount(rng.integers(0, len(gu), len(gu)), minlength=len(gu)).astype(float)) for _ in range(1000)])
    lo, hi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    num = np.zeros((len(gu), len(years))); den = np.zeros_like(num)
    np.add.at(num, (gc, yi), p.prot.to_numpy()); np.add.at(den, (gc, yi), 1)
    rng = np.random.default_rng(SEED + 5)
    bs = []
    for _ in range(1000):
        w = np.bincount(rng.integers(0, len(gu), len(gu)), minlength=len(gu))
        bs.append((w @ num) / (w @ den))
    plo, phi = np.nanpercentile(bs, [2.5, 97.5], axis=0)
    result = pd.DataFrame({'papers_with_journal': p.groupby('publication_year').size(), 'top10_share': estimate, 'top10_lo': lo, 'top10_hi': hi, 'proteomics_titled_share': num.sum(0) / den.sum(0), 'proteomics_lo': plo, 'proteomics_hi': phi, 'distinct_journals': p.groupby('publication_year').jn.nunique(), 'effective_journals': p.groupby('publication_year').jn.apply(lambda x: 1 / (x.value_counts(normalize=True) ** 2).sum())}, index=years)
    result.index.name = 'paper_year'
    compare('37-journals-concentration.csv', result.round(4).reset_index())
    report = {'checks': checks, 'scope': 'Recomputes downstream summaries from frozen derived classifications; does not rerun classification of raw archive listings or file headers.'}
    (out / 'Validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    validate(args.data, args.output)
