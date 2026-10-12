"""Recompute Figure 5a / Table S15 from public dates, including 1,000-draw intervals."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def validate(data, output):
    output.mkdir(parents=True, exist_ok=False)
    d = pd.read_csv(data / 'Dataset-clocks.csv')
    f = pd.read_csv(data / 'Dataset-features.csv', low_memory=False).set_index('accession')
    d['grp'] = d.accession.map(f.bootstrap_group_order)
    assert d.grp.notna().all() and d.accession.is_unique
    expected = json.loads((data / 'expected/Clocks-summary.json').read_text())
    dates = ['submission_date', 'release_date', 'paper_date', 'acq_date', 'db_date', 'sw_date']
    for col in dates:
        d[col] = pd.to_datetime(d[col])
    pairs = {
        'acq_to_sub': ('submission_date', 'acq_date'), 'acq_to_rel': ('release_date', 'acq_date'),
        'db_to_sub': ('submission_date', 'db_date'), 'db_to_rel': ('release_date', 'db_date'),
        'sw_to_sub': ('submission_date', 'sw_date'), 'sw_to_rel': ('release_date', 'sw_date'),
        'sub_to_paper': ('paper_date', 'submission_date'), 'paper_to_rel': ('release_date', 'paper_date'),
        'sub_to_rel': ('release_date', 'submission_date'), 'db_minus_acq': ('db_date', 'acq_date'),
        'sw_minus_acq': ('sw_date', 'acq_date'),
    }
    frozen = d[list(pairs)].copy()
    for col, (a, b) in pairs.items():
        d[col] = (d[a] - d[b]).dt.days.astype(float)
    for col in ['acq_to_sub', 'db_to_sub', 'sw_to_sub']:
        d.loc[d[col] < -31, [col, col.replace('_sub', '_rel')]] = np.nan
    d.loc[d.sub_to_rel < 0, ['sub_to_rel', 'paper_to_rel']] = np.nan
    pd.testing.assert_frame_equal(d[list(pairs)], frozen, check_dtype=False)
    rng = np.random.default_rng(20261001)

    def median_ci(col, frame):
        x = frame[[col, 'grp']].dropna()
        groups = list(x.groupby('grp').indices.values())
        v = x[col].to_numpy()
        bs = [np.median(v[np.concatenate([groups[k] for k in rng.integers(0, len(groups), len(groups))])]) for _ in range(1000)]
        return {'datasets': len(x), 'groups': len(groups), 'median_days': float(np.median(v)),
                'median_years': round(float(np.median(v)) / 365.25, 2),
                'ci_years': [round(np.percentile(bs, 2.5) / 365.25, 2), round(np.percentile(bs, 97.5) / 365.25, 2)],
                'share_gt_1y': round(float((v > 365.25).mean()), 3), 'share_gt_2y': round(float((v > 730.5).mean()), 3)}

    clocks = {col: median_ci(col, d) for col in list(pairs)[:9]}
    by_tool = {tool: {col: median_ci(col, d[d.software == tool]) for col in ['sw_to_sub', 'sw_to_rel']}
               for tool in ['MaxQuant', 'DIA-NN', 'FragPipe', 'Skyline', 'Mascot']}
    assert clocks == expected['clocks'], 'Clock summary differs'
    assert by_tool == expected['software_by_tool'], 'Software-specific clock summary differs'
    result = {'status': 'passed', 'date_differences_checked': len(pairs), 'dataset_rows': len(d),
              'median_summaries_checked': 19, 'bootstrap_replicates_per_summary': 1000,
              'scope': 'Dates to summary estimates and intervals; extraction of dates from original filenames is upstream.'}
    (output / 'Clocks-summary.json').write_text(json.dumps({'clocks': clocks, 'software_by_tool': by_tool}, indent=2) + '\n')
    (output / 'Validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('data/late'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    validate(args.data, args.output)
