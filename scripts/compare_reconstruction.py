"""Compare rebuilt tables with a reference snapshot; never change the reference."""
import argparse
import json
from pathlib import Path
import pandas as pd


def compare(reference, rebuilt):
    checks = []
    for path in sorted(rebuilt.rglob('*.csv')):
        relative = path.relative_to(rebuilt)
        old = reference / relative
        if not old.exists():
            continue
        a, b = pd.read_csv(old), pd.read_csv(path)
        row = {'path': str(relative), 'reference_rows': len(a), 'rebuilt_rows': len(b)}
        if set(a.columns) != set(b.columns) or len(a) != len(b):
            row.update(status='changed', reason='schema or row count', reference_columns=list(a), rebuilt_columns=list(b))
        else:
            b = b[a.columns]
            keys = [k for k in a if not pd.api.types.is_numeric_dtype(a[k])]
            # Include common numeric keys so row order alone cannot create a change.
            keys += [k for k in a if k not in keys and k in ['year', 'release_year', 'paper_year', 'submission_year', 'age_years', 'utc_hour']]
            if keys:
                a = a.sort_values(keys, kind='stable').reset_index(drop=True)
                b = b.sort_values(keys, kind='stable').reset_index(drop=True)
            changes = {}
            for col in a:
                try:
                    pd.testing.assert_series_equal(a[col], b[col], check_dtype=False,
                        check_exact=pd.api.types.is_integer_dtype(a[col]), rtol=1e-8, atol=1e-9)
                except AssertionError:
                    bad = ~(a[col].eq(b[col]) | (a[col].isna() & b[col].isna()))
                    changes[col] = {'different_cells': int(bad.sum()), 'examples': [
                        {'row': int(i), 'reference': None if pd.isna(a.at[i,col]) else str(a.at[i,col]),
                         'rebuilt': None if pd.isna(b.at[i,col]) else str(b.at[i,col])} for i in a.index[bad][:4]]}
            row.update(status='changed' if changes else 'match', changes=changes)
        checks.append(row)
    return checks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--rebuilt', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--require-match', action='store_true')
    a = p.parse_args(); checks = compare(a.reference, a.rebuilt)
    assert checks, 'No comparable tables found'
    report = {'tables': len(checks), 'matched': sum(r['status'] == 'match' for r in checks), 'checks': checks}
    a.report.parent.mkdir(parents=True, exist_ok=True)
    a.report.write_text(json.dumps(report, indent=2) + '\n')
    print(f"Compared {len(checks)} tables; {report['matched']} match.")
    if a.require_match:
        assert report['matched'] == len(checks), 'Numerical differences: see comparison report'


if __name__ == '__main__':
    main()
