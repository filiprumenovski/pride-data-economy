"""Run original analyses offline against the reduced, checksummed source bundle.

Output must be outside the repository's frozen outputs. Reuse it across steps.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import sys

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HERE = ROOT / 'work/stack-and-reuse'


def prepare(source, output):
    frozen = ROOT / 'outputs/Stack-and-reuse'
    if output == frozen or output.is_relative_to(frozen):
        raise ValueError('Use a separate reconstruction directory')
    manifest = json.loads((ROOT / 'data/reconstruction/Source-manifest.json').read_text())
    for entry in manifest['files']:
        path = (source / entry['path']).resolve()
        assert path.is_relative_to(source)
        with path.open('rb') as f:
            assert hashlib.file_digest(f, 'sha256').hexdigest() == entry['sha256'], path
    output.mkdir(parents=True, exist_ok=True)
    feature = output / 'Deposit-features.parquet'
    if not feature.exists():
        f = pd.read_parquet(ROOT / 'data/early/Deposit-features.parquet')
        order = pd.read_csv(source / 'Group-order.csv').set_index('accession')
        f['pi_group'] = f.accession.map(order.sampling_group)
        assert f.pi_group.notna().equals(f.contributor_group_known)
        f.to_parquet(feature, index=False)
    for path in (source / 'inputs').rglob('*'):
        if path.is_file():
            dest = output / path.relative_to(source / 'inputs')
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--font-dir', type=Path, help='Optional directory containing Liberation Sans TTFs')
    parser.add_argument('script', help='Original analysis module, without .py')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    source, output = args.source.resolve(), args.output.resolve()
    manifest = prepare(source, output)
    os.environ.update(PRIDE_RECONSTRUCTION_SOURCE=str(source), PRIDE_RECONSTRUCTION_OUTPUT=str(output))
    if args.font_dir:
        os.environ['PRIDE_RECONSTRUCTION_FONTS'] = str(args.font_dir.resolve())
    if args.script.startswith('gallery'):
        from matplotlib import font_manager
        if args.font_dir:
            for path in args.font_dir.glob('LiberationSans-*.ttf'):
                font_manager.fontManager.addfont(str(path))
        font_manager.findfont('Liberation Sans', fallback_to_default=False)
    script = (HERE / (args.script + '.py')).resolve()
    assert script.parent == HERE and script.is_file()
    sys.path.insert(0, str(HERE))
    sys.path.insert(0, str(ROOT / 'scripts'))
    import common
    import query
    common.O = output
    common.features = lambda: pd.read_parquet(output / 'Deposit-features.parquet')
    common.software_flags = lambda c=None: common.features()[['accession'] + ['sw_' + k for k in common.SOFTWARE]]

    def connect():
        c = duckdb.connect()
        # Date casts must not inherit the reader's local timezone.
        c.execute("set threads=1; set memory_limit='3GB'; set TimeZone='UTC'")
        c.execute('create schema g; create schema lit')
        for entry in manifest['files']:
            if 'table' in entry:
                path = str(source / entry['path']).replace("'", "''")
                c.execute(f"create view {entry['table']} as select * from read_parquet('{path}')")
        study = str(ROOT / 'data/early/Study-cohort.parquet').replace("'", "''")
        c.execute(f"create view study as select s.*, t.* exclude(accession) from read_parquet('{study}') s left join study_text t using(accession)")
        return c

    query.connect = connect
    sys.argv = [str(script)] + args.arguments
    if args.script == 'gallery' and args.arguments:
        namespace = runpy.run_path(str(script), run_name='reconstruction')
        for name in args.arguments:
            assert name.startswith('fig_') and callable(namespace.get(name)), name
            namespace[name]()
    else:
        runpy.run_path(str(script), run_name='__main__')
    marker = output / 'Completed'
    marker.mkdir(exist_ok=True)
    (marker / (args.script + '-' + '-'.join(args.arguments) + '.json')).write_text(json.dumps({
        'script': args.script, 'arguments': args.arguments,
        'source_checksums_verified': len(manifest['files']), 'status': 'executed',
        'note': 'Execution is not a numerical comparison; see the comparison report.'}, indent=2) + '\n')


if __name__ == '__main__':
    main()
