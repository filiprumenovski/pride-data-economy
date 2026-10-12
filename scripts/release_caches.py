"""Reclassify the October 11–12 retrieval offline; compare with frozen paper inputs.

Download and extract Extraction-caches.zip from the matching GitHub release first.
Use a fresh output directory on a disk with several GB available for DuckDB spill.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import runpy
import shutil
import sys

import pandas as pd
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def compare_flags(frozen, current):
    a = pd.read_csv(frozen).set_index('accession').sort_index()
    b = pd.read_csv(current).set_index('accession').sort_index()
    pd.testing.assert_index_equal(a.index, b.index)
    assert a.index.is_unique
    changes = {}
    for col in a.columns:
        x, y = a[col], b[col]
        different = ~(x.eq(y) | (x.isna() & y.isna()))
        if different.any():
            changes[col] = {'datasets': int(different.sum()),
                            'accessions': a.index[different].tolist()}
    return changes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    cache, output = args.cache_dir.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    inputs = ROOT / 'data/extraction'
    manifest = json.loads((inputs / 'Cache-manifest.json').read_text())
    for entry in manifest['outputs'].values():
        path = (cache / entry['path']).resolve()
        assert path.is_relative_to(cache), path
        assert path.stat().st_size == entry['bytes'], path
        assert sha256(path) == entry['sha256'], path
    print('Three extraction-cache checksums verified.', flush=True)
    (output / 'Layers').mkdir()
    shutil.copy2(ROOT / 'data/early/Dataset-layers.parquet', output / 'Layers/Dataset-layers.parquet')
    for name in ['Archive-contents', 'Header-sniff']:
        (output / name).mkdir()
        target = inputs / (name + '-Targets.csv.gz')
        assert sha256(target) == manifest['target_sha256'][name]
        shutil.copy2(target, output / name / 'Targets.csv.gz')

    ac = output / 'Archive-contents'
    (ac / 'Members.parquet').symlink_to(cache / 'Archive-contents/Members.parquet')
    targets = pd.read_csv(ac / 'Targets.csv.gz')
    wanted = set(targets.loc[(targets.fmt == 'zip') & (targets.file_size_bytes > 0), 'url'])
    seen, members = set(), 0
    fields = ['accession', 'url', 'file_name', 'layer', 'cat', 'size', 'status', 'n_members']
    with gzip.open(cache / 'Archive-contents/Listings.jsonl.gz', 'rt') as source, \
            gzip.open(ac / 'Archives.csv.gz', 'wt') as dest:
        writer = csv.DictWriter(dest, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for line in source:
            record = json.loads(line)
            assert record['url'] not in seen, record['url']
            seen.add(record['url'])
            writer.writerow(record)
            if record['status'] in ['ok', 'partial']:
                members += len(record['members'])
    assert seen == wanted, 'Archive target coverage differs'
    assert members == pq.ParquetFile(ac / 'Members.parquet').metadata.num_rows
    print(f'Reclassifying {members:,} archive members.', flush=True)
    sys.path.insert(0, str(ROOT / 'work/stack-and-reuse'))
    import common
    common.O = output
    runpy.run_path(str(ROOT / 'work/stack-and-reuse/archive_summary.py'))

    hs = output / 'Header-sniff'
    targets = pd.read_csv(hs / 'Targets.csv.gz')
    wanted = set(targets.loc[targets.file_size_bytes > 0, 'url'])
    seen = set()
    with gzip.open(cache / 'Header-sniff/Headers.jsonl.gz', 'rt') as source, \
            (hs / 'Headers.jsonl').open('w') as dest:
        for line in source:
            record = json.loads(line)
            assert record['url'] not in seen, record['url']
            seen.add(record['url'])
            dest.write(line)
    assert seen == wanted, 'Header target coverage differs'
    import header_sniff
    header_sniff.header_summary()
    report = {
        'scope': 'Offline reclassification of fresh retrieval; not an exact reconstruction of the October 2 extraction.',
        'cache_checksums_verified': True,
        'archive_targets': len(set(pd.read_csv(ac / 'Archives.csv.gz').url)),
        'archive_member_rows': members,
        'header_targets': len(seen),
        'archive_flag_changes': compare_flags(ROOT / 'data/late/Dataset-archive-flags.csv', ac / 'Dataset-before-after.csv.gz'),
        'header_flag_changes': compare_flags(ROOT / 'data/late/Dataset-header-flags.csv', hs / 'Dataset-headers.csv.gz'),
    }
    expected = json.loads((inputs / 'Expected-reclassification.json').read_text())
    assert report == expected, 'Reclassification differs from the checked fresh retrieval'
    (output / 'Validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
