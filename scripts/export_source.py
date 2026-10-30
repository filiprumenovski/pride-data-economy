"""Export the reduced reconstruction bundle from the author's archived databases.

Readers use the checksummed release bundle; this optional export needs the original
private snapshot and the historical MassIVE seed list. Contact strings stay local.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def quoted(path):
    return str(path).replace("'", "''")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study-root', type=Path, required=True)
    p.add_argument('--massive-list', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); study, output = a.study_root.resolve(), a.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    c = duckdb.connect(str(study / 'data/PRIDE.duckdb'), read_only=True)
    c.execute("set threads=1; set memory_limit='3GB'; set timezone='UTC'")
    c.execute(f"attach '{quoted(study / 'outputs/Global-proteomics-study/Global-proteomics.duckdb')}' as g (read_only)")
    c.execute(f"attach '{quoted(study / 'outputs/Global-proteomics-study/Literature-assets/Literature-enrichment.duckdb')}' as lit (read_only)")
    c.execute(f"create temp view study as select * from read_parquet('{quoted(study / 'data/Study-cohort.parquet')}')")
    manifest = json.loads((ROOT / 'data/reconstruction/Source-manifest.json').read_text())
    original = pd.read_parquet(study / 'outputs/Stack-and-reuse/Deposit-features.parquet').set_index('accession')
    public = pd.read_parquet(ROOT / 'data/early/Deposit-features.parquet').set_index('accession')
    old = original.pi_group.reindex(public.index)
    known = old.notna()
    assert known.equals(public.contributor_group_known)
    # Preserve order relative to the three historical singleton fallbacks.
    order = pd.DataFrame({'public_group': public.pi_group, 'known': known})
    order['sampling_group'] = [None if pd.isna(v) else ('0' if v < 'P' else 'a0' if v < 'acc:' else 'f') + public.pi_group.loc[k]
                               for k, v in old.items()]
    for prefix in ['', 'acc:', 'solo:']:
        assert old.fillna(pd.Series(prefix + old.index, index=old.index)).sort_values().index.tolist() == order.sampling_group.fillna(pd.Series(prefix + old.index, index=old.index)).sort_values().index.tolist()
    replace = dict(zip(old[known], order.loc[known, 'sampling_group']))
    for entry in manifest['files']:
        dest = output / entry['path']; dest.parent.mkdir(parents=True, exist_ok=True)
        if 'query' in entry:
            c.execute(f"copy ({entry['query']}) to '{quoted(dest)}' (format parquet, compression zstd)")
        elif entry['path'] == 'Group-order.csv':
            order.to_csv(dest, index_label='accession')
        elif entry['path'] == 'Geography-evidence.parquet':
            subprocess.run([sys.executable, str(ROOT / 'scripts/export_geography_evidence.py'), '--pride-db',
                str(study / 'data/PRIDE.duckdb'), '--massive-list', str(a.massive_list), '--output', str(dest)], check=True)
        elif entry['path'].startswith('inputs/'):
            src = study / 'outputs/Stack-and-reuse' / Path(entry['path']).relative_to('inputs')
            if dest.name == 'H1-dataset-scale.csv.gz':
                df = pd.read_parquet(ROOT / 'data/early/Dataset-scale.parquet')
                df[['accession', 'runs_lb', 'runs_exact', 'status']].to_csv(dest, index=False)
            elif '.csv' in dest.name:
                df = pd.read_csv(src)
                for col in ['pi_group', 'grp', 'g']:
                    if col in df:
                        df[col] = df[col].map(lambda v: replace.get(v, v))
                df.to_csv(dest, index=False)
            else:
                shutil.copy2(src, dest)
        else:
            raise ValueError(entry['path'])
        entry.update(bytes=dest.stat().st_size, sha256=hashlib.file_digest(dest.open('rb'), 'sha256').hexdigest())
    (output / 'Source-manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print('Exported source bundle. Inspect its manifest before adopting it as a new release snapshot.')


if __name__ == '__main__':
    main()
