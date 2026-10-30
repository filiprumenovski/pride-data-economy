"""Rebuild database-date selection and dependent figures without unrelated analyses."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--cache-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--font-dir', type=Path)
    a = p.parse_args()
    out = a.output.resolve()
    if out.is_relative_to(ROOT):
        raise ValueError('Choose an output directory outside the repository')
    out.mkdir(parents=True, exist_ok=False)
    logs = out / 'logs'; logs.mkdir()
    rebuilt = out / 'rebuilt'; rebuilt.mkdir()
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MPLBACKEND': 'Agg'}
    stages = []

    def run(label, script, *args, render=False, offset=0):
        start = time.monotonic()
        call_env = dict(env)
        if render:
            call_env.update(GAL_JOURNAL='1', GAL_FROZEN_DATA='0', GAL_PANEL_OFFSET=str(offset))
        print('Running:', label, flush=True)
        with (logs / (label + '.log')).open('w') as log:
            r = subprocess.run([sys.executable, str(ROOT / 'scripts' / script), *map(str, args)],
                               cwd=ROOT, env=call_env, stdout=log, stderr=subprocess.STDOUT)
        stages.append({'stage': label, 'exit_code': r.returncode, 'seconds': round(time.monotonic()-start, 2)})
        (out / 'Stages.json').write_text(json.dumps(stages, indent=2)+'\n')
        if r.returncode:
            raise RuntimeError(f'{label} failed; see {logs}')

    run('caches', 'release_caches.py', '--cache-dir', a.cache_dir.resolve(), '--output', out/'caches')
    for name in ['Archive-contents', 'Header-sniff']:
        (rebuilt / name).symlink_to(out/'caches'/name, target_is_directory=True)
    options = ['--source', a.source.resolve(), '--output', rebuilt]
    if a.font_dir:
        options += ['--font-dir', a.font_dir.resolve()]
    for module in ['layers', 'tool_outputs', 'deposit_habits', 'arsenal_time',
                   'arsenal_databases_storage', 'arsenal_clocks']:
        run(module, 'release_full.py', *options, module)
    run('database-figure', 'release_full.py', *options, 'gallery3_databases', render=True, offset=1)
    run('clock-figure', 'release_full.py', *options, 'gallery', 'fig_clocks', render=True)
    clocks = json.loads((rebuilt/'Arsenal/H11-clocks/Summary.json').read_text())['clocks']
    table = pd.DataFrame([{'interval': k, **{key: str(v) if isinstance(v, list) else v
                                           for key, v in row.items()}} for k, row in clocks.items()])
    tables = rebuilt/'Supplementary-Tables'; tables.mkdir()
    table.to_csv(tables/'Table-S15.csv', index=False)


if __name__ == '__main__':
    main()
