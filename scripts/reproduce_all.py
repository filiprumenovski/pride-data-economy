"""Rebuild analyses, figure statistics, 600-dpi figures and table CSVs offline.

Requires the Source-metadata and Extraction-caches release bundles. Outputs go to
one new directory; the repository's frozen manuscript artifacts are never edited.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--cache-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--font-dir', type=Path)
    a = p.parse_args()
    source, cache, output = a.source.resolve(), a.cache_dir.resolve(), a.output.resolve()
    if output == ROOT or output.is_relative_to(ROOT):
        raise ValueError('Choose an output directory outside the repository')
    output.mkdir(parents=True, exist_ok=False)
    logs = output / 'logs'; logs.mkdir()
    rebuilt = output / 'rebuilt'
    env = {**os.environ, 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1', 'MPLBACKEND': 'Agg', 'GAL_FROZEN_DATA': '0'}
    stages = []

    def run(label, script, *args):
        start = time.monotonic()
        print('Running:', label, flush=True)
        with (logs / (label + '.log')).open('w') as log:
            r = subprocess.run([sys.executable, str(ROOT / 'scripts' / script), *map(str, args)], cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        stages.append({'stage': label, 'exit_code': r.returncode, 'seconds': round(time.monotonic() - start, 2)})
        (output / 'Stages.json').write_text(json.dumps(stages, indent=2) + '\n')
        if r.returncode:
            raise RuntimeError(f'{label} failed; see {logs / (label + ".log")}')

    def analysis(module, *args):
        options = ['--source', source, '--output', rebuilt]
        if a.font_dir:
            options += ['--font-dir', a.font_dir.resolve()]
        run(module, 'release_full.py', *options, module, *args)

    run('metadata', 'release_metadata.py', '--source', source, '--output', rebuilt)
    run('geography', 'release_geography.py', '--source', source, '--output', rebuilt)
    run('early-validation', 'release_early.py', 'validate', '--data-dir', ROOT / 'data/early', '--output-dir', output / 'early')
    run('cache-reclassification', 'release_caches.py', '--cache-dir', cache, '--output', output / 'caches')
    for name in ['Archive-contents', 'Header-sniff']:
        (rebuilt / name).symlink_to(output / 'caches' / name, target_is_directory=True)
    for module in ['layers', 'tool_outputs', 'deposit_habits', 'arsenal_scale_instruments', 'attack_runs',
                   'arsenal_time', 'attack_release_timing', 'arsenal_databases_storage',
                   'arsenal_strategies_software_docs', 'arsenal_clocks', 'software_geography',
                   'figure_data', 'standards_by_stack', 'within_group']:
        analysis(module)
    run('late-validation', 'release_late.py', 'validate', '--data', ROOT / 'data/late', '--output', output / 'late')
    run('late-support-validation', 'release_late_support.py', '--data', ROOT / 'data/late', '--output', output / 'late-support')
    for relative in ['Reuse-mentions', 'Attacks/reuse-mentions']:
        shutil.copytree(output / 'late' / relative, rebuilt / relative)
    analysis('compose_journal')
    analysis('build_supp_tables', '--csv-only')
    run('numerical-comparison', 'compare_reconstruction.py',
        '--reference', ROOT / 'outputs/Stack-and-reuse', '--rebuilt', rebuilt,
        '--expected-list', ROOT / 'data/reconstruction/Expected-tables.json',
        '--report', output / 'Numerical-comparison.json', '--require-match')
    print(f'Completed {len(stages)} stages. Rebuilt artifacts: {rebuilt}', flush=True)


if __name__ == '__main__':
    main()
