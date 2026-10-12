"""Reduce original PRIDE contacts and the historical MassIVE seed cache to public evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys

import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'work/stack-and-reuse'))
from geography_rules import norm, address_pattern, digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pride-db', type=Path, required=True)
    p.add_argument('--massive-list', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    c = duckdb.connect(str(args.pride_db), read_only=True)
    b = c.execute('select b.accession,b.submitters,b.labPIs,p.countries from pride_bronze_projects b join pride_projects p on p.project_accession=b.accession').df()
    rows = []

    def add(repo, acc, role, person, project=None):
        aff = person.get('affiliation') or person.get('institution') or ''
        email = (person.get('email') or '').lower().strip(); pc = person.get('country')
        ma = address_pattern.search(re.sub(r'\s*\(\s*lab head\s*\)\s*$', '', aff, flags=re.I))
        ph = (digest(email) if '@' in email else digest(person.get('orcid'))
              if re.search(r'\d{4}-\d{4}-\d{4}-', person.get('orcid') or '') else None)
        rows.append(dict(repository=repo, accession=acc, role=role,
            affinity=norm(aff) if norm(aff) not in ('', 'na', 'null', 'none', 'unknown') else None,
            email=email if '@' in email else None, person=ph,
            raw_country=pc or (project if role == 'submitter' else None),
            person_raw_country=pc, project_country=project,
            affiliation_suffix=ma.group(1) if ma else None))

    for r in b.itertuples():
        for role, people in [('submitter', r.submitters), ('pi', r.labPIs)]:
            for person in people:
                add('PRIDE', r.accession, role, person, r.countries[0] if len(r.countries) else None)
    live = json.loads(args.massive_list.read_text(errors='replace'))
    assert len(live['row_data']) == int(live['total_rows'])
    for r in live['row_data']:
        for person in r.get('pis', []) or []:
            add('MassIVE', r['dataset'], 'pi', person)
    d = pd.DataFrame(rows)
    for col, prefix in [('affinity', 'A'), ('email', 'E'), ('person', 'P')]:
        mapping = {x: f'{prefix}{i:06d}' for i, x in enumerate(sorted(d[col].dropna().unique()))}
        d[col] = d[col].map(mapping)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    d.to_parquet(args.output, index=False)
    print(json.dumps({'rows': len(d), 'PRIDE': int(d.repository.eq('PRIDE').sum()),
        'MassIVE_seed_rows': int(d.repository.eq('MassIVE').sum()),
        'sha256': hashlib.file_digest(args.output.open('rb'), 'sha256').hexdigest()}, indent=2))


if __name__ == '__main__':
    main()
