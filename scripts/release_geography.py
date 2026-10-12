"""Reconstruct geography from deidentified contact evidence and validate the cohort.

Affiliation equality and matched country suffixes are retained; full affiliation and
contact strings are excluded. This tests inference and grouping, not raw ingestion.
"""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'work/stack-and-reuse'))
from geography_rules import country, address_country
from release_full import prepare


def reconstruct(d):
    aff = defaultdict(lambda: defaultdict(set))
    email = defaultdict(lambda: defaultdict(set))
    d = d.astype(object).where(d.notna(), None)
    for r in d.itertuples():
        ct = (country(r.person_raw_country) if r.repository == 'MassIVE' else
              (country(r.person_raw_country) or country(r.project_country)) if r.role == 'submitter'
              else address_country(r.affiliation_suffix))
        if ct and r.affinity:
            aff[r.affinity][ct].add(r.repository + ':' + r.accession)
        if ct and r.email and (r.role == 'submitter' or r.repository == 'MassIVE'):
            email[r.email][ct].add(r.repository + ':' + r.accession)

    def infer(evidence):
        if len(evidence) == 1:
            ct, accessions = next(iter(evidence.items()))
            if len(accessions) >= 2:
                return ct
        return None

    rows = []
    for r in d.itertuples():
        if r.repository != 'PRIDE':
            continue
        candidates = [(country(r.raw_country), 'reported_country'),
                      (address_country(r.affiliation_suffix), 'country_in_affiliation'),
                      (infer(aff.get(r.affinity, {})), 'exact_affiliation_crossrecord'),
                      (infer(email.get(r.email, {})), 'exact_contact_crossrecord')]
        ct, tier = next(((ct, tier) for ct, tier in candidates if ct), (None, 'unresolved'))
        rows.append({'accession': r.accession, 'role': r.role, 'country': ct, 'tier': tier, 'person': r.person})
    contacts = pd.DataFrame(rows)
    rows = []
    for acc, g in contacts.groupby('accession'):
        pi = g[g.role == 'pi']; sub = g[g.role == 'submitter']
        countries = pi.country.dropna().unique()
        submitters = sub.loc[sub.tier == 'reported_country', 'country'].dropna().unique()
        persons = pi.person.dropna().unique()
        rows.append({'accession': acc, 'country': countries[0] if len(countries) == 1 else None,
                     'submitter_country': submitters[0] if len(submitters) == 1 else None,
                     'pi_geo_inferred_only': not (pi.tier == 'reported_country').any() and len(countries) > 0,
                     'person': persons[0] if len(persons) == 1 else None})
    return contacts, pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); source, output = a.source.resolve(), a.output.resolve()
    prepare(source, output)
    contacts, rebuilt = reconstruct(pd.read_parquet(source / 'Geography-evidence.parquet'))
    expected = pd.read_parquet(ROOT / 'data/early/Deposit-features.parquet')
    x = expected.merge(rebuilt, on='accession', suffixes=('_expected', '_rebuilt'), validate='one_to_one')
    report = {'datasets': len(x), 'scope': __doc__, 'differences': {}}
    for col in ['country', 'submitter_country', 'pi_geo_inferred_only']:
        old, new = x[col + '_expected'], x[col + '_rebuilt']
        bad = ~(old.eq(new) | (old.isna() & new.isna()))
        report['differences'][col] = x.loc[bad, 'accession'].tolist()
    known = x.person.notna()
    report['group_known_mismatches'] = int((known != x.contributor_group_known).sum())
    z = x[known]
    report['group_partition_matches'] = bool(z.groupby('person').pi_group.nunique().eq(1).all()
        and z.groupby('pi_group').person.nunique().eq(1).all())
    rebuilt.to_csv(output / 'Rebuilt-geography.csv', index=False)
    contacts.groupby(['role', 'tier']).size().to_csv(output / 'Geography-evidence-tiers.csv')
    (output / 'Geography-validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({**report, 'differences': {k: len(v) for k, v in report['differences'].items()}}, indent=2))
    assert not any(report['differences'].values()) and not report['group_known_mismatches'] and report['group_partition_matches']


if __name__ == '__main__':
    main()
