"""Recompute cohort membership and annotation flags from reduced source metadata."""
import argparse
import json
from pathlib import Path
import sys
import warnings
import duckdb
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'work/stack-and-reuse'))
from common import SOFTWARE
from release_full import prepare


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    a = p.parse_args(); source, output = a.source.resolve(), a.output.resolve(); prepare(source, output)
    f = pd.read_parquet(ROOT / 'data/early/Deposit-features.parquet').set_index('accession')
    projects = pd.read_parquet(source / 'pride_projects.parquet').set_index('project_accession')
    dates = pd.to_datetime(projects.publication_date)
    assert set(projects.index[dates.dt.year.between(2013, 2025)]) == set(f.index)
    assert pd.to_datetime(projects.reindex(f.index).submission_date).eq(f.submission_date).all()
    assert dates.reindex(f.index).eq(f.release_date).all()
    t = pd.read_parquet(source / 'study_text.parquet').set_index('accession').reindex(f.index).fillna('')
    annotation = t.software.str.lower() + ' ' + t.processing.str.lower()
    result = {'cohort_accessions': len(f), 'submission_and_release_dates_match': True, 'software_flags': {}}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        for key, pattern in SOFTWARE.items():
            values = annotation.str.contains(pattern, regex=True)
            assert values.eq(f['sw_' + key]).all(), key
            result['software_flags'][key] = int(values.sum())
    sp = pd.read_parquet(source / 'pride_species.parquet')
    sp['human'] = sp.species_accession.isin(['NEWT:9606', 'NCBITaxon:9606'])
    species = sp.groupby('project_accession').agg(human=('human', 'any'), species_n=('species_name', 'nunique')).reindex(f.index)
    for col in species:
        assert (species[col].eq(f[col]) | (species[col].isna() & f[col].isna())).all(), col
    result['species_flags_match'] = True
    c = duckdb.connect()
    pp = pd.read_parquet(source / 'lit.project_publications.parquet')
    works = pd.read_parquet(source / 'lit.openalex_works.parquet')
    c.register('pp', pp); c.register('works', works)
    first = c.execute('''select accession,work_id from pp join works using(work_id)
        qualify row_number() over(partition by accession order by publication_date asc nulls last,work_id)=1''').df().set_index('accession')
    linked = first.work_id.reindex(f.index)
    assert (linked.eq(f.work_id) | (linked.isna() & f.work_id.isna())).all()
    result['primary_publication_links_match'] = True
    c.close()
    (output / 'Metadata-validation.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
