"""Read-only SQL entry point for the tracked cohort and optional full snapshot."""
from pathlib import Path
import argparse, duckdb
R = Path(__file__).resolve().parents[1]

def connect():
    db = R / 'data/PRIDE.duckdb'
    c = duckdb.connect(str(db), read_only=True) if db.exists() else duckdb.connect()
    cohort = str(R/'data/Study-cohort.parquet').replace("'", "''")
    c.execute(f"CREATE TEMP VIEW study AS SELECT * FROM read_parquet('{cohort}')")
    for alias, path in [('g', R/'outputs/Global-proteomics-study/Global-proteomics.duckdb'), ('lit', R/'outputs/Global-proteomics-study/Literature-assets/Literature-enrichment.duckdb')]:
        if path.exists():
            c.execute(f"ATTACH '{str(path).replace(chr(39), chr(39)*2)}' AS {alias} (READ_ONLY)")
    return c

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('--sql');p.add_argument('--file',type=Path);p.add_argument('--check',action='store_true');a=p.parse_args()
    c=connect()
    if a.check:
        n,unique,geo,known=c.execute('SELECT count(*),count(distinct accession),count(country),count(standards) FROM study').fetchone()
        assert (n,unique,geo,known)==(36249,36249,28225,36249),(n,unique,geo,known)
        assert c.execute('SELECT count(*) FROM study WHERE standards != nongenerated_expanded').fetchone()[0] == 0
        for t in ['pride_projects','pride_files','pride_bronze_projects']:
            if (R/'data/PRIDE.duckdb').exists():
                print(t,c.execute('SELECT count(*) FROM '+t).fetchone()[0])
        print('Cohort check passed:',n,'deposits;',geo,'country resolved.')
    else:
        sql=a.file.read_text() if a.file else a.sql
        if not sql:p.error('Provide --sql, --file or --check')
        print(c.execute(sql).df().to_string(index=False))
