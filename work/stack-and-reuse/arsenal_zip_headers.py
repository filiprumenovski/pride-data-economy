"""H10: header rows of text tables inside 'tables only' results zips, for datasets still classed
'ad hoc tables only' after top-level header reading (deposit_habits.py).

Each zip's central directory is read to locate the member; the first 64 KiB of the member is
inflated and its header row matched with the same rules as header_sniff.py (PRIDE header
definitions + tool signatures). At most 20 tables per dataset (sorted, then seeded sample).
Nested xlsx inside zips cannot be read this way and are not targeted.
Output: Arsenal/H10-zip-headers/ (Zip-headers.jsonl gitignored, Summary.json, CSVs).
"""
import json
import pandas as pd, duckdb
from common import O, SEED
from header_sniff import split_header, classify
from version_scan import central_directory, member_head, ungz
from fastfetch import run_jsonl

OUT = O / 'Arsenal/H10-zip-headers'; OUT.mkdir(parents=True, exist_ok=True)
REC = OUT / 'Zip-headers.jsonl'
AC = O / 'Archive-contents'


def targets():
    A = pd.read_csv(AC / 'Archive-content.csv.gz', usecols=['url', 'accession', 'content'])
    F = pd.read_csv(O / 'Habits/Dataset-form.csv.gz', usecols=['accession', 'form'])
    A = A.merge(F, on='accession'); A = A[(A.content == 'tables only') & (A.form == 'ad hoc tables only')]
    c = duckdb.connect(); c.execute("set threads=1; set memory_limit='1500MB'"); c.register('tab', A[['url', 'accession']])
    t = c.execute(r"""select accession, url, member from read_parquet(?) p join tab using (url)
        where regexp_matches(lower(member), '\.(txt|tsv|csv)(\.gz)?$') and not regexp_matches(member, '(^|/)__MACOSX/|(^|/)\._')""",
                  [str(AC / 'Members.parquet')]).df()
    t = t.sort_values(['accession', 'url', 'member']).sample(frac=1, random_state=SEED).groupby('accession').head(20)
    return t, int(A.accession.nunique())


async def work(get, row):
    rec = {'url': row.url, 'member': row.member, 'accession': row.accession}
    try:
        cd = await central_directory(get, row.url)
        if row.member not in cd:
            return rec | {'status': 'member_not_found', 'columns': []}
        data = ungz(await member_head(get, row.url, cd[row.member]), row.member)
        return rec | {'status': 'ok', 'columns': split_header(data.decode('utf-8', 'replace'))}
    except Exception as ex:
        return rec | {'status': f'error: {type(ex).__name__}: {str(ex)[:60]}', 'columns': []}


def summary(n_ds):
    H = pd.DataFrame([json.loads(l) for l in REC.open()])
    cl = H['columns'].apply(classify)
    H['pride_def'] = cl.str[0].str.join('; '); H['signature'] = cl.str[1].str.join('; ')
    H['identified'] = (H.pride_def != '') | (H.signature != '')
    H.drop(columns=['columns']).to_csv(OUT / 'Zip-header-matches.csv.gz', index=False)
    ds = H[H.status == 'ok'].groupby('accession').identified.any()
    ds.rename('zip_header_tool').to_csv(OUT / 'Dataset-zip-headers.csv')
    out = {'datasets_ad_hoc_with_tables_only_zip': n_ds, 'datasets_with_text_tables_read': int(len(ds)),
           'tables_read': int((H.status == 'ok').sum()), 'status': H.status.str.split(':').str[0].value_counts().to_dict(),
           'datasets_identified_pct': round(100 * ds.mean(), 1), 'datasets_identified': int(ds.sum()),
           'by_signature_tables': H[H.identified].signature.str.split('; ').explode().value_counts().drop('', errors='ignore').to_dict()}
    (OUT / 'Summary.json').write_text(json.dumps(out, indent=2)); print(json.dumps(out, indent=1))


if __name__ == '__main__':
    t, n_ds = targets(); t.to_csv(OUT / 'Targets.csv.gz', index=False)
    done = {(r['url'], r['member']) for r in map(json.loads, REC.open())} if REC.exists() else set()
    todo = [r for r in t.itertuples(index=False) if (r.url, r.member) not in done]
    print(len(t), 'tables in', t.accession.nunique(), 'datasets;', len(todo), 'to fetch', flush=True)
    run_jsonl(todo, work, REC, concurrency=96)
    summary(n_ds)
