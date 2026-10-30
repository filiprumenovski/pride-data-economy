"""Read the header row of deposited tables that our filename rules could not identify.

Targets: top-level, non-generated files classified as `unidentified_tables` (txt/tsv/csv, also
.gz, and xlsx), at most 20 per dataset (fixed seed). Only the first 32 KiB of a text table is
requested (HTTP range); for an xlsx (itself a zip) only the central directory, the start of the
shared-strings part and the start of the first worksheet. Legacy .xls and parquet are skipped
and counted.

Each header is then matched against
  (A) PRIDE's own header definitions (`content_verification.required_columns` in
      reference/pride-submission-validation.yaml v2.2.2): all required columns present;
  (B) column signatures of tool outputs that the YAML does not define (matrices, protein tables
      re-saved after processing), listed in SIGNATURES below. The list was written first, then
      extended once after a 16-file pilot (seed 5) showed misses: DIA-NN gene matrices, MaxQuant
      sites/mzRange tables, PD protein tables, PLGS, Crux/Tide and mzTab sections saved as .txt.
      After the full run, a model-assisted precision review (Attacks/headers/, 389 tables, 382 correct)
      led to three changes: the mzTab-section rule was dropped (0/7: ProteomeXchange submission
      summaries), statistics tables (limma/Perseus) are no longer counted as tool output, and the
      run-column rule is labelled as a DIA-NN/Spectronaut-style matrix.
Output: Header-sniff/Headers.jsonl (resumable) and summaries written by header_summary().

  python header_sniff.py            # fetch, then summarise
  python header_sniff.py --summary  # summarise existing headers only
"""
import argparse, asyncio, json, re, struct, zlib
import numpy as np, pandas as pd, yaml
from common import R, O, SEED
from layer_rules import RULES
from archive_contents import https_url
from fastfetch import run_jsonl
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Header-sniff'; OUT.mkdir(exist_ok=True)
HEAD = OUT / 'Headers.jsonl'
PER_DATASET = 20
Y = yaml.safe_load(open(R / 'work/stack-and-reuse/reference/pride-submission-validation.yaml'))

# (B) signatures: every listed column (lower case) must be present
SIGNATURES = {
    'DIA-NN matrix': [['protein.group', 'protein.ids', 'genes'], ['protein.group', 'protein.names'], ['precursor.id', 'stripped.sequence']],
    'DIA-NN report (not PRIDE-complete)': [['file.name', 'protein.group', 'precursor.id']],
    'Spectronaut report/pivot': [['pg.proteingroups'], ['pg.proteinaccessions'], ['eg.precursorid'], ['r.filename', 'pg.genes']],
    'MaxQuant table': [['majority protein ids'], ['razor + unique peptides'], ['leading razor protein', 'raw file'], ['protein ids', 'peptide counts (all)'],
                       ['positions within proteins', 'localization prob'], ['raw file', 'peaks / da']],
    'Proteome Discoverer export': [['master', 'accession'], ['annotated sequence', '# protein groups'],
                                   ['confidence', 'annotated sequence'], ['identifying node']],
    'FragPipe output': [['protein id', 'protein probability'], ['spectrum file', 'modified peptide'], ['protein id', 'combined spectral count']],
    'PEAKS export': [['-10lgp']],
    'Scaffold export': [['exclusive unique peptide count'], ['protein identification probability']],
    'Mascot export': [['prot_acc', 'pep_seq'], ['prot_hit_num']],
    'PLGS export': [['protein.key', 'protein.accession'], ['peptide.seq', 'protein.accession']],
    'Crux/Tide output': [['file_idx', 'scan', 'spectrum precursor m/z']],
}


def pride_defs():
    out = []
    for t in Y['tools']:
        for f in t.get('files', []):
            cv = f.get('content_verification') or {}
            if 'formats' in cv:
                cv = cv['formats'].get('tsv') or list(cv['formats'].values())[0]
            if cv.get('required_columns'):
                out.append((f"{t['name']}: {f['name']}", [c.lower() for c in cv['required_columns']]))
    return out


PRIDE = pride_defs()


def targets():
    case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
    t = connect().execute(f"""
    with f as (select project_accession a, file_name, lower(file_name) n,
      lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
      coalesce(file_category_value, '') cat, file_size_bytes, primary_public_location_value ploc
      from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
      and project_accession in (select accession from study))
    select a accession, file_name, file_size_bytes, ploc, {case} layer from f""").df()
    t = t[t.layer == 'unidentified_tables'].copy()
    k = t.file_name.str.lower().str.extract(r'\.(txt|tsv|csv|xlsx|xls|parquet)(\.gz)?$')
    t['kind'] = k[0].fillna('other') + k[1].fillna('')
    t['url'] = t.ploc.map(lambda x: https_url(x or ''))
    t = t.sort_values(['accession', 'file_name']).sample(frac=1, random_state=SEED).groupby('accession').head(PER_DATASET)  # sort first: query order is not fixed
    return t


def split_header(text):
    for line in text.replace('\r\n', '\n').replace('\r', '\n').split('\n'):
        s = line.strip().lstrip('﻿')
        if not s or s.startswith('#'):
            continue
        d = '\t' if '\t' in s else (',' if s.count(',') >= s.count(';') else ';')
        return [c.strip().strip('"').strip("'") for c in s.split(d)][:500]
    return []


async def read_text(get, url, gz):
    data, _ = await get(url, 'bytes=0-32767')
    if gz:
        data = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(data)
    return split_header(data.decode('utf-8', 'replace'))


async def zip_entries(get, url):
    tail, total = await get(url, 'bytes=-65558')
    e = tail.rfind(b'PK\x05\x06')
    n, cd_size, cd_off = struct.unpack_from('<HII', tail, e + 10)
    ts = total - len(tail)
    cd = tail[cd_off - ts:cd_off - ts + cd_size] if cd_off >= ts else (await get(url, f'bytes={cd_off}-{cd_off + cd_size - 1}'))[0]
    out, p = {}, 0
    while p + 46 <= len(cd) and cd[p:p + 4] == b'PK\x01\x02':
        method, = struct.unpack_from('<H', cd, p + 10); csize, = struct.unpack_from('<I', cd, p + 20)
        nlen, xlen, clen = struct.unpack_from('<HHH', cd, p + 28); off, = struct.unpack_from('<I', cd, p + 42)
        out[cd[p + 46:p + 46 + nlen].decode('utf-8', 'replace')] = (method, csize, off)
        p += 46 + nlen + xlen + clen
    return out, tail, ts


async def zip_part(get, url, entry, tail, ts, nbytes=131072):
    """Start of one member, decompressed. One request: the local header and the data are read
    together (the local header's name/extra lengths are not known in advance, so 1 KiB slack)."""
    method, csize, off = entry
    want = 30 + 1024 + min(csize, nbytes)
    if off >= ts:                                  # small workbooks: already inside the tail read
        buf = tail[off - ts:off - ts + want]
    else:
        buf, _ = await get(url, f'bytes={off}-{off + want - 1}')
    nlen, xlen = struct.unpack_from('<HH', buf, 26)
    data = buf[30 + nlen + xlen:30 + nlen + xlen + min(csize, nbytes)]
    return zlib.decompressobj(-15).decompress(data) if method == 8 else data


async def read_xlsx(get, url):
    ents, tail, ts = await zip_entries(get, url)
    sheet = 'xl/worksheets/sheet1.xml' if 'xl/worksheets/sheet1.xml' in ents else \
        next((k for k in sorted(ents) if k.startswith('xl/worksheets/sheet')), None)
    if sheet is None:
        return [], 'no_sheet'
    has_ss = 'xl/sharedStrings.xml' in ents
    parts = await asyncio.gather(zip_part(get, url, ents[sheet], tail, ts),
                                 *([zip_part(get, url, ents['xl/sharedStrings.xml'], tail, ts)] if has_ss else []))
    xml = parts[0].decode('utf-8', 'replace')
    m = re.search(r'<(?:\w+:)?row[^>]*>(.*?)</(?:\w+:)?row>', xml, re.S)
    if not m:
        return [], 'no_row'
    cells = re.findall(r'<(?:\w+:)?c ([^>]*?)(?:/>|>(.*?)</(?:\w+:)?c>)', m.group(1), re.S)
    shared = []
    if has_ss:
        ss = parts[1].decode('utf-8', 'replace')
        shared = [''.join(re.findall(r'<(?:\w+:)?t[^>]*>([^<]*)</(?:\w+:)?t>', si)) for si in re.findall(r'<(?:\w+:)?si>(.*?)</(?:\w+:)?si>', ss, re.S)]
    cols, partial = [], False
    for a, body in cells:
        v = re.search(r'<(?:\w+:)?v>([^<]*)</', body or ''); t = re.search(r'<(?:\w+:)?t[^>]*>([^<]*)</', body or '')
        if 't="s"' in a and v:
            i = int(v.group(1)); partial |= i >= len(shared); cols.append(shared[i] if i < len(shared) else '')
        elif t:
            cols.append(t.group(1))
        elif v:
            cols.append(v.group(1))
    return [c.strip() for c in cols][:500], 'partial' if partial else 'ok'


async def fetch_one(get, row):
    try:
        if row.kind in ('txt', 'tsv', 'csv', 'txt.gz', 'tsv.gz', 'csv.gz'):
            cols, st = await read_text(get, row.url, row.kind.endswith('.gz')), 'ok'
        elif row.kind == 'xlsx':
            cols, st = await read_xlsx(get, row.url)
        else:
            cols, st = [], 'skipped_format'
    except Exception as ex:
        cols, st = [], f'error: {type(ex).__name__}: {str(ex)[:60]}'
    return {'url': row.url, 'accession': row.accession, 'file_name': row.file_name, 'kind': row.kind, 'status': st, 'columns': cols}


RUNCOL = re.compile(r'\.(raw|d|wiff|wiff2|mzml|dia)$', re.I)


def classify(cols):
    s = {c.lower().strip() for c in cols}
    a = [name for name, req in PRIDE if all(r in s for r in req)]
    b = [name for name, alts in SIGNATURES.items() if any(all(x in s for x in alt) for alt in alts)]
    # per-run quantity matrix: an identifier column plus >= 2 columns named after raw runs (DIA-NN/Spectronaut style)
    if not b and s & {'genes', 'protein.group', 'protein.ids', 'precursor.id', 'modified.sequence', 'stripped.sequence'} \
            and sum(bool(RUNCOL.search(c.strip())) for c in cols) >= 2:
        b = ['Per-run quantity matrix (DIA-NN/Spectronaut style)']
    return a, b


def header_summary():
    H = pd.DataFrame([json.loads(l) for l in HEAD.open()]).drop_duplicates('url', keep='last')
    H = H[H.url.isin(pd.read_csv(OUT / 'Targets.csv.gz').url)]   # only the fixed sample
    cl = H['columns'].apply(classify)
    H['pride_def'] = cl.str[0].str.join('; '); H['signature'] = cl.str[1].str.join('; ')
    H['read'] = H.status.isin(['ok', 'partial']) & (H['columns'].str.len() > 0)
    H['identified'] = (H.pride_def != '') | (H.signature != '')
    H.drop(columns=['columns']).to_csv(OUT / 'Header-matches.csv.gz', index=False)
    L = pd.read_parquet(O / 'Layers/Dataset-layers.parquet').set_index('accession')
    forms = O / 'Archive-contents/Dataset-before-after.csv.gz'
    D = pd.read_csv(forms).set_index('accession')
    ds = H[H.read].groupby('accession').agg(any_identified=('identified', 'any'), any_pride=('pride_def', lambda x: (x != '').any()),
                                            files_read=('url', 'size'))
    D = D.join(ds)
    D['ad_hoc_only'] = D.adhoc_only_after
    sub = D[D.ad_hoc_only & D.files_read.notna()]
    out = {'files_targeted': int(len(H)), 'files_read': int(H.read.sum()),
           'status': H.status.str.split(':').str[0].value_counts().to_dict(),
           'files_identified_pct_of_read': round(100 * H[H.read].identified.mean(), 1),
           'by_signature_files': H[H.read].signature.str.split('; ').explode().value_counts().drop('', errors='ignore').to_dict(),
           'by_pride_definition_files': H[H.read].pride_def.str.split('; ').explode().value_counts().drop('', errors='ignore').to_dict(),
           'ad_hoc_only_datasets_with_a_read_table': int(len(sub)),
           'of_which_tool_output_by_header_pct': round(100 * sub.any_identified.fillna(False).mean(), 1),
           'of_which_meets_a_pride_header_definition_pct': round(100 * sub.any_pride.fillna(False).mean(), 1)}
    yr = sub.groupby('release_year').any_identified.mean().mul(100).round(1)
    yr.to_csv(OUT / 'Ad-hoc-identified-by-year.csv')
    D['recognisable_with_headers'] = D.recognisable_after | D.any_identified.fillna(False).astype(bool)
    D['adhoc_after_headers'] = D.adhoc_only_after & ~D.any_identified.fillna(False).astype(bool)
    by = D.groupby('release_year')[['recognisable_before', 'recognisable_after', 'recognisable_with_headers', 'adhoc_only_before',
                                    'adhoc_only_after', 'adhoc_after_headers']].mean().mul(100).round(1)
    by.to_csv(OUT / 'Recognisable-by-year-three-ways.csv')
    D[['release_year', 'dia_expanded', 'recognisable_with_headers', 'adhoc_after_headers', 'any_identified', 'any_pride']].to_csv(OUT / 'Dataset-headers.csv.gz')
    (OUT / 'Summary.json').write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=1)); print(by.to_string()); print(yr.to_string())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--summary', action='store_true'); ap.add_argument('--workers', type=int, default=96)
    a = ap.parse_args()
    if not a.summary:
        t = targets(); t.drop(columns=['ploc']).to_csv(OUT / 'Targets.csv.gz', index=False)
        done = {json.loads(l)['url'] for l in HEAD.open()} if HEAD.exists() else set()
        todo = t[t.url.notna() & ~t.url.isin(done) & (t.file_size_bytes > 0)]
        print(len(t), 'targets;', len(done), 'done;', len(todo), 'to fetch', flush=True)
        run_jsonl(todo.itertuples(index=False), fetch_one, HEAD, concurrency=a.workers)
    header_summary()


if __name__ == '__main__':
    main()
