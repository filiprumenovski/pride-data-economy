"""Which software versions produced the deposited results? Read from the files themselves.

Version-bearing files (top level, or inside listed results zips):
  mq_parameters   MaxQuant parameters.txt       'Version<TAB>1.6.17.0'
  mqpar           MaxQuant mqpar.xml            <maxQuantVersion>
  diann_log       DIA-NN report log             'DIA-NN 1.8.1 ...'
  fragpipe        fragpipe.workflow             FragPipe version in the header comment
  skyline         Skyline .sky                  software_version="Skyline ..."
  mascot_dat      Mascot result .dat            version= in the header section (x-Mascot files only)
  mzid / mztab    PSI standard files            AnalysisSoftware / MTD software[n]: the tool that WROTE the file
At most 2 files per dataset and type (top-level files first, fixed order). Only the first 64 KiB
of a file is requested; inside a zip, the central directory is read to locate the member, then
the start of the member (inflated). Nothing is downloaded in full.

Derived (descriptive): per tool, version in use by release year; 'version lag' = deposit
submission date minus the date that version FIRST appears anywhere in this archive (internal
frontier; no external release dates are assumed).

  python version_scan.py            # fetch then summarise
  python version_scan.py --summary
"""
import argparse, asyncio, json, re, struct, zlib
import numpy as np, pandas as pd, duckdb
from common import R, O, SEED
from archive_contents import https_url
from fastfetch import run_jsonl
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Arsenal/H9-versions'; OUT.mkdir(parents=True, exist_ok=True)
REC = OUT / 'Versions.jsonl'
AC = O / 'Archive-contents'
PER = 2
TYPE_SQL = r"""case when bn like '%parameters.txt' then 'mq_parameters' when bn like '%mqpar%.xml' then 'mqpar'
  when bn like '%log.txt' and (bn like '%report%' or bn like '%diann%' or bn like '%dia-nn%') then 'diann_log'
  when bn like '%fragpipe.workflow' then 'fragpipe' when bn like '%.sky' then 'skyline'
  when regexp_matches(bn, '\.mzid(\.gz)?$') then 'mzid' when regexp_matches(bn, '\.mztab(\.gz)?$') then 'mztab'
  when regexp_matches(bn, '(^|[^a-z0-9])f\d{5,6}\.dat$') then 'mascot_dat' end"""


def targets():
    c = connect()
    top = c.execute(f"""
    with f as (select project_accession accession, file_name, primary_public_location_value ploc, file_size_bytes,
                      lower(regexp_replace(file_name, '\\.(zip|bz2)$', '', 'i')) bn
               from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
               and project_accession in (select accession from study) and file_size_bytes > 0
               and not regexp_matches(lower(file_name), '\\.(zip|bz2|7z|rar|tgz|tar)$'))
    select accession, file_name, ploc, {TYPE_SQL} as "type" from f""").df().dropna(subset=['type'])
    top['url'] = top.ploc.map(lambda x: https_url(x or '')); top['member'] = None
    con = duckdb.connect(); con.execute("set threads=1; set memory_limit='1500MB'; set preserve_insertion_order=false")
    inz = con.execute(f"""
    with p as (select url, member from read_parquet('{AC / 'Members.parquet'}')     -- cheap name filter first (36M rows)
               where regexp_matches(lower(member), '(parameters\\.txt|mqpar[^/]*\\.xml|log\\.txt|fragpipe\\.workflow|\\.sky|\\.mzid|\\.mztab|\\.dat)(\\.gz)?$')),
    m as (select a.accession, p.url, p.member, lower(regexp_extract(p.member, '[^/]*$')) bn
          from p join read_csv_auto('{AC / 'Archives.csv.gz'}') a using (url))
    select accession, url, member, {TYPE_SQL} as "type" from m""").df().dropna(subset=['type'])
    inz['file_name'] = inz.url.str.rsplit('/', n=1).str[1]
    top['src'], inz['src'] = 'top', 'zip'
    t = pd.concat([top[['accession', 'type', 'src', 'url', 'member', 'file_name']], inz[['accession', 'type', 'src', 'url', 'member', 'file_name']]])
    t = t[t.url.notna()].sort_values(['accession', 'type', 'src', 'file_name', 'member'], na_position='first')
    t = t[~t.member.fillna('').str.contains(r'(^|/)__MACOSX/|(^|/)\._', regex=True)]
    return t.groupby(['accession', 'type']).head(PER).reset_index(drop=True)


async def central_directory(get, url):
    """name -> (method, csize, local_header_offset); handles zip64."""
    tail, total = await get(url, 'bytes=-65558')
    e = tail.rfind(b'PK\x05\x06')
    if e < 0:
        raise RuntimeError('no_zip_directory')
    n, cd_size, cd_off = struct.unpack_from('<HII', tail, e + 10)
    ts = total - len(tail)
    if 0xFFFFFFFF in (cd_size, cd_off) or n == 0xFFFF:
        loc = tail.rfind(b'PK\x06\x07', 0, e)
        z64, = struct.unpack_from('<Q', tail, loc + 8)
        rec = tail[z64 - ts:z64 - ts + 56] if z64 >= ts else (await get(url, f'bytes={z64}-{z64 + 55}'))[0]
        n, cd_size, cd_off = struct.unpack_from('<QQQ', rec, 32)
    cd = tail[cd_off - ts:cd_off - ts + cd_size] if cd_off >= ts else (await get(url, f'bytes={cd_off}-{cd_off + cd_size - 1}'))[0]
    out, p = {}, 0
    while p + 46 <= len(cd) and cd[p:p + 4] == b'PK\x01\x02':
        method, = struct.unpack_from('<H', cd, p + 10); csize, usize = struct.unpack_from('<II', cd, p + 20)
        nlen, xlen, clen = struct.unpack_from('<HHH', cd, p + 28); off, = struct.unpack_from('<I', cd, p + 42)
        flags, = struct.unpack_from('<H', cd, p + 8)
        name = cd[p + 46:p + 46 + nlen].decode('utf-8' if flags & 0x800 else 'cp437', 'replace')
        x, xe = p + 46 + nlen, p + 46 + nlen + xlen
        while x + 4 <= xe:                        # zip64 extra: usize, csize, offset in that order, if saturated
            hid, hl = struct.unpack_from('<HH', cd, x)
            if hid == 1:
                q = x + 4
                if usize == 0xFFFFFFFF: q += 8
                if csize == 0xFFFFFFFF: csize, = struct.unpack_from('<Q', cd, q); q += 8
                if off == 0xFFFFFFFF: off, = struct.unpack_from('<Q', cd, q)
                break
            x += 4 + hl
        out[name] = (method, csize, off)
        p += 46 + nlen + xlen + clen
    return out


async def member_head(get, url, entry, nbytes=65536):
    method, csize, off = entry
    want = 30 + 2048 + min(csize, nbytes)
    buf, _ = await get(url, f'bytes={off}-{off + want - 1}')
    nlen, xlen = struct.unpack_from('<HH', buf, 26)
    data = buf[30 + nlen + xlen:30 + nlen + xlen + min(csize, nbytes)]
    return zlib.decompressobj(-15).decompress(data) if method == 8 else data


def ungz(b, name):
    return zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(b) if name.lower().endswith('.gz') else b


def parse(kind, text):
    """-> list of (software, version)"""
    t = text
    if kind == 'mq_parameters':
        m = re.search(r'^Version\t([^\t\r\n]+)', t, re.M); return [('MaxQuant', m.group(1).strip())] if m else []
    if kind == 'mqpar':
        m = re.search(r'<maxQuantVersion>\s*([^<]+?)\s*<', t); return [('MaxQuant', m.group(1))] if m else []
    if kind == 'diann_log':   # '*report*log.txt' also catches Spectronaut '*_AnalysisLog.txt' files
        m = re.search(r'DIA-NN\s+(\d+(?:\.\d+)+(?:\s*\(?beta\s*\d+\)?)?)', t)
        if m:
            return [('DIA-NN', m.group(1).strip())]
        if 'Initializing Pipeline' in t or 'Spectronaut' in t:
            m = re.search(r'Spectronaut\D{0,20}?(\d+(?:\.\d+)+)', t)
            return [('Spectronaut', m.group(1) if m else '')]
        return []
    if kind == 'fragpipe':
        m = re.search(r'FragPipe\s*\(?v?(\d+(?:\.\d+)+(?:-[\w.]+)?)', t); return [('FragPipe', m.group(1))] if m else []
    if kind == 'skyline':
        m = re.search(r'software_version="([^"]+)"', t)
        if m:
            v = re.search(r'(\d+(?:\.\d+)+)', m.group(1)); return [('Skyline', v.group(1) if v else m.group(1))]
        return []
    if kind == 'mascot_dat':
        if 'x-Mascot' not in t and 'Mascot' not in t[:2000]:
            return []
        m = re.search(r'^version=(\d+(?:\.\d+)+)', t, re.M); return [('Mascot', m.group(1))] if m else []
    if kind == 'mzid':
        out = []
        for a in re.findall(r'<AnalysisSoftware\b([^>]*)>(.*?)</AnalysisSoftware>', t, re.S):
            attrs, body = a
            name = re.search(r'\bname="([^"]*)"', attrs); ver = re.search(r'\bversion="([^"]*)"', attrs)
            cv = re.search(r'<SoftwareName>.*?name="([^"]+)"', body, re.S)
            nm = (cv.group(1) if cv else None) or (name.group(1) if name else None)
            if nm:
                out.append((nm.strip(), ver.group(1).strip() if ver else ''))
        return out
    if kind == 'mztab':
        out = []
        for m in re.finditer(r'^MTD\tsoftware\[\d+\]\t\[([^\]]*)\]', t, re.M):
            parts = [p.strip() for p in m.group(1).split(',')]
            if len(parts) >= 4:
                out.append((parts[2], parts[3]))
            elif len(parts) >= 3:
                out.append((parts[2], ''))
        return out
    return []


async def work(get, row):
    rec = {'url': row.url, 'member': row.member, 'accession': row.accession, 'type': row.type, 'src': row.src}
    try:
        if row.member is None or (isinstance(row.member, float) and np.isnan(row.member)):
            data, _ = await get(row.url, 'bytes=0-65535'); data = ungz(data, row.file_name)
        else:
            cd = await central_directory(get, row.url)
            if row.member not in cd:
                return rec | {'status': 'member_not_found'}
            data = ungz(await member_head(get, row.url, cd[row.member]), row.member)
        sw = parse(row.type, data.decode('utf-8', 'replace'))
        return rec | {'status': 'ok' if sw else 'no_version', 'software': sw}
    except Exception as ex:
        return rec | {'status': f'error: {type(ex).__name__}: {str(ex)[:60]}'}


NORM = [('Mascot', r'mascot(?! ?percolator)'), ('MaxQuant', r'maxquant|andromeda'), ('Proteome Discoverer', r'proteome ?discoverer|\bpd\b'),
        ('Scaffold', r'scaffold'), ('MS-GF+', r'ms-?gf'), ('PeptideShaker', r'peptide ?shaker'), ('SearchGUI', r'searchgui'),
        ('OpenMS', r'openms|topp'), ('quantms', r'quantms'), ('Byonic', r'byonic'), ('PEAKS', r'\bpeaks'), ('ProteinPilot', r'protein ?pilot|paragon'),
        ('X!Tandem', r'tandem'), ('Comet', r'\bcomet'), ('SEQUEST', r'sequest'), ('Percolator', r'percolator'), ('ProteoWizard', r'proteowizard|msconvert'),
        ('Spectrum Mill', r'spectrum ?mill'), ('PIA', r'\bpia\b'), ('MSFragger', r'msfragger|fragpipe'), ('Crux', r'\bcrux|tide'),
        ('MetaMorpheus', r'metamorpheus'), ('Progenesis', r'progenesis'), ('PLGS', r'plgs|proteinlynx'), ('jmzTab / PRIDE converter', r'jmztab|pride ?converter|pride-converter'),
        ('MyriMatch/IDPicker', r'myrimatch|idpicker'), ('Skyline', r'skyline')]


def norm_sw(name, version=''):
    s = f'{name} {version}'.lower()
    for k, rx in NORM:
        if re.search(rx, s):
            return k
    return 'other: ' + name.strip()[:40]


def vkey(v):
    return tuple(int(x) for x in re.findall(r'\d+', v)[:4])


def summary():
    V = [json.loads(l) for l in REC.open()]
    rows = [{'accession': r['accession'], 'type': r['type'], 'src': r['src'], 'software': s, 'version': v}
            for r in V if r.get('status') == 'ok' for s, v in r['software']]
    S = pd.DataFrame(rows)
    Vd = pd.DataFrame(V); Vd['st'] = Vd.status.str.split(':').str[0]
    st = pd.crosstab([Vd.type, Vd.src], Vd.st); st.to_csv(OUT / 'Status-by-type.csv')
    F = pd.read_parquet(O / 'Deposit-features.parquet', columns=['accession', 'submission_date', 'release_year', 'country', 'pi_group', 'dia_expanded'])
    F['submission_date'] = pd.to_datetime(F.submission_date)
    S = S.merge(F, on='accession')
    S.to_csv(OUT / 'Dataset-software-versions.csv.gz', index=False)
    # tools that WROTE the standard files
    W = S[S.type.isin(['mzid', 'mztab'])].copy(); W['product'] = [norm_sw(a, b) for a, b in zip(W.software, W.version)]
    W.to_csv(OUT / 'Standard-file-software.csv.gz', index=False)
    writers = W.groupby(['type', 'product']).accession.nunique().sort_values(ascending=False)
    wy = W.groupby(['product', 'release_year']).accession.nunique().unstack(fill_value=0)
    wy.loc[wy.sum(1).sort_values(ascending=False).index].to_csv(OUT / 'Standard-file-software-by-year.csv')
    writers.to_csv(OUT / 'Standard-file-writers.csv', header=['datasets'])
    # native tool versions: one version per dataset x tool (most common, then highest)
    nat = S[S.type.isin(['mq_parameters', 'mqpar', 'diann_log', 'fragpipe', 'skyline', 'mascot_dat'])].copy()
    nat = nat[nat.version.str.contains(r'\d')]
    nat['vk'] = nat.version.map(vkey); nat['major_minor'] = nat.vk.map(lambda k: '.'.join(map(str, k[:2])))
    nat['release'] = nat.vk.map(lambda k: '.'.join(map(str, k[:3])))     # e.g. MaxQuant 1.6.17, DIA-NN 1.8.1
    one = nat.sort_values('vk').groupby(['accession', 'software']).tail(1)
    # version age: submission date minus the first submission in this archive that used the same release
    # (an internal lower bound on time since release; no external release dates assumed)
    first_seen = one.groupby(['software', 'release']).submission_date.min().rename('first_seen')
    one = one.join(first_seen, on=['software', 'release'])
    one['lag_years'] = (one.submission_date - one.first_seen).dt.days / 365.25
    one.to_csv(OUT / 'Dataset-tool-version.csv.gz', index=False)
    by_year = one.groupby(['software', 'release_year', 'major_minor']).size().rename('datasets').reset_index()
    by_year.to_csv(OUT / 'Version-by-year.csv', index=False)
    lag = one.groupby(['software', 'release_year']).lag_years.agg(['size', 'median', lambda x: x.quantile(.75)]).round(2)
    lag.columns = ['datasets', 'median_lag_years', 'p75_lag_years']; lag.to_csv(OUT / 'Version-lag-by-year.csv')
    rel = one.groupby(['software', 'release_year', 'release']).size().rename('datasets').reset_index()
    rel.to_csv(OUT / 'Release-by-year.csv', index=False)
    top_rel = rel.sort_values('datasets').groupby(['software', 'release_year']).tail(1).set_index(['software', 'release_year'])
    tot = rel.groupby(['software', 'release_year']).datasets.sum()
    top_rel['share_pct'] = (100 * top_rel.datasets / tot).round(1); top_rel.to_csv(OUT / 'Most-used-release-by-year.csv')
    out = {'files_attempted': len(V), 'status_by_type': {f'{k[0]}|{k[1]}': r.to_dict() for k, r in st.iterrows()},
           'datasets_with_native_version': int(one.accession.nunique()),
           'native_versions_by_tool_datasets': one.groupby('software').accession.nunique().to_dict(),
           'median_lag_years_by_tool': one.groupby('software').lag_years.median().round(2).to_dict(),
           'standard_writers_top': {f'{k[0]}|{k[1]}': int(v) for k, v in writers.head(25).items()}}
    (OUT / 'Summary.json').write_text(json.dumps(out, indent=2, default=str))
    pd.set_option('display.width', 250)
    print(json.dumps(out, indent=1, default=str)); print(lag.to_string())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--summary', action='store_true'); ap.add_argument('--workers', type=int, default=128)
    a = ap.parse_args()
    if not a.summary:
        t = targets(); t.to_csv(OUT / 'Targets.csv.gz', index=False)
        done = set()
        if REC.exists():
            done = {(r['url'], r['member']) for r in map(json.loads, REC.open())}
        todo = [r for r in t.itertuples(index=False) if (r.url, r.member if isinstance(r.member, str) else None) not in done]
        print(len(t), 'targets;', t.groupby('type').size().to_dict(), ';', len(todo), 'to fetch', flush=True)
        run_jsonl(todo, work, REC, concurrency=a.workers)
    summary()


if __name__ == '__main__':
    main()
