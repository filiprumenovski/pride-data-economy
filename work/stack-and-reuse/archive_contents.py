"""List what is inside deposited results archives without downloading them.

A zip file ends with a central directory: the name and size of every member. Two HTTP range
requests (the last 64 KiB, then the directory itself if it is not in that tail) recover the
full member list of an archive of any size. Only the directory bytes are transferred.

Targets: every non-generated .zip in a study dataset classified by layers.py as tool_results,
archived_results or other (vendor-raw zips such as *.raw.zip / *.d.zip are not results and are
skipped). tar.gz, tgz, rar and 7z cannot be listed this way (compressed or no trailing index)
and are recorded as unreadable, so the denominator stays explicit.

Resumable: listings are appended to Archive-contents/Listings.jsonl, one line per archive, and
archives already listed are skipped on rerun. Requests are latency-bound, so they run
concurrently (fastfetch.py: asyncio, 128 in flight by default), each a small range read, with
retries and a User-Agent naming the study.

  python archive_contents.py               # all target zips
  python archive_contents.py --datasets 500 # random sample of datasets (fixed seed)
  python archive_contents.py --probe        # check the host answers range requests, then stop
"""
import argparse, json, struct
import pandas as pd
from common import R, O, SEED
from layer_rules import RULES
from fastfetch import run_jsonl
import sys; sys.path.insert(0, str(R / 'scripts'))
from query import connect

OUT = O / 'Archive-contents'; OUT.mkdir(exist_ok=True)
LISTINGS = OUT / 'Listings.jsonl'
HOST = 'https://ftp.pride.ebi.ac.uk/'
UA = 'pride-archive-anatomy-study/1.0 (research; zip central-directory listing only)'
TAIL = 65536 + 22
MAX_CD = 64 * 2**20          # do not fetch a central directory larger than 64 MiB
TARGET_LAYERS = ('tool_results', 'archived_results', 'other')


def https_url(loc):
    """ftp://ftp.pride.ebi.ac.uk/pride/... or prd_ascp@fasp.ebi.ac.uk:pride/... -> https URL."""
    if loc.startswith('ftp://ftp.pride.ebi.ac.uk/'):
        return HOST + loc[len('ftp://ftp.pride.ebi.ac.uk/'):]
    if ':pride/' in loc:
        return HOST + loc.split(':', 1)[1]
    return None


def targets():
    case = 'case ' + ' '.join(f"when {cond} then '{layer}'" for layer, cond in RULES) + " else 'other' end"
    c = connect()
    t = c.execute(f"""
    with f as (select project_accession a, file_name, lower(file_name) n,
      lower(regexp_replace(file_name, '\\.(gz|zip|bz2|7z|tar|rar|xz|tgz)$', '', 'i')) b,
      coalesce(file_category_value, '') cat, file_size_bytes, primary_public_location_value ploc
      from pride_files where coalesce(primary_public_location_value, '') not like '%/generated/%'
      and project_accession in (select accession from study))
    select a accession, file_name, cat, file_size_bytes, ploc, {case} layer,
           regexp_extract(n, '\\.(zip|tgz|tar\\.gz|tar|rar|7z)$', 1) fmt
    from f where regexp_matches(n, '\\.(zip|tgz|tar|rar|7z)$|\\.tar\\.gz$')""").df()
    t = t[t.layer.isin(TARGET_LAYERS)].copy()
    t['url'] = t.ploc.map(lambda x: https_url(x or ''))
    return t


def parse_cd(buf, n_expected):
    """Parse central-directory records -> list of (name, uncompressed_size)."""
    out, p = [], 0
    while p + 46 <= len(buf) and buf[p:p + 4] == b'PK\x01\x02':
        flags, = struct.unpack_from('<H', buf, p + 8)
        usize, = struct.unpack_from('<I', buf, p + 24)
        nlen, xlen, clen = struct.unpack_from('<HHH', buf, p + 28)
        raw = buf[p + 46:p + 46 + nlen]
        name = raw.decode('utf-8' if flags & 0x800 else 'cp437', 'replace')
        if usize == 0xFFFFFFFF:                       # zip64 extra field carries the real size
            x, e = p + 46 + nlen, p + 46 + nlen + xlen
            while x + 4 <= e:
                hid, hlen = struct.unpack_from('<HH', buf, x)
                if hid == 1 and hlen >= 8:
                    usize, = struct.unpack_from('<Q', buf, x + 4); break
                x += 4 + hlen
        out.append((name, int(usize)))
        p += 46 + nlen + xlen + clen
    return out, len(out) == n_expected


async def list_zip(get, url):
    tail, total = await get(url, f'bytes=-{TAIL}')
    e = tail.rfind(b'PK\x05\x06')
    if e < 0:
        return {'status': 'no_zip_directory'}
    n, cd_size, cd_off = struct.unpack_from('<HII', tail, e + 10)
    if 0xFFFFFFFF in (cd_size, cd_off) or n == 0xFFFF:
        loc = tail.rfind(b'PK\x06\x07', 0, e)
        if loc < 0:
            return {'status': 'zip64_locator_missing'}
        z64_off, = struct.unpack_from('<Q', tail, loc + 8)
        tail_start = total - len(tail)
        if z64_off >= tail_start:
            rec = tail[z64_off - tail_start:z64_off - tail_start + 56]
        else:
            rec, _ = await get(url, f'bytes={z64_off}-{z64_off + 55}')
        if rec[:4] != b'PK\x06\x06':
            return {'status': 'zip64_record_bad'}
        n, cd_size, cd_off = struct.unpack_from('<QQQ', rec, 32)
    if cd_size > MAX_CD:
        return {'status': 'directory_too_large', 'n_members': int(n)}
    tail_start = total - len(tail)
    if cd_off >= tail_start:
        cd = tail[cd_off - tail_start:cd_off - tail_start + cd_size]
    else:
        cd, _ = await get(url, f'bytes={cd_off}-{cd_off + cd_size - 1}')
    members, complete = parse_cd(cd, n)
    return {'status': 'ok' if complete else 'partial', 'n_members': int(n), 'members': members}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--datasets', type=int, default=0, help='random sample of N datasets (0 = all)')
    ap.add_argument('--workers', type=int, default=128)
    ap.add_argument('--probe', action='store_true')
    a = ap.parse_args()
    t = targets()
    t.to_csv(OUT / 'Targets.csv.gz', index=False)
    summary = t.groupby(['layer', 'fmt']).agg(files=('url', 'size'), datasets=('accession', 'nunique')).reset_index()
    summary.to_csv(OUT / 'Target-denominators.csv', index=False)
    z = t[(t.fmt == 'zip') & t.url.notna() & (t.file_size_bytes > 0)]   # 4 empty zips skipped
    if a.datasets:
        keep = pd.Series(sorted(z.accession.unique())).sample(a.datasets, random_state=SEED)
        z = z[z.accession.isin(keep)]
    done = set()
    if LISTINGS.exists():
        done = {l[9:l.index('"', 9)] for l in LISTINGS.open() if l.startswith('{"url": "')}   # url is the first key
    todo = z[~z.url.isin(done)]
    print(f'{len(t)} candidate archives ({t.accession.nunique()} datasets); {len(z)} zips selected; '
          f'{len(done)} already listed; {len(todo)} to fetch')
    if a.probe:
        row = z.sort_values('file_size_bytes').iloc[0]
        print(row.url)
        run_jsonl([row], work, OUT / 'Probe.jsonl', concurrency=1); return
    run_jsonl(todo.itertuples(index=False), work, LISTINGS, concurrency=a.workers)


async def work(get, row):
    try:
        res = await list_zip(get, row.url)
    except Exception as ex:
        res = {'status': str(ex)}
    return dict(url=row.url, accession=row.accession, file_name=row.file_name, layer=row.layer,
                cat=row.cat, size=int(row.file_size_bytes or 0), **res)


if __name__ == '__main__':
    main()
