"""Concurrent HTTP range reads for the archive and header scans.

The scans are latency-bound (0.4-1 s to first byte per request from EBI, microseconds of local
work), so throughput comes from keeping many small requests in flight: asyncio + aiohttp with a
shared connection pool, through the environment's proxy (trust_env) and CA bundle.

  get(url, 'bytes=-65558') -> (bytes, total_size)   # raises RuntimeError on 404/416/no-range
  run_jsonl(rows, work, path, concurrency)           # work(get, row) -> dict, appended as JSON lines
"""
import asyncio, json, os, ssl, time
import aiohttp

UA = 'pride-archive-anatomy-study/1.0 (research; ranged reads of file indexes and table headers only)'


class RangeClient:
    def __init__(self, session):
        self.s = session

    async def get(self, url, rng):
        for i in range(5):
            try:
                async with self.s.get(url, headers={'Range': rng}) as r:
                    if r.status == 206:
                        total = int(r.headers.get('Content-Range', '/0').rsplit('/', 1)[1] or 0)
                        return await r.read(), total
                    if r.status == 200:            # range ignored: never download a whole file
                        raise RuntimeError('range_not_supported')
                    if r.status in (404, 410, 416):
                        raise RuntimeError(f'http_{r.status}')
            except (aiohttp.ClientError, asyncio.TimeoutError):
                pass
            await asyncio.sleep(2 ** i)
        raise RuntimeError('http_failed')


async def _run(rows, work, path, concurrency, every):
    ctx = ssl.create_default_context(cafile=os.environ.get('SSL_CERT_FILE') or None)
    conn = aiohttp.TCPConnector(limit=concurrency, limit_per_host=concurrency, ssl=ctx, ttl_dns_cache=3600)
    timeout = aiohttp.ClientTimeout(total=120, sock_read=60)
    queue = asyncio.Queue(maxsize=concurrency * 4)
    done, t0 = 0, time.time()
    async with aiohttp.ClientSession(connector=conn, timeout=timeout, trust_env=True, headers={'User-Agent': UA}) as s:
        client = RangeClient(s)
        with open(path, 'a') as fh:
            async def worker():
                nonlocal done
                while True:
                    row = await queue.get()
                    if row is None:
                        return
                    try:
                        rec = await work(client.get, row)
                    except Exception as ex:   # a failure is recorded, never fatal for the run
                        rec = {'url': getattr(row, 'url', None), 'status': f'error: {type(ex).__name__}: {str(ex)[:60]}'}
                    fh.write(json.dumps(rec) + '\n'); done += 1
                    if done % every == 0:
                        fh.flush(); r = done / (time.time() - t0)
                        print(f'{done} of {len(rows)}  {r:.0f}/s  eta {(len(rows) - done) / r / 60:.1f} min', flush=True)
            ws = [asyncio.create_task(worker()) for _ in range(concurrency)]
            for row in rows:
                await queue.put(row)
            for _ in ws:
                await queue.put(None)
            await asyncio.gather(*ws)
    print(f'{done} done in {(time.time() - t0) / 60:.1f} min', flush=True)


def run_jsonl(rows, work, path, concurrency=96, every=2000):
    asyncio.run(_run(list(rows), work, path, concurrency, every))
