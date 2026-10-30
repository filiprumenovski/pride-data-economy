"""Verify the release inventory and public analytical input manifests offline."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]


def verify(manifest, base):
    entries = json.loads(manifest.read_text())['files']
    for entry in entries:
        path = (base / entry['path']).resolve()
        if not path.is_relative_to(base.resolve()):
            raise ValueError(f'Path outside release: {entry["path"]}')
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != entry['sha256']:
            raise ValueError(f'Checksum mismatch: {entry["path"]}')
    print(f'{manifest.relative_to(ROOT)}: {len(entries)} checksums verified')


if __name__ == '__main__':
    verify(ROOT / 'MANIFEST.json', ROOT)
    verify(ROOT / 'data/early/Data-manifest.json', ROOT / 'data/early')
    verify(ROOT / 'data/late/Input-manifest.json', ROOT / 'data/late')
