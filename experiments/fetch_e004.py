"""Fetch the pinned external snapshot into ignored local cache, verify SHA-256."""
from pathlib import Path
import hashlib
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[1]


def fetch(config):
    cache = ROOT/config['cache_directory']
    cache.mkdir(parents=True, exist_ok=True)
    paths = {}
    for key, source in config['files'].items():
        path = cache/Path(source['path']).name
        if not path.exists():
            url = ('https://raw.githubusercontent.com/'+config['source_repository']
                   +'/'+config['source_revision']+'/'+source['path'])
            temporary = path.with_suffix(path.suffix+'.part')
            try:
                with urllib.request.urlopen(url, timeout=60) as response, temporary.open('wb') as output:
                    total = 0
                    while chunk := response.read(1024*1024):
                        total += len(chunk)
                        if total > 100_000_000:
                            raise ValueError('External input exceeds declared 100 MB fetch limit')
                        output.write(chunk)
                if hashlib.sha256(temporary.read_bytes()).hexdigest() != source['sha256']:
                    raise ValueError(f'Checksum mismatch for downloaded {key}')
                temporary.replace(path)
            finally:
                temporary.unlink(missing_ok=True)
        if hashlib.sha256(path.read_bytes()).hexdigest() != source['sha256']:
            raise ValueError(f'Checksum mismatch for cached {key}: {path}')
        paths[key] = path
    return paths


if __name__ == '__main__':
    cfg = json.loads((ROOT/'experiments/e004.json').read_text())
    for key, path in fetch(cfg).items():
        print(f'Verified external {key}: {path.relative_to(ROOT)}')
