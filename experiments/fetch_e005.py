"""Read locked PARAM members from a pinned Git archive, without extracting paths."""
from pathlib import Path, PurePosixPath
import hashlib
import json
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def verified_bytes(data, entry):
    """Verify the exact Git object, including its blob header, and byte size."""
    if len(data) != entry['size_bytes']:
        raise ValueError(f"Size mismatch: {entry['path']}")
    digest = hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if digest != entry['git_blob_sha1']:
        raise ValueError(f"Git blob mismatch: {entry['path']}")
    return hashlib.sha256(data).hexdigest()


def fetch(config):
    lock = json.loads((ROOT/config['input_lock']).read_text())
    for key in ['source_repository', 'source_revision']:
        if lock[key] != config[key]:
            raise ValueError(f'Input lock mismatch: {key}')
    if lock['tree_truncated'] is not False:
        raise ValueError('Input lock must come from a complete Git tree')
    entries = {item['path']: item for item in lock['files']}
    if len(entries) != len(lock['files']) or len(entries) != config['expected_snapshot_count']:
        raise ValueError('Duplicate or unexpected input paths')
    for name in entries:
        parts = PurePosixPath(name).parts
        if len(parts) != 3 or parts[0] != 'CMC' or any(p in ('.', '..') for p in parts):
            raise ValueError('Unexpected catalog member path')
    if len({name.split('/')[1] for name in entries}) != config['expected_model_count']:
        raise ValueError('Unexpected number of model histories')
    cache = ROOT/config['cache_directory']
    cache.mkdir(parents=True, exist_ok=True)
    paths = {name: cache/name for name in entries}
    missing = []
    for name, path in paths.items():
        if path.exists():
            verified_bytes(path.read_bytes(), entries[name])
        else:
            missing.append(name)
    if missing:
        archive = cache/'upstream.tar.gz'
        if not archive.exists():
            url = ('https://codeload.github.com/'+config['source_repository']+
                   '/tar.gz/'+config['source_revision'])
            partial = archive.with_suffix('.part')
            try:
                # curl uses the OS certificate store; TLS verification stays enabled.
                subprocess.run(['curl', '--fail', '--location', '--silent', '--show-error',
                                '--retry', '2', '--max-time', '180', '--max-filesize',
                                str(config['maximum_archive_bytes']), url, '--output', str(partial)],
                               check=True)
                if partial.stat().st_size > config['maximum_archive_bytes']:
                    raise ValueError('Archive exceeds declared size limit')
                partial.replace(archive)
            finally:
                partial.unlink(missing_ok=True)
        if archive.stat().st_size > config['maximum_archive_bytes']:
            raise ValueError('Cached archive exceeds declared size limit')
        found = set()
        wanted = set(missing)
        with tarfile.open(archive, 'r:gz') as tar:
            for member in tar:
                parts = PurePosixPath(member.name).parts
                name = '/'.join(parts[1:])
                if name not in wanted:
                    continue
                if name in found or not member.isfile() or member.size != entries[name]['size_bytes']:
                    raise ValueError(f'Unexpected archive member: {name}')
                # Never execute source files or extract arbitrary archive paths.
                with tar.extractfile(member) as stream:
                    data = stream.read()
                verified_bytes(data, entries[name])
                paths[name].parent.mkdir(parents=True, exist_ok=True)
                paths[name].write_bytes(data)
                found.add(name)
        if found != wanted:
            raise ValueError(f'Archive missing {len(wanted-found)} locked inputs')
    return paths


if __name__ == '__main__':
    config = json.loads((ROOT/'experiments/e005.json').read_text())
    print(f'Verified {len(fetch(config))} pinned PARAM files; no population analysis performed.')
