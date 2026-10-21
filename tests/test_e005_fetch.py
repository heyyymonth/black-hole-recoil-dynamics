"""Small archive fixtures test provenance and selective, non-executing reads."""
from pathlib import Path
import hashlib
import io
import json
import sys
import tarfile

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'experiments'))
import fetch_e005


def fixture_archive(tmp_path, monkeypatch):
    monkeypatch.setattr(fetch_e005, 'ROOT', tmp_path)
    data = b'age 12\n'
    entry = {'path': 'CMC/model/PARAM0001.txt', 'size_bytes': len(data),
             'git_blob_sha1': hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()}
    lock = {'source_repository': 'owner/repo', 'source_revision': 'pinned',
            'tree_truncated': False, 'files': [entry]}
    (tmp_path/'inputs.json').write_text(json.dumps(lock))
    cache = tmp_path/'cache'
    cache.mkdir()
    with tarfile.open(cache/'upstream.tar.gz', 'w:gz') as tar:
        for name, content in [('repo-pinned/'+entry['path'], data),
                              ('repo-pinned/../../must-not-extract', b'untrusted')]:
            info = tarfile.TarInfo(name)
            info.size = len(content)
            tar.addfile(info, io.BytesIO(content))
    config = {**{k: lock[k] for k in ['source_repository', 'source_revision']},
              'input_lock': 'inputs.json', 'cache_directory': 'cache',
              'expected_snapshot_count': 1, 'expected_model_count': 1,
              'maximum_archive_bytes': 100000}
    return config, entry, data


def test_reads_only_locked_verified_member(tmp_path, monkeypatch):
    config, entry, data = fixture_archive(tmp_path, monkeypatch)
    result = fetch_e005.fetch(config)
    assert result[entry['path']].read_bytes() == data
    assert not (tmp_path/'must-not-extract').exists()
    (tmp_path/'cache/upstream.tar.gz').unlink()
    assert fetch_e005.fetch(config) == result  # subsequent fetch is offline


def test_changed_cached_bytes_are_rejected(tmp_path, monkeypatch):
    config, entry, data = fixture_archive(tmp_path, monkeypatch)
    result = fetch_e005.fetch(config)
    result[entry['path']].write_bytes(data.replace(b'12', b'13'))
    with pytest.raises(ValueError, match='Git blob mismatch'):
        fetch_e005.fetch(config)


def test_incomplete_tree_is_rejected(tmp_path, monkeypatch):
    config, _, _ = fixture_archive(tmp_path, monkeypatch)
    path = tmp_path/'inputs.json'
    lock = json.loads(path.read_text())
    lock['tree_truncated'] = True
    path.write_text(json.dumps(lock))
    with pytest.raises(ValueError, match='complete Git tree'):
        fetch_e005.fetch(config)


def test_mismatched_revision_is_rejected(tmp_path, monkeypatch):
    config, _, _ = fixture_archive(tmp_path, monkeypatch)
    config['source_revision'] = 'another'
    with pytest.raises(ValueError, match='lock mismatch'):
        fetch_e005.fetch(config)


def test_git_blob_includes_header():
    data = b'age 12\n'
    entry = {'path': 'x', 'size_bytes': len(data), 'git_blob_sha1': hashlib.sha1(data).hexdigest()}
    with pytest.raises(ValueError, match='Git blob mismatch'):
        fetch_e005.verified_bytes(data, entry)
