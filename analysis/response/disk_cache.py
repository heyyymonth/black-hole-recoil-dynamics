"""Lossless disk backing for the unchanged compact encounter sampler.

Only numerical storage changes. Float64 CDFs, packed indices, all scalar types,
and exact bound dictionaries are restored without fitting or quantization.
Entries are uncompressed NPZ files loaded with allow_pickle=False, paired with
SHA256 manifests. An incomplete or corrupt entry is rejected, never rebuilt
silently. The cache directory is disposable local acceleration, not evidence.
One writing process per cache directory is supported.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile
import time

import numpy as np

from analysis.response.cache import CompactMatchedKernel, CompactVertexSampler, _ExactValues

ARRAY_FIELDS = ('cdf', 'positive', 'relative_nodes', 'mass_values', 'rotation')
SCALAR_FIELDS = ('node_count', 'remnant', 'G', 'velocity_scale', 'kind')
FORMAT = 'exact-compact-sampler-npz-1'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()


def pack_sampler(sampler):
    arrays = {key: np.asarray(getattr(sampler, key)) for key in ARRAY_FIELDS}
    scalars = {}
    for key in SCALAR_FIELDS:
        value = getattr(sampler, key)
        array = np.asarray(value)
        if array.ndim != 0 or array.dtype.hasobject:
            raise ValueError('Sampler scalar cannot be serialized safely: '+key)
        arrays[key] = array
        scalars[key] = ('numpy' if isinstance(value, np.generic) else type(value).__name__)
    bounds = {}
    for key in ('lower', 'upper'):
        value = getattr(sampler, key)
        arrays[key+'_values'] = value.values
        if value.indices is not None:
            arrays[key+'_indices'] = value.indices
        bounds[key] = {'size': int(value.size), 'constant': bool(value.constant),
                       'has_indices': value.indices is not None}
    if any(value.dtype.hasobject for value in arrays.values()):
        raise ValueError('Object arrays are forbidden in the disk cache')
    return arrays, {'scalar_types': scalars, 'bounds': bounds,
                    'array_inventory': {key: {'shape': list(value.shape), 'dtype': value.dtype.str}
                                        for key, value in arrays.items()},
                    'resident_nbytes': int(sampler.nbytes)}


def unpack_sampler(arrays, metadata):
    inventory = metadata['array_inventory']
    if set(arrays) != set(inventory):
        raise ValueError('Disk sampler array inventory mismatch')
    for key, value in arrays.items():
        if (value.dtype.hasobject or list(value.shape) != inventory[key]['shape']
                or value.dtype.str != inventory[key]['dtype']):
            raise ValueError('Disk sampler array shape/type mismatch: '+key)
    sampler = CompactVertexSampler.__new__(CompactVertexSampler)
    for key in ARRAY_FIELDS:
        setattr(sampler, key, arrays[key])
    for key in SCALAR_FIELDS:
        value = arrays[key]
        if value.ndim != 0:
            raise ValueError('Disk sampler scalar is not scalar')
        kind = metadata['scalar_types'][key]
        if kind == 'numpy':
            scalar = value[()]
        elif kind in ('int', 'float', 'str', 'bool'):
            scalar = value.item()
            if type(scalar).__name__ != kind:
                raise ValueError('Disk sampler scalar type mismatch')
        else:
            raise ValueError('Unsupported disk sampler scalar type')
        setattr(sampler, key, scalar)
    for key in ('lower', 'upper'):
        info = metadata['bounds'][key]
        bound = _ExactValues.__new__(_ExactValues)
        bound.size = info['size']
        bound.constant = info['constant']
        bound.values = arrays[key+'_values']
        bound.indices = arrays[key+'_indices'] if info['has_indices'] else None
        if (type(bound.size) is not int or bound.size < 1 or type(bound.constant) is not bool
                or bound.values.ndim != 1 or bound.values.dtype != np.dtype('float64')
                or not np.isfinite(bound.values).all()):
            raise ValueError('Invalid disk bound dictionary')
        if bound.constant and (len(bound.values) != 1 or bound.indices is not None):
            raise ValueError('Invalid constant disk bound')
        if bound.indices is not None and (bound.indices.shape != (bound.size,)
                or bound.indices.dtype.kind != 'u' or np.any(bound.indices >= len(bound.values))):
            raise ValueError('Invalid disk bound indices')
        if not bound.constant and bound.indices is None and len(bound.values) != bound.size:
            raise ValueError('Invalid unpacked disk bound')
        setattr(sampler, key, bound)
    count = len(sampler.cdf)
    if (sampler.cdf.dtype != np.dtype('float64') or sampler.cdf.shape != (count,)
            or count < 1 or not np.isfinite(sampler.cdf).all() or np.any(np.diff(sampler.cdf) < 0)
            or sampler.cdf[0] <= 0 or sampler.cdf[-1] != 1.
            or sampler.positive.shape != (count,) or sampler.positive.dtype.kind != 'u'
            or sampler.lower.size != count or sampler.upper.size != count
            or sampler.relative_nodes.shape != (sampler.node_count, 3)
            or sampler.rotation.shape != (3, 3) or sampler.mass_values.ndim != 1
            or np.any(sampler.positive >= sampler.node_count*len(sampler.mass_values))):
        raise ValueError('Invalid restored sampler structure')
    if sampler.nbytes != metadata['resident_nbytes']:
        raise ValueError('Restored sampler storage differs from manifest')
    return sampler


class DiskMatchedKernel(CompactMatchedKernel):
    """Unchanged in-memory policy with checksum-verified persistent samplers."""
    def __init__(self, host_path, table_path, protocol_path, *, disk_cache_dir,
                 namespace, grid, cache_bytes=2*1024**3, cache_vertices=768,
                 disk_max_bytes=15*1024**3):
        if (not isinstance(namespace, str) or len(namespace) != 64
                or any(character not in '0123456789abcdef' for character in namespace)):
            raise ValueError('namespace must be a lowercase SHA256 of the declared run identity')
        if grid not in ('coarse', 'fine'):
            raise ValueError('grid must be coarse or fine')
        if type(disk_max_bytes) is not int or disk_max_bytes < 1:
            raise ValueError('disk_max_bytes must be a positive integer')
        super().__init__(host_path, table_path, protocol_path,
                         cache_bytes=cache_bytes, cache_vertices=cache_vertices)
        self.disk_root = Path(disk_cache_dir)
        self.disk_root.mkdir(parents=True, exist_ok=True)
        self.disk_max_bytes = disk_max_bytes
        self.disk_identity = {'format': FORMAT, 'namespace': namespace, 'grid': grid,
            'adapter_sha256': sha256(__file__),
            'compact_cache_sha256': sha256(Path(__file__).with_name('cache.py')),
            'host_sha256': sha256(host_path), 'table_sha256': sha256(table_path),
            'table_record_sha256': sha256(Path(table_path).with_suffix('.json')),
            'matched_protocol_sha256': sha256(protocol_path)}
        key = hashlib.sha256(_json_bytes(self.disk_identity)).hexdigest()
        self.disk_path = self.disk_root/f'{grid}-{key}'
        self.disk_path.mkdir(exist_ok=True)
        manifest = self.disk_path/'identity.json'
        if manifest.exists():
            if json.loads(manifest.read_text()) != self.disk_identity:
                raise ValueError('Disk namespace identity mismatch')
        else:
            if any(self.disk_path.iterdir()):
                raise ValueError('Disk namespace has entries without an identity manifest')
            self._atomic_bytes(manifest, _json_bytes(self.disk_identity))
        self._disk_bytes = sum(path.stat().st_size for path in self.disk_root.rglob('*') if path.is_file())
        if self._disk_bytes > self.disk_max_bytes:
            raise MemoryError('Existing disk cache exceeds its declared byte budget')
        self._disk_metrics = {'loads': 0, 'builds': 0, 'read_seconds': 0.,
                              'write_seconds': 0., 'build_seconds': 0.,
                              'read_bytes': 0, 'written_bytes': 0}

    @staticmethod
    def _atomic_bytes(path, data):
        descriptor, temporary = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=path.parent)
        try:
            with os.fdopen(descriptor, 'wb') as handle:
                handle.write(data)
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    def _entry_paths(self, index):
        path = self.disk_path/f'vertex-{int(index):05d}.npz'
        return path, path.with_suffix('.json')

    def _read_sampler(self, path, manifest, index):
        started = time.perf_counter()
        record = json.loads(manifest.read_text())
        if (record.get('identity') != self.disk_identity or record.get('vertex') != int(index)
                or record.get('npz_sha256') != sha256(path)):
            raise ValueError('Disk sampler checksum/identity mismatch: '+str(path))
        with np.load(path, allow_pickle=False) as loaded:
            arrays = {key: loaded[key] for key in loaded.files}
        sampler = unpack_sampler(arrays, record['sampler'])
        self._disk_metrics['loads'] += 1
        self._disk_metrics['read_bytes'] += path.stat().st_size+manifest.stat().st_size
        self._disk_metrics['read_seconds'] += time.perf_counter()-started
        return sampler

    def _write_sampler(self, path, manifest, index, sampler):
        started = time.perf_counter()
        arrays, metadata = pack_sampler(sampler)
        descriptor, temporary = tempfile.mkstemp(prefix=path.name+'.', suffix='.tmp', dir=path.parent)
        try:
            with os.fdopen(descriptor, 'wb') as handle:
                np.savez(handle, **arrays)
            record = {'identity': self.disk_identity, 'vertex': int(index),
                      'npz_sha256': sha256(temporary), 'sampler': metadata}
            payload = _json_bytes(record)
            size = os.path.getsize(temporary)+len(payload)
            if self._disk_bytes+size > self.disk_max_bytes:
                raise MemoryError('Writing this sampler would exceed the declared disk budget')
            if path.exists() or manifest.exists():
                raise FileExistsError('Concurrent disk cache write; only one writer is supported')
            os.replace(temporary, path)
            self._atomic_bytes(manifest, payload)
            self._disk_bytes += size
            self._disk_metrics['written_bytes'] += size
            self._disk_metrics['write_seconds'] += time.perf_counter()-started
        finally:
            if os.path.exists(temporary):
                os.unlink(temporary)

    def _physical_sampler(self, index):
        path, manifest = self._entry_paths(index)
        if path.exists() or manifest.exists():
            if not (path.exists() and manifest.exists()):
                raise ValueError('Incomplete disk sampler entry; explicit cache repair required: '+str(path))
            return self._read_sampler(path, manifest, index)
        started = time.perf_counter()
        sampler = super()._physical_sampler(index)
        self._disk_metrics['build_seconds'] += time.perf_counter()-started
        self._disk_metrics['builds'] += 1
        self._write_sampler(path, manifest, index, sampler)
        return sampler

    def disk_statistics(self):
        return {**self._disk_metrics, 'cache_path': str(self.disk_path),
                'bytes_in_cache_root_at_last_write': self._disk_bytes,
                'maximum_disk_bytes': self.disk_max_bytes, 'identity': self.disk_identity,
                'scope': 'Lossless storage only; disk budget covers this single-writer cache root; in-memory policy is inherited unchanged'}
