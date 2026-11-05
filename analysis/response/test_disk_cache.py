"""Cold/warm disk reuse must preserve every sampled bit and RNG state."""
import json
from pathlib import Path

import numpy as np
import pytest

from analysis.response.cache import CompactMatchedKernel, _ExactValues
from analysis.response.disk_cache import DiskMatchedKernel, pack_sampler, unpack_sampler

ROOT = Path(__file__).resolve().parents[2]


def kernel(directory, **kwargs):
    p = ROOT/'analysis/matched/protocol.json'
    specification = json.loads(p.read_text())
    return DiskMatchedKernel(ROOT/specification['host_path'],
        ROOT/'results/matched-residence/kernel-fine.npz', p,
        disk_cache_dir=directory, namespace='ab'*32, grid='fine',
        cache_bytes=64*1024**2, cache_vertices=2, **kwargs)


def assert_same_draws(first, second):
    left, right = np.random.default_rng(30219), np.random.default_rng(30219)
    for count in (0, 1, 64, 4096):
        np.testing.assert_array_equal(first.sample(count, left), second.sample(count, right))
        assert left.bit_generator.state == right.bit_generator.state
    for name in ('node_count', 'remnant', 'G', 'kind', 'velocity_scale'):
        assert type(getattr(first, name)) is type(getattr(second, name))
    assert first.nbytes == second.nbytes


def test_actual_vertices_cold_warm_and_eviction_are_bitwise_identical(tmp_path):
    cold = kernel(tmp_path)
    coordinates = [(0, 0, 0), (2, 3, 0), (6, 4, 4), (8, 2, 8)]
    indices = [int(np.ravel_multi_index(point, cold.shape)) for point in coordinates]
    for index in indices:
        fresh = CompactMatchedKernel._physical_sampler(cold, index)
        stored = cold._sampler(index)
        assert_same_draws(fresh, stored)
    assert cold.disk_statistics()['builds'] == 4
    assert cold.cache_statistics()['evictions'] == 2
    warm = kernel(tmp_path)
    for index in indices:
        fresh = CompactMatchedKernel._physical_sampler(warm, index)
        restored = warm._sampler(index)
        assert_same_draws(fresh, restored)
    assert warm.disk_statistics()['loads'] == 4
    assert warm.disk_statistics()['builds'] == 0
    assert warm.cache_statistics()['resident_bytes'] <= 64*1024**2


@pytest.mark.parametrize('values', [
    np.full(100, .125), np.tile([.123456789012345, .3], 50), np.arange(100, dtype=float),
])
def test_bound_dictionary_encoding_and_numpy_scalar_types_roundtrip(tmp_path, values):
    instance = kernel(tmp_path)
    original = CompactMatchedKernel._physical_sampler(instance, 0)
    # Keep a structurally valid sampler with each dictionary representation.
    count = len(original.cdf)
    values = np.resize(values, count)
    original.lower = _ExactValues(values)
    original.upper = _ExactValues(values+1.)
    original.G = np.float64(original.G)
    arrays, metadata = pack_sampler(original)
    restored = unpack_sampler(arrays, metadata)
    for name in ('lower', 'upper'):
        before, after = getattr(original, name), getattr(restored, name)
        assert before.constant == after.constant
        assert before.size == after.size
        np.testing.assert_array_equal(before.values, after.values)
        assert (before.indices is None) == (after.indices is None)
        if before.indices is not None:
            assert before.indices.dtype == after.indices.dtype
            np.testing.assert_array_equal(before.indices, after.indices)
    assert type(restored.G) is np.float64


def test_corrupt_npz_is_rejected_without_rebuilding(tmp_path):
    cold = kernel(tmp_path)
    cold._sampler(0)
    path, _ = cold._entry_paths(0)
    with path.open('ab') as handle:
        handle.write(b'corrupt')
    warm = kernel(tmp_path)
    with pytest.raises(ValueError, match='checksum'):
        warm._sampler(0)
    assert warm.disk_statistics()['builds'] == 0


def test_missing_manifest_is_rejected_without_rebuilding(tmp_path):
    cold = kernel(tmp_path)
    cold._sampler(0)
    _, metadata = cold._entry_paths(0)
    metadata.unlink()
    warm = kernel(tmp_path)
    with pytest.raises(ValueError, match='Incomplete'):
        warm._sampler(0)
    assert warm.disk_statistics()['builds'] == 0


def test_disk_budget_refuses_oversized_write(tmp_path):
    instance = kernel(tmp_path, disk_max_bytes=10000)
    with pytest.raises(MemoryError, match='disk budget'):
        instance._sampler(0)
    assert not list(instance.disk_path.glob('vertex-*'))


def test_identity_key_separates_different_run_namespaces(tmp_path):
    first = kernel(tmp_path)
    p = ROOT/'analysis/matched/protocol.json'
    specification = json.loads(p.read_text())
    second = DiskMatchedKernel(ROOT/specification['host_path'],
        ROOT/'results/matched-residence/kernel-fine.npz', p,
        disk_cache_dir=tmp_path, namespace='cd'*32, grid='fine')
    assert first.disk_path != second.disk_path
