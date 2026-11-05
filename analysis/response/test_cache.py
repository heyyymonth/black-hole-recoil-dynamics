"""Exact-law and bounded-storage checks for the independent cache adapter."""
from collections import OrderedDict
import json

import numpy as np
import pytest

from analysis.matched.kernel import (ROOT, Host, MatchedKernel, VertexSampler,
                                    host_scales, vertex_law)
from analysis.response.cache import (CompactMatchedKernel, CompactVertexSampler,
                                     _ExactValues)


@pytest.fixture(scope='module')
def physical_inputs():
    p = json.loads((ROOT/'analysis/matched/protocol.json').read_text())
    e = json.loads((ROOT/p['encounter_protocol']).read_text())
    h = Host(dict(np.load(ROOT/p['host_path'], allow_pickle=False)))
    return h, host_scales(h), p, e


@pytest.mark.parametrize('coordinates', [
    (0., 0., 0.), (0., .83, .61), (.93, 0., 0.),
    (.37, 1.41, 0.), (.63, 1.17, .31), (1.83, .37, 1.),
])
def test_compact_vertex_draws_and_rng_are_bitwise_unchanged(physical_inputs, coordinates):
    h, scales, p, e = physical_inputs
    _, law = vertex_law(h, coordinates, scales, p, e, 28)
    ordinary = VertexSampler(*law, scales[1])
    compact = CompactVertexSampler(*law, scales[1])
    np.testing.assert_array_equal(compact.cdf, ordinary.cdf)
    left, right = np.random.default_rng(33017), np.random.default_rng(33017)
    for count in [0, 1, 64, 8000]:
        np.testing.assert_array_equal(ordinary.sample(count, left), compact.sample(count, right))
        assert left.bit_generator.state == right.bit_generator.state
    old_bytes = sum(value.nbytes for value in vars(ordinary).values()
                    if isinstance(value, np.ndarray))
    assert compact.nbytes < old_bytes/2


@pytest.mark.parametrize('values', [
    np.zeros(100), np.tile([.123456789012345, .3], 50), np.arange(100, dtype=float),
])
def test_dictionary_has_no_quantization(values):
    compact = _ExactValues(values)
    order = np.random.default_rng(23).permutation(len(values))
    np.testing.assert_array_equal(compact.gather(order), values[order])
    assert compact.nbytes <= values.nbytes


def test_cache_evicts_by_bytes_and_preserves_lru():
    cache = CompactMatchedKernel.__new__(CompactMatchedKernel)
    cache._cache = OrderedDict()
    cache._resident_bytes = 0
    cache.cache_bytes, cache.cache_size = 25, 3
    cache._cache_metrics = {'hits': 0, 'misses': 0, 'evictions': 0,
                            'build_seconds': 0., 'peak_resident_bytes': 0}
    class Sampler:
        nbytes = 10
    cache._sampler_factory = lambda index: Sampler()
    first = cache._sampler(0)
    cache._sampler(1)
    assert cache._sampler(0) is first
    cache._sampler(2)
    assert list(cache._cache) == [0, 2]
    assert cache.cache_statistics()['resident_bytes'] == 20
    assert cache.cache_statistics()['evictions'] == 1
    cache.cache_bytes = 5
    with pytest.raises(MemoryError):
        cache._sampler(3)


def test_full_mixture_and_event_stream_are_bitwise_unchanged():
    p = ROOT/'analysis/matched/protocol.json'
    spec = json.loads(p.read_text())
    args = (ROOT/spec['host_path'], ROOT/'results/matched-residence/kernel-fine.npz', p)
    old = MatchedKernel(*args)
    compact = CompactMatchedKernel(*args, cache_bytes=128*1024**2, cache_vertices=32)
    states = np.tile(np.array([[.63, 0., 0., -.81, .43, .21],
                               [0., 0., 0., .1, -.4, .7],
                               [0., .93, 0., 0., 0., 0.]]), (300, 1))
    for name, value in old.moments(states).items():
        np.testing.assert_array_equal(value, compact.moments(states)[name])
    left, right = np.random.default_rng(774), np.random.default_rng(774)
    for _ in range(2):
        np.testing.assert_array_equal(old.sample_events(states, left),
                                      compact.sample_events(states, right))
        assert left.bit_generator.state == right.bit_generator.state
    assert compact.cache_statistics()['hits'] > 0
