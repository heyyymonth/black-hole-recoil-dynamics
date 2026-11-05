"""Lossless, byte-bounded cache for the frozen matched encounter kernel.

The original float64 CDF and ordering of its positive species/node entries
are retained. Duplicate node vectors, masses and impact bounds are stored
once and indexed without quantization. Sampling uses the same random calls
and arithmetic as ``VertexSampler``. This changes storage, not the numerical
jump measure, its symmetry averaging, or the interpolation table.
"""
from collections import OrderedDict
import time

import numpy as np

from analysis.matched.kernel import (
    MatchedKernel, scatter_impulse, symmetry_samples, vertex_law,
)


def _unsigned_dtype(maximum):
    for dtype in (np.uint8, np.uint16, np.uint32, np.uint64):
        if maximum <= np.iinfo(dtype).max:
            return dtype
    raise OverflowError('Index exceeds unsigned 64-bit range')


class _ExactValues:
    """Lossless dictionary when useful; otherwise retain the float64 vector."""
    def __init__(self, values):
        values = np.asarray(values, dtype=np.float64)
        unique, inverse = np.unique(values, return_inverse=True)
        self.size = len(values)
        self.values, self.indices = values, None
        self.constant = len(unique) == 1
        if self.constant:
            self.values = unique
        else:
            inverse = inverse.astype(_unsigned_dtype(len(unique)-1))
            if unique.nbytes + inverse.nbytes < values.nbytes:
                self.values, self.indices = unique, inverse

    @property
    def nbytes(self):
        return self.values.nbytes + (0 if self.indices is None else self.indices.nbytes)

    def gather(self, indices):
        if self.constant:
            return np.full(len(indices), self.values[0], dtype=np.float64)
        if self.indices is None:
            return self.values[indices]
        return self.values[self.indices[indices]]


class CompactVertexSampler:
    """The same positive-quadrature CDF as VertexSampler, stored compactly."""
    def __init__(self, bath, velocity, kind, rotation, velocity_scale):
        relative, lower, upper = bath._state(velocity)
        rates = (bath.number_weights * np.linalg.norm(relative, axis=-1)
                 * np.pi * (upper-lower) * (upper+lower))
        positive = np.flatnonzero(rates.ravel() > 0)
        if not len(positive):
            raise ValueError('A sampled vertex must have positive allowed rate')
        species, nodes = np.unravel_index(positive, rates.shape)
        self.cdf = np.cumsum(rates.ravel()[positive])
        self.cdf /= self.cdf[-1]
        self.cdf[-1] = 1.
        self.positive = positive.astype(_unsigned_dtype(int(positive[-1])))
        self.node_count = rates.shape[1]
        self.relative_nodes = np.array(relative[0], copy=True)
        self.mass_values = np.array(bath.masses, copy=True)
        self.lower = _ExactValues(lower[species, nodes])
        self.upper = _ExactValues(upper[species, nodes])
        self.remnant, self.G = bath.remnant_mass, bath.G
        self.kind = kind
        self.rotation = np.array(rotation, copy=True)
        self.velocity_scale = velocity_scale

    @property
    def nbytes(self):
        """Resident numerical-array bytes; Python object overhead is excluded."""
        return (self.cdf.nbytes + self.positive.nbytes + self.relative_nodes.nbytes
                + self.mass_values.nbytes + self.lower.nbytes + self.upper.nbytes
                + self.rotation.nbytes)

    def sample(self, count, rng):
        index = np.searchsorted(self.cdf, rng.random(count), side='right')
        packed = self.positive[index]
        species, nodes = packed // self.node_count, packed % self.node_count
        low, high = self.lower.gather(index), self.upper.gather(index)
        impact = np.sqrt(low*low + rng.random(count)*(high-low)*(high+low))
        increments = scatter_impulse(self.relative_nodes[nodes], self.mass_values[species],
                                     self.remnant, self.G, impact,
                                     rng.uniform(0., 2*np.pi, count))
        increments = symmetry_samples(increments, self.kind, rng)
        return increments @ self.rotation.T / self.velocity_scale


class CompactMatchedKernel(MatchedKernel):
    """Unchanged matched-law API with a byte- and vertex-bounded LRU cache.

    ``cache_bytes`` limits retained sampler arrays, not temporary construction
    arrays or process memory. Callers must budget for those and for concurrent
    kernels. Table/source validation is inherited from MatchedKernel. The
    adapter's source hash and performance policy belong in the new analysis
    provenance; the existing kernel/table files need no changes.
    """
    def __init__(self, host_path, table_path, protocol_path,
                 cache_bytes=2*1024**3, cache_vertices=768):
        for name, value in [('cache_bytes', cache_bytes), ('cache_vertices', cache_vertices)]:
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
                raise ValueError(f'{name} must be a positive integer')
        super().__init__(host_path, table_path, protocol_path, cache_size=int(cache_vertices))
        self.cache_bytes = int(cache_bytes)
        self._cache = OrderedDict()
        self._resident_bytes = 0
        self._cache_metrics = {'hits': 0, 'misses': 0, 'evictions': 0,
                               'build_seconds': 0., 'peak_resident_bytes': 0}

    def _physical_sampler(self, index):
        coordinates = [self.axes[i][j] for i, j in
                       enumerate(np.unravel_index(index, self.shape))]
        scales = (self.r_scale_pc, self.v_scale_kms, self.time_scale_Myr)
        answer, law = vertex_law(self.host, coordinates, scales, self.protocol,
                                 self.encounter_protocol, int(self.orders[index]))
        for key, factor in [('rate', self.time_scale_Myr),
                            ('excluded', self.time_scale_Myr),
                            ('A', self.time_scale_Myr/self.v_scale_kms),
                            ('D', self.time_scale_Myr/self.v_scale_kms**2)]:
            if not np.allclose(np.asarray(answer[key])*factor, getattr(self, key)[index],
                               rtol=2e-11, atol=2e-13):
                raise ValueError('Rebuilt vertex law differs from its recorded moments')
        return CompactVertexSampler(*law, self.v_scale_kms)

    def _sampler(self, index):
        index = int(index)
        if index in self._cache:
            self._cache.move_to_end(index)
            self._cache_metrics['hits'] += 1
            return self._cache[index]
        self._cache_metrics['misses'] += 1
        started = time.perf_counter()
        sampler = self._sampler_factory(index)
        self._cache_metrics['build_seconds'] += time.perf_counter()-started
        if sampler.nbytes > self.cache_bytes:
            raise MemoryError('A single sampler exceeds the declared cache byte budget')
        while self._cache and (self._resident_bytes + sampler.nbytes > self.cache_bytes
                               or len(self._cache) >= self.cache_size):
            _, removed = self._cache.popitem(last=False)
            self._resident_bytes -= removed.nbytes
            self._cache_metrics['evictions'] += 1
        self._cache[index] = sampler
        self._resident_bytes += sampler.nbytes
        self._cache_metrics['peak_resident_bytes'] = max(
            self._cache_metrics['peak_resident_bytes'], self._resident_bytes)
        return sampler

    def cache_statistics(self):
        return {**self._cache_metrics, 'resident_bytes': self._resident_bytes,
                'resident_vertices': len(self._cache), 'maximum_bytes': self.cache_bytes,
                'maximum_vertices': self.cache_size,
                'bytes_scope': 'Retained sampler arrays; excludes Python overhead and temporary construction arrays',
                'law': 'Original float64 CDF, exact index/bound encoding, unchanged draw order'}
