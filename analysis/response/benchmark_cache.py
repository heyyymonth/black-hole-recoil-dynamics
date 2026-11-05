"""Reproducible storage/runtime check; no trajectory or physical result.

The fixed cyclic workload deliberately exceeds the old 128-vertex cache.
It measures this workload only, not an asserted production speedup.
"""
import gc
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis.matched.kernel import MatchedKernel, digest
from analysis.response.cache import CompactMatchedKernel


def run():
    protocol = ROOT/'analysis/matched/protocol.json'
    specification = json.loads(protocol.read_text())
    table = ROOT/'results/matched-residence/kernel-fine.npz'
    host = ROOT/specification['host_path']
    metadata = json.loads(table.with_suffix('.json').read_text())
    source_paths = [protocol, table, table.with_suffix('.json'), host, Path(__file__),
                    ROOT/'analysis/response/cache.py', ROOT/specification['encounter_protocol']]
    source_paths.extend(ROOT/p for p in metadata['implementation_sha256'])
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in sorted(set(source_paths))}
    data = dict(np.load(table, allow_pickle=False))
    shape = tuple(len(data[name]) for name in ('radius', 'speed', 'pitch'))
    vertices = [int(np.ravel_multi_index((r, s, p), shape))
                for r in (2, 4, 6, 8) for s in (2, 4, 6, 8, 10)
                for p in (0, 1, 2, 3, 4, 5, 6, 8)]
    output = ROOT/'results/response'
    output.mkdir(parents=True, exist_ok=True)
    design = {'revision': 'compact-cache-benchmark-1',
              'source_sha256': hashes, 'vertices': vertices,
              'passes': 2, 'draws_per_visit': 16, 'seed': 2026110201,
              'ordinary_cache_vertices': 128,
              'compact_cache_vertices': 768, 'compact_cache_bytes': 2*1024**3,
              'scope': 'Synthetic cyclic cache workload; not a trajectory experiment or an estimate of production speedup'}
    design_path = output/'cache-benchmark-design.json'
    design_path.write_text(json.dumps(design, indent=2)+'\n')
    records = []
    for name in ('ordinary', 'compact'):
        factory = MatchedKernel if name == 'ordinary' else CompactMatchedKernel
        kwargs = ({} if name == 'ordinary' else
                  {'cache_bytes': design['compact_cache_bytes'],
                   'cache_vertices': design['compact_cache_vertices']})
        kernel = factory(host, table, protocol, **kwargs)
        rng = np.random.default_rng(design['seed'])
        increment_hash = hashlib.sha256()
        hits, misses = 0, 0
        pass_seconds = []
        first_sampler_bytes = None
        for sweep in range(design['passes']):
            started = time.perf_counter()
            for j, index in enumerate(vertices):
                present = index in kernel._cache
                hits += int(present)
                misses += int(not present)
                sampler = kernel._sampler(index)
                if first_sampler_bytes is None:
                    first_sampler_bytes = (sampler.nbytes if name == 'compact' else
                        sum(a.nbytes for a in vars(sampler).values() if isinstance(a, np.ndarray)))
                increment_hash.update(sampler.sample(design['draws_per_visit'], rng).tobytes())
                if (j+1) % 40 == 0:
                    print(name, 'pass', sweep+1, 'visited', j+1, flush=True)
            pass_seconds.append(time.perf_counter()-started)
        retained_bytes = (kernel.cache_statistics()['resident_bytes'] if name == 'compact' else
            sum(a.nbytes for s in kernel._cache.values() for a in vars(s).values()
                if isinstance(a, np.ndarray)))
        record = {'name': name, 'pass_seconds': pass_seconds,
                  'total_seconds': sum(pass_seconds), 'hits': hits, 'misses': misses,
                  'first_sampler_bytes': first_sampler_bytes,
                  'retained_sampler_bytes': retained_bytes,
                  'increment_sha256': increment_hash.hexdigest(),
                  'rng_state': rng.bit_generator.state}
        if name == 'compact':
            record['cache_statistics'] = kernel.cache_statistics()
        records.append(record)
        del sampler, kernel
        gc.collect()
    changed = [path for path, value in hashes.items() if digest(ROOT/path) != value]
    result = {'revision': design['revision'], 'design_sha256': digest(design_path),
              'source_sha256': hashes, 'sources_changed_during_run': changed,
              'records': records,
              'bitwise_equal_increment_stream': records[0]['increment_sha256'] == records[1]['increment_sha256'],
              'equal_final_rng_state': records[0]['rng_state'] == records[1]['rng_state'],
              'scope': design['scope']}
    result_path = output/'cache-benchmark.json'
    result_path.write_text(json.dumps(result, indent=2)+'\n')
    if changed or not result['bitwise_equal_increment_stream'] or not result['equal_final_rng_state']:
        raise SystemExit('Cache benchmark exactness/provenance failed')
    print(json.dumps(result['records'], indent=2), flush=True)


if __name__ == '__main__':
    run()
