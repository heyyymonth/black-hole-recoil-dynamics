"""Continue the fixed experiment with lossless on-disk sampler reuse.

The original protocol, numerical implementation, seeds, sample allocations and
completed records are immutable. This driver adds a separately hashed storage
adapter, with its exact remaining job inventory fixed before use. It wraps the
record writer only to attach adapter provenance before atomic publication;
all trajectories execute through the original Runner and transport functions.
"""
from pathlib import Path
import argparse
import gc
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from analysis.response import run as execution
from analysis.response.disk_cache import DiskMatchedKernel

OUT = execution.OUT
MANIFEST = OUT/'storage-adapter.json'
SOURCE_FILES = ('analysis/response/execute_cached.py', 'analysis/response/disk_cache.py',
                'analysis/response/test_disk_cache.py')


def adapter_manifest(p, ids, counts):
    """Record storage-only implementation changes without rewriting earlier evidence."""
    all_jobs = [job for stage in ('pilot', 'prediction', 'validation')
                for job in execution.jobs(p, stage, counts)]
    names = {job['name'] for job in all_jobs}
    sources = {name: execution.digest(ROOT/name) for name in SOURCE_FILES}
    inputs = {str(path.relative_to(ROOT)): execution.digest(path)
              for path in (OUT/'execution-design.json', OUT/'sample-count-lock.json')}
    if MANIFEST.exists():
        manifest = execution.read_json(MANIFEST)
        if (manifest['frozen_identity'] != ids or manifest['source_sha256'] != sources
                or manifest['input_sha256'] != inputs):
            raise ValueError('Storage adapter sources or frozen inputs changed')
        existing = manifest['existing_record_sha256']
        if any(execution.digest(ROOT/name) != sha for name, sha in existing.items()):
            raise ValueError('A pre-adapter record changed')
        old_names = {Path(name).stem for name in existing}
        adapted_names = set(manifest['adapted_jobs'])
        if (old_names | adapted_names != names or old_names & adapted_names
                or len(adapted_names) != len(manifest['adapted_jobs'])):
            raise ValueError('Storage adapter job inventory changed')
        evidence = manifest['test_evidence']
        if execution.digest(ROOT/evidence['path']) != evidence['sha256'] or evidence['returncode'] != 0:
            raise ValueError('Storage equivalence evidence changed')
        return manifest
    existing, adapted = {}, []
    for job in all_jobs:
        path = OUT/(job['name']+'.json')
        if path.exists():
            record = execution.read_json(path)
            if (not record['complete'] or record['identity'] != ids
                    or execution.digest(path.with_suffix('.npz')) != record['samples_sha256']):
                raise ValueError('Invalid existing result before adapter initialization')
            existing[str(path.relative_to(ROOT))] = execution.digest(path)
        else:
            if path.with_suffix('.npz').exists():
                raise ValueError('Unregistered arrays exist; resolve before continuation')
            adapted.append(job['name'])
    completed = subprocess.run([sys.executable, '-m', 'pytest', '-q',
        'analysis/response/test_disk_cache.py'], cwd=ROOT, text=True, capture_output=True)
    evidence = OUT/'storage-adapter-validation.txt'
    if evidence.exists():
        raise FileExistsError('Storage test evidence exists without its manifest')
    evidence.write_text(completed.stdout+completed.stderr)
    if completed.returncode:
        raise RuntimeError(completed.stdout+completed.stderr)
    if {name: execution.digest(ROOT/name) for name in SOURCE_FILES} != sources:
        raise ValueError('Storage adapter sources changed during equivalence checks')
    manifest = {'revision': 'lossless-disk-storage-1', 'created_utc': execution.utc_now(),
        'frozen_identity': ids, 'source_sha256': sources, 'input_sha256': inputs,
        'existing_record_sha256': existing, 'adapted_jobs': adapted,
        'test_evidence': {'path': str(evidence.relative_to(ROOT)),
            'sha256': execution.digest(evidence), 'returncode': completed.returncode},
        'scope': 'Post-pilot storage optimization only. Exact numeric arrays and sampler RNG calls are preserved; no transition law, seed, sample size or statistical threshold changes.',
        'storage': {'retained_array_bytes': p['cache']['retained_array_bytes'],
            'vertex_cap': p['cache']['vertex_cap'], 'disk_max_bytes_per_grid': 15*1024**3,
            'location': 'tmp/response-sampler-cache; disposable, reconstructed from frozen physical inputs'}}
    execution.write_new_json(MANIFEST, manifest)
    return manifest


class CachedRunner(execution.Runner):
    def __init__(self, p, ids, resume, manifest):
        super().__init__(p, ids, resume)
        self.adapter = manifest
        self.adapter_ref = {'path': str(MANIFEST.relative_to(ROOT)),
                            'sha256': execution.digest(MANIFEST)}

    def get_kernel(self, grid):
        if self.grid != grid:
            self.kernel = None
            gc.collect()
            table = ROOT/self.p['kernel_tables'][grid]
            metadata = execution.read_json(table.with_suffix('.json'))
            if not (metadata['completed'] and metadata['all_vertex_quadrature_passed']
                    and metadata['all_holdout_quadrature_passed']):
                raise ValueError('Accepted quadrature and holdouts required')
            self.kernel = DiskMatchedKernel(ROOT/self.p['host_path'], table,
                ROOT/self.p['matched_protocol'],
                disk_cache_dir=ROOT/'tmp/response-sampler-cache'/grid,
                namespace=execution.digest(OUT/'execution-design.json'), grid=grid,
                cache_bytes=self.p['cache']['retained_array_bytes'],
                cache_vertices=self.p['cache']['vertex_cap'], disk_max_bytes=15*1024**3)
            self.grid = grid
        return self.kernel

    def execute(self, job, dependencies):
        path = OUT/(job['name']+'.json')
        adapted = job['name'] in self.adapter['adapted_jobs']
        if path.exists():
            r = execution.read_json(path)
            if (adapted and r.get('storage_adapter') != self.adapter_ref) or (not adapted and 'storage_adapter' in r):
                raise ValueError('Result does not match the declared storage-adapter inventory')
        elif not adapted:
            raise ValueError('A pre-adapter result is missing')
        original_writer = execution.write_new_json

        def attach_storage_provenance(target, record):
            if Path(target) == path:
                if execution.digest(MANIFEST) != self.adapter_ref['sha256']:
                    raise ValueError('Storage adapter manifest changed during execution')
                if any(execution.digest(ROOT/name) != sha
                       for name, sha in self.adapter['source_sha256'].items()):
                    raise ValueError('Storage adapter source changed during execution')
                storage = self.kernel.disk_statistics()
                storage['cache_path'] = str(Path(storage['cache_path']).relative_to(ROOT))
                record = {**record, 'storage_adapter': self.adapter_ref,
                          'disk_cache': storage}
            original_writer(target, record)

        execution.write_new_json = attach_storage_provenance
        try:
            super().execute(job, dependencies)
        finally:
            execution.write_new_json = original_writer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['prediction', 'validation', 'all'], default='all')
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    p = execution.read_json(execution.PROTOCOL); ids = execution.identity(p)
    design_sha = execution.ensure_design(p, ids)
    counts, count_sha = execution.lock_counts(p, ids, design_sha)
    manifest = adapter_manifest(p, ids, counts)
    runner = CachedRunner(p, ids, args.resume, manifest)
    deps = {'execution_design_sha256': design_sha, 'sample_count_lock_sha256': count_sha}
    if args.stage in ('prediction', 'all'):
        for job in execution.jobs(p, 'prediction', counts):
            runner.execute(job, deps)
    register_sha = execution.lock_predictions(p, ids, counts, design_sha, count_sha)
    if args.stage in ('validation', 'all'):
        for job in execution.jobs(p, 'validation'):
            runner.execute(job, {**deps, 'prediction_register_sha256': register_sha})
    print('completed cached stage', args.stage, flush=True)


if __name__ == '__main__':
    main()
