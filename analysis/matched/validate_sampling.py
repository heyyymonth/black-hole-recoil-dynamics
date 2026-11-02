"""Independent off-grid sampling check of the frozen positive-mixture law.

Write the declared states and source hashes before drawing any increments.
This checks the event sampler against its own numerical mixture moments;
it does not establish convergence to the continuous host encounter law.
"""
from pathlib import Path
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from analysis.matched.kernel import MatchedKernel, canonical_velocity, digest


SEED = 2026110101
SAMPLES_PER_STATE = 40000
OUTPUT = ROOT / 'results/matched-residence'
PROTOCOL = ROOT / 'analysis/matched/protocol.json'


def declared_states():
    """Fixed interior and symmetry controls, chosen without sampled outcomes."""
    rotation = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]])
    settings = [
        ('generic_outward', .63, 1.17, .31, 1., np.eye(3)),
        ('generic_inward_rotated', 1.17, .83, .61, -1., rotation),
        ('near_radial_rotated', .37, 1.41, .0001, 1., rotation),
        ('near_tangential', 1.83, .37, .9999, 1., np.eye(3)),
        ('center_moving_rotated', 0., .83, 0., 1., rotation),
        ('stationary_off_center', .93, 0., 0., 1., np.eye(3)),
    ]
    records = []
    for label, radius, speed, pitch, sign, frame in settings:
        position = frame @ np.array([0., 0., radius])
        velocity = canonical_velocity(speed, pitch)
        velocity[2] *= sign
        state = np.r_[position, frame @ velocity]
        records.append({'label': label, 'state': state.tolist(),
                        'radius': radius, 'speed': speed, 'pitch': pitch,
                        'radial_sign': sign})
    return records


def compare(sample, expected):
    estimate = sample.mean(axis=0)
    se = sample.std(axis=0, ddof=1) / np.sqrt(len(sample))
    residual = estimate - expected
    threshold = 6. * se + 1e-13
    return {'expected': np.asarray(expected).tolist(),
            'estimate': estimate.tolist(),
            'standard_error': se.tolist(),
            'residual': residual.tolist(),
            'standardized_residual': np.divide(residual, se,
                out=np.zeros_like(residual), where=se > 0).tolist(),
            'passed': bool(np.all(abs(residual) < threshold))}


def run():
    protocol = json.loads(PROTOCOL.read_text())
    host = ROOT / protocol['host_path']
    sources = [PROTOCOL, Path(__file__), host,
               ROOT / protocol['encounter_protocol']]
    for grid in ('coarse', 'fine'):
        table = OUTPUT / f'kernel-{grid}.npz'
        metadata = json.loads(table.with_suffix('.json').read_text())
        sources.extend([table, table.with_suffix('.json')])
        sources.extend(ROOT / path for path in metadata['implementation_sha256'])
    hashes = {str(path.relative_to(ROOT)): digest(path)
              for path in sorted(set(sources))}
    design = {
        'revision': 'matched-sampling-validation-1',
        'seed': SEED, 'samples_per_state_per_grid': SAMPLES_PER_STATE,
        'grids': ['coarse', 'fine'], 'states': declared_states(),
        'units': 'dimensionless kernel r_rho, v_dyn, t_dyn',
        'random_streams': 'Independent SeedSequence([seed, grid index, state index]); independent of trajectory streams',
        'acceptance': 'Every mean and raw-second component within six estimated Monte Carlo standard errors plus 1e-13 absolute roundoff floor',
        'interpretation': 'A reproducible sampler-versus-table diagnostic, not a calibrated simultaneous hypothesis test or a claim of continuous-law or trajectory convergence',
        'source_sha256': hashes,
    }
    design_path = OUTPUT / 'sampling-validation-design.json'
    design_path.write_text(json.dumps(design, indent=2) + '\n')
    # The complete design is on disk before any random draws occur.
    records = []
    start = time.perf_counter()
    for grid_index, grid in enumerate(design['grids']):
        kernel = MatchedKernel(host, OUTPUT / f'kernel-{grid}.npz', PROTOCOL)
        for state_index, record in enumerate(design['states']):
            state = np.asarray(record['state'])[None, :]
            moment = kernel.moments(state)
            rate = moment['rate'][0]
            rng = np.random.default_rng(np.random.SeedSequence([SEED, grid_index, state_index]))
            draws = kernel.sample_events(np.repeat(state, SAMPLES_PER_STATE, axis=0), rng)
            checks = {
                'mean': compare(draws, moment['A'][0] / rate),
                'raw_second': compare(draws[:, :, None] * draws[:, None, :],
                                      moment['D'][0] / rate),
            }
            item = {'grid': grid, 'label': record['label'],
                    'rate_per_t_dyn': float(rate), **checks}
            item['passed'] = all(check['passed'] for check in checks.values())
            records.append(item)
            print(grid, record['label'], 'pass' if item['passed'] else 'FAIL', flush=True)
        del kernel
    changed = [path for path, sha in hashes.items() if digest(ROOT / path) != sha]
    result = {'revision': design['revision'], 'design_sha256': digest(design_path),
              'source_sha256': hashes, 'sources_changed_during_run': changed,
              'records': records, 'all_checks_passed': not changed and all(r['passed'] for r in records),
              'elapsed_seconds': time.perf_counter() - start}
    output = OUTPUT / 'sampling-validation.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    if not result['all_checks_passed']:
        raise SystemExit('Sampler validation failed; retain diagnostics')
    print(output.relative_to(ROOT), flush=True)


if __name__ == '__main__':
    run()
