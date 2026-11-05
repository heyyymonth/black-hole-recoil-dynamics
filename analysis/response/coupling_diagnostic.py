"""Deterministic exploratory checks of Gaussian common-noise factor continuity.

Generation queries existing local moment tables but never evolves trajectories.
Verification uses NumPy and retained matrices only. Eigenbasis factors are not
unique, so verification checks their retained eigendecomposition and covariance
rather than requiring a platform to reproduce arbitrary eigenvector signs.
"""
from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/response/coupling-diagnostic.json'
PROTOCOL = ROOT/'analysis/response/protocol.json'
SYNTHETIC_EPS = [1e-3, 1e-6, 1e-9]
LOCAL_EPS = [1e-4, 1e-8, 1e-12]
DIRECTIONS = [('v_y', 4), ('v_z', 5)]
SOURCE_FILES = [
    'analysis/response/coupling_diagnostic.py', 'analysis/matched/kernel.py',
    'analysis/matched/transport.py', 'analysis/host_reconstruction/host.py',
    'analysis/encounters/physics.py', 'analysis/encounters/run.py',
]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def factors(matrix):
    matrix = np.asarray(matrix, dtype=float)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    if np.any(eigenvalues <= 0):
        raise ValueError('This diagnostic requires positive definite matrices')
    factor = eigenvectors*np.sqrt(eigenvalues)[None, :]
    root = factor@eigenvectors.T
    return {'D': matrix.tolist(), 'eigenvalues': eigenvalues.tolist(),
            'eigenvectors': eigenvectors.tolist(), 'eigenfactor': factor.tolist(),
            'symmetric_root': root.tolist()}


def contrast(first, second):
    D0, D1 = [np.asarray(record['D']) for record in (first, second)]
    F0, F1 = [np.asarray(record['eigenfactor']) for record in (first, second)]
    R0, R1 = [np.asarray(record['symmetric_root']) for record in (first, second)]
    return {'D_frobenius_difference': float(np.linalg.norm(D1-D0)),
            'D_relative_frobenius_difference': float(np.linalg.norm(D1-D0)/np.linalg.norm(D0)),
            'eigenfactor_frobenius_difference': float(np.linalg.norm(F1-F0)),
            'symmetric_root_frobenius_difference': float(np.linalg.norm(R1-R0)),
            'eigenfactor_common_noise_mse_per_unit_dt': float(np.sum((F1-F0)**2)),
            'symmetric_root_common_noise_mse_per_unit_dt': float(np.sum((R1-R0)**2)),
            'eigenfactor_common_noise_mse_over_trace_D0': float(np.sum((F1-F0)**2)/np.trace(D0)),
            'symmetric_root_common_noise_mse_over_trace_D0': float(np.sum((R1-R0)**2)/np.trace(D0))}


def initial_state(kernel, matched, launch):
    radius = matched['launch']['radius_over_density_radius']
    apocenter = matched['launch']['radial_apocenter_over_density_radius']
    speed = np.sqrt(2*(kernel.host.psi(radius*kernel.r_scale_pc)
                      - kernel.host.psi(apocenter*kernel.r_scale_pc)))/kernel.v_scale_kms
    angle = launch['angle_from_radial_degrees']
    if angle == 0:
        direction = np.array([1., 0., 0.])
    elif angle == 90:
        direction = np.array([0., 1., 0.])
    else:
        direction = np.array([np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle)), 0.])
    return np.r_[radius, 0., 0., speed*direction]


def generate():
    # Keep SciPy-dependent reconstruction imports out of the verification path.
    sys.path.insert(0, str(ROOT))
    from analysis.matched.kernel import MatchedKernel
    p = json.loads(PROTOCOL.read_text())
    matched = json.loads((ROOT/p['matched_protocol']).read_text())
    table = ROOT/p['kernel_tables']['fine']
    kernel = MatchedKernel(ROOT/p['host_path'], table, ROOT/p['matched_protocol'])
    input_names = [str(PROTOCOL.relative_to(ROOT)), p['matched_protocol'], p['host_path'],
                   matched['encounter_protocol'], p['kernel_tables']['fine'],
                   str(table.with_suffix('.json').relative_to(ROOT))]
    # These hashes establish that the local check followed the retained pilot.
    pilot_names = [f'results/response/pilot-fine-dt0.010-{launch["name"]}-weighted.json'
                   for launch in p['launches']]
    if not all((ROOT/name).exists() for name in pilot_names):
        raise ValueError('The declared post-pilot diagnostic requires completed pilot records')
    ids = {name: digest(ROOT/name) for name in input_names+pilot_names+SOURCE_FILES}
    synthetic = []
    for eps in SYNTHETIC_EPS:
        minus = factors(np.diag([1-eps, 1+eps, 2.]))
        plus = factors(np.diag([1+eps, 1-eps, 2.]))
        synthetic.append({'epsilon': eps, 'minus': minus, 'plus': plus,
                          'contrast': contrast(minus, plus)})
    local = []
    for launch in p['launches']:
        state = initial_state(kernel, matched, launch)
        base = factors(kernel.moments(state[None, :])['D'][0])
        probes = []
        for label, component in DIRECTIONS:
            for eps in LOCAL_EPS:
                perturbed = state.copy()
                perturbed[component] += eps
                other = factors(kernel.moments(perturbed[None, :])['D'][0])
                probes.append({'direction': label, 'component': component, 'epsilon': eps,
                               'state': perturbed.tolist(), 'factors': other,
                               'contrast': contrast(base, other)})
        local.append({'launch': launch['name'], 'state': state.tolist(),
                      'base': base, 'probes': probes})
    if any(digest(ROOT/name) != sha for name, sha in ids.items()):
        raise ValueError('Diagnostic sources or inputs changed during evaluation')
    return {'revision': 'gaussian-coupling-local-diagnostic-1',
            'created_utc': datetime.now(timezone.utc).isoformat(),
            'scope': 'Post-pilot exploratory deterministic local-factor check; no stochastic trajectories; probe choices followed preliminary local inspection and are not blind validation',
            'design': {'synthetic_eps': SYNTHETIC_EPS, 'local_eps': LOCAL_EPS,
                       'directions': [{'label': label, 'component': component} for label, component in DIRECTIONS],
                       'local_grid': 'fine', 'launches': [launch['name'] for launch in p['launches']]},
            'factor_definition': 'F=V diag(sqrt(lambda)); R=V diag(sqrt(lambda)) V^T; F F^T=R R^T=D',
            'coupling_identity': 'For a common standard normal Z, E||sqrt(dt)*(F1-F0)*Z||^2=dt*||F1-F0||_F^2; deterministic drift is excluded from this noise-only metric',
            'synthetic': synthetic, 'local': local, 'input_and_source_sha256': ids,
            'numpy_version': np.__version__,
            'interpretation': 'Eigenvector signs and bases can change abruptly while D changes little. This can weaken common-random-number coupling without changing either marginal Gaussian transition law. Local continuity diagnostics do not attribute the measured residence-score variance to this mechanism or prove an improved estimator meets precision targets.',
            'verification_scope': 'Check retained matrix decompositions, equal marginal covariance, analytic common-noise MSE, probe inventory and source hashes; verification does not reconstruct cluster moments or evolve trajectories.'}


def check_factor(record):
    D = np.asarray(record['D'])
    values = np.asarray(record['eigenvalues'])
    V = np.asarray(record['eigenvectors'])
    F = np.asarray(record['eigenfactor'])
    R = np.asarray(record['symmetric_root'])
    if any(x.shape != (3, 3) or not np.isfinite(x).all() for x in (D, V, F, R)):
        raise ValueError('Invalid retained matrix')
    if values.shape != (3,) or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('Invalid retained spectrum')
    for actual, expected in [(D, D.T), (V.T@V, np.eye(3)), (D@V, V*values),
                             (F, V*np.sqrt(values)), (R, F@V.T),
                             (F@F.T, D), (R@R.T, D), (R, R.T)]:
        np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=2e-14)
    # The symmetric positive root is unique even where eigenvector bases are not.
    values2, V2 = np.linalg.eigh(D)
    np.testing.assert_allclose(R, (V2*np.sqrt(values2))@V2.T, rtol=2e-11, atol=2e-14)


def check_contrast(record, first, second):
    expected = contrast(first, second)
    if set(record) != set(expected):
        raise ValueError('Contrast field inventory differs')
    for key, value in expected.items():
        np.testing.assert_allclose(record[key], value, rtol=2e-11, atol=1e-30)
    for label, field in [('eigenfactor', 'eigenfactor'), ('symmetric_root', 'symmetric_root')]:
        difference = np.asarray(second[field])-np.asarray(first[field])
        # Covariance trace independently evaluates the common-noise MSE identity.
        mse = float(np.trace(difference@difference.T))
        np.testing.assert_allclose(record[label+'_common_noise_mse_per_unit_dt'], mse,
                                   rtol=2e-11, atol=1e-30)


def verify(record):
    if record['revision'] != 'gaussian-coupling-local-diagnostic-1':
        raise ValueError('Unexpected diagnostic revision')
    for name, sha in record['input_and_source_sha256'].items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts or digest(ROOT/path) != sha:
            raise ValueError('Stale or unsafe diagnostic provenance: '+name)
    p = json.loads(PROTOCOL.read_text())
    expected_design = {'synthetic_eps': SYNTHETIC_EPS, 'local_eps': LOCAL_EPS,
                       'directions': [{'label': label, 'component': component} for label, component in DIRECTIONS],
                       'local_grid': 'fine', 'launches': [launch['name'] for launch in p['launches']]}
    if record['design'] != expected_design:
        raise ValueError('Declared diagnostic design differs')
    if [row['epsilon'] for row in record['synthetic']] != SYNTHETIC_EPS:
        raise ValueError('Synthetic probe inventory differs')
    for row in record['synthetic']:
        eps = row['epsilon']
        np.testing.assert_array_equal(row['minus']['D'], np.diag([1-eps, 1+eps, 2.]))
        np.testing.assert_array_equal(row['plus']['D'], np.diag([1+eps, 1-eps, 2.]))
        check_factor(row['minus']); check_factor(row['plus'])
        check_contrast(row['contrast'], row['minus'], row['plus'])
        np.testing.assert_allclose(row['contrast']['eigenfactor_frobenius_difference'], 2., rtol=2e-11)
        analytic_root_difference = np.sqrt(2.)*2*eps/(np.sqrt(1+eps)+np.sqrt(1-eps))
        np.testing.assert_allclose(row['contrast']['symmetric_root_frobenius_difference'],
                                   analytic_root_difference, rtol=2e-7)
    if [row['launch'] for row in record['local']] != [launch['name'] for launch in p['launches']]:
        raise ValueError('Local launch inventory differs')
    expected_probes = [(label, component, eps) for label, component in DIRECTIONS for eps in LOCAL_EPS]
    for row in record['local']:
        check_factor(row['base'])
        if [(probe['direction'], probe['component'], probe['epsilon']) for probe in row['probes']] != expected_probes:
            raise ValueError('Local perturbation inventory differs')
        for probe in row['probes']:
            expected_state = np.asarray(row['state']).copy()
            expected_state[probe['component']] += probe['epsilon']
            np.testing.assert_array_equal(probe['state'], expected_state)
            check_factor(probe['factors'])
            check_contrast(probe['contrast'], row['base'], probe['factors'])
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.verify:
        verify(json.loads(OUT.read_text()))
        print('Verified retained covariance and local common-noise contrasts; no physical reconstruction or trajectories.')
    else:
        result = generate()
        verify(result)
        OUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
        print('Recorded 3 synthetic and 18 local deterministic coupling comparisons; no trajectories.')


if __name__ == '__main__':
    main()
