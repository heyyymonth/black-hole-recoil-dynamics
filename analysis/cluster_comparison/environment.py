"""Six-cluster environmental controls; no remnant dynamics or new host fit.

Reuses the registered Galactic force, orbit solver, coordinate transform, and
enclosed-host representation. Source rows and per-cluster hosts remain distinct.
Run --catalog-only before reconstructed hosts exist; the default adds their
conditional boundary and tidal diagnostics. Missing or rejected hosts are explicit.
"""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import html
import json
import re
import sys

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
for directory in ['analysis/host_reconstruction', 'analysis/boundary', 'analysis/orbit']:
    sys.path.insert(0, str(ROOT / directory))
from host import Host
from model import EnclosedHost, max_force_ratio
from galaxy import Galaxy
from orbits import K, initial_state, integrate, conserved, angular_fields

CLUSTERS = ['NGC0104', 'NGC1851', 'NGC3201', 'NGC5139', 'NGC6121', 'NGC6397']
DATA = ROOT / 'data/cluster-comparison-environment'
OUT = ROOT / 'results/cluster-comparison-environment'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def extract_source(raw, cluster, anchor):
    """Extract exact table row by source ID, never by a substring of its name."""
    matches = [row for row in raw.splitlines()
               if re.search(r"data-sort=['\"]\s*" + cluster + r"['\"]", row)]
    if len(matches) != 1:
        raise ValueError(f'{cluster}: expected one source row, got {len(matches)}')
    row = matches[0]
    cells = [html.unescape(re.sub('<[^>]*>', '', cell)).strip()
             for cell in re.findall(r'<td[^>]*>(.*?)(?=<td|$)', row)]
    if len(cells) != 19:
        raise ValueError(f'{cluster}: expected 19 source columns')
    values = [float(cell.split('±')[0]) for cell in cells[1:]]
    if not np.isfinite(values).all():
        raise ValueError('Nonfinite source value')
    source = {k: anchor[k] for k in ['url', 'captured_date', 'raw_path', 'raw_sha256',
                                    'columns', 'frame', 'galaxy_source', 'limitations']}
    source.update(cluster=cluster, row_cells=cells,
                  baseline_point_kpc_kms=values[10:16], RA_deg=values[0],
                  DEC_deg=values[1], distance_kpc=values[4], R_GC_kpc=values[5],
                  radial_velocity_kms=values[6], pmra_cosdec_masyr=values[7],
                  pmdec_masyr=values[8], pm_correlation=values[9],
                  catalog_mean_pericenter_kpc=values[16],
                  catalog_mean_apocenter_kpc=values[17])
    for key, index in [('distance_error_kpc', 5), ('radial_velocity_error_kms', 7),
                       ('pmra_error_masyr', 8), ('pmdec_error_masyr', 9)]:
        source[key] = float(cells[index].split('±')[1])
    if source['distance_kpc'] <= source['distance_error_kpc'] or abs(source['pm_correlation']) > 1:
        raise ValueError('Invalid distance or proper-motion correlation')
    return source, row


def prepare_sources():
    anchor = json.loads((ROOT / 'data/orbit/ngc3201-source.json').read_text())
    rawpath = ROOT / anchor['raw_path']
    if rawpath.exists():
        if digest(rawpath) != anchor['raw_sha256']:
            raise ValueError('Pinned orbit-source snapshot changed')
        raw = rawpath.read_text()
    else:
        # A clean checkout contains the exact six rows, not the entire cached
        # upstream page. Validate both extracts and parsed metadata against the
        # recorded provenance before using this equivalent source route.
        pins = json.loads((OUT / 'provenance.json').read_text())['input_sha256']
        raw = None
    sources = {}
    for cluster in CLUSTERS:
        rowpath, sourcepath = DATA / f'{cluster}-row.html', DATA / f'{cluster}.json'
        if raw is None:
            for path in (rowpath, sourcepath):
                if digest(path) != pins[str(path.relative_to(ROOT))]:
                    raise ValueError(f'Pinned orbit-source extract changed: {path.name}')
            source, row = extract_source(rowpath.read_text(), cluster, anchor)
            if source != json.loads(sourcepath.read_text()):
                raise ValueError(f'Orbit-source row and metadata disagree: {cluster}')
        else:
            source, row = extract_source(raw, cluster, anchor)
            save_json(sourcepath, source)
            rowpath.write_text(row + '\n')
        sources[cluster] = source
    return sources


def protocol():
    original = json.loads((ROOT / 'analysis/orbit/protocol.json').read_text())
    value = {k: original[k] for k in ['galaxy', 'units', 'integration', 'force_accuracy',
                                    'acceptance', 'tensor_derivative_step_Myr']}
    value.update(revision='six-cluster-environment-1', clusters=CLUSTERS,
                 source_snapshot='data/orbit/raw/uq_orbits.html',
                 host_path_pattern='results/cluster-comparison/{cluster}/reference/host.npz',
                 cases=['catalog_reference', 'host_distance',
                        'catalog_distance_minus_sigma', 'catalog_distance_plus_sigma'],
                 static_scenarios={'R': 'each source row R_GC_kpc, not host distance',
                                   'vc_kms': [220.0, 198.0, 242.0],
                                   'tolerances': [0.01, 0.03, 0.1]},
                 interpretation='Conditional environmental sensitivity controls; no posterior interval, remnant path, escape classification or tidally re-equilibrated host.',
                 outer_root_policy='Spherical host gravity continues as the point-mass field beyond its finite DF support. Root brackets may extend beyond rt, avoiding a spurious failure for compact hosts.',
                 catalog_comparison_policy='The 5-percent peri/apo comparison is a descriptive model-consistency flag, independent of numerical integration acceptance; a mismatch is retained, never fitted away.',
                 host_gate='Require an exported host with all_numerical_checks_passed=true; absent or rejected hosts retain catalog-only orbits.',
                 diagnostic_retry='Halve output cadence and derivative step up to four times; compare successive sufficient-radius minima against unchanged 0.5-percent tolerance; record every level.',
                 integration_retry='On a failed numerical check, halve integration max steps and tighten relative/absolute tolerances by ten, at most three times; preserve preceding attempts. Acceptance limits never change.',
                 uncertainty='Quoted distance errors define sensitivity cases only. Full phase-space, Galactic-model and host/orbit cross-covariances are unavailable.')
    save_json(DATA / 'protocol.json', value)
    return value


def bracketed_root(function, lower, upper):
    lo, hi = float(function(lower)), float(function(upper))
    if not np.isfinite([lo, hi]).all() or lo * hi >= 0:
        raise ValueError(f'No root bracket: endpoints {lo:g}, {hi:g}')
    return brentq(function, lower, upper, xtol=1e-8, rtol=1e-12)


def static_diagnostics(host, R, vc, tolerances):
    """Same static formulas, with valid point-mass continuation outside rt."""
    lower, upper = 1e-3 * host.host.rh, .2 * R
    tau = vc ** 2 / R ** 2
    rj = bracketed_root(lambda r: host.G * host.enclosed(r) / r ** 3 - 2 * tau,
                        lower, upper)
    result = {'R_pc': R, 'vc_kms': vc, 'rJ_extended_pc': rj,
              'rJ_point_mass_pc': (host.G * host.M / (2 * tau)) ** (1 / 3),
              'rJ_over_rh': rj / host.host.rh,
              'mass_outside_rJ_fraction': 1 - float(host.enclosed(rj)) / host.M,
              'root_relative_residual': abs(host.G * host.enclosed(rj) / rj ** 3 / (2 * tau) - 1),
              'eta_at_rh': float(max_force_ratio(host.host.rh, host.enclosed, host.G, R, vc)),
              'eta_at_3rh': float(max_force_ratio(3 * host.host.rh, host.enclosed, host.G, R, vc)),
              'hill_central_opening_speed_kms': np.sqrt(2 * (host.potential(rj) - tau * rj * rj - host.potential(0)))}
    for side, sign in [('inward', -1), ('outward', 1)]:
        fn = lambda r: host.G * host.enclosed(r) / r ** 3 - vc ** 2 / (R * (R + sign * r)) - tau
        rs = bracketed_root(fn, lower, upper)
        result[f'exact_{side}_saddle_pc'] = rs
        result[f'exact_{side}_relative_residual'] = abs(fn(rs) / (2 * tau))
    for eps in tolerances:
        rv = bracketed_root(lambda r: max_force_ratio(r, host.enclosed, host.G, R, vc) - eps,
                            lower, upper)
        result[f'r_valid_{eps:g}_pc'] = rv
        result[f'valid_{eps:g}_relative_residual'] = abs(max_force_ratio(rv, host.enclosed, host.G, R, vc) / eps - 1)
    return {k: float(v) for k, v in result.items()}


def host_diagnostics(galaxy, host, states, derivative_step):
    """Registered diagnostics; brackets extend outside rt using exterior gravity."""
    x, v = states[:, :3], states[:, 3:]
    tensor = galaxy.fields(x)[2]
    norm = np.max(np.abs(np.linalg.eigvalsh(tensor)), axis=1) / 1e6
    omega, omega_dot = angular_fields(x, v, galaxy.fields(x)[1])
    effective = tensor + np.sum(omega ** 2, axis=1)[:, None, None] * np.eye(3) - omega[:, :, None] * omega[:, None, :]
    lam = np.linalg.eigvalsh(effective)[:, -1] / 1e6
    shift = K * v * derivative_step
    plus, minus = galaxy.fields(x + shift)[2], galaxy.fields(x - shift)[2]
    variation = np.linalg.norm((plus - minus) / (2 * derivative_step), ord=2, axis=(1, 2)) / 1e6 / norm
    amplitude = np.abs(np.linalg.norm(plus, ord=2, axis=(1, 2)) - np.linalg.norm(minus, ord=2, axis=(1, 2))) / (2 * derivative_step) / 1e6 / norm
    radii, scales = [], []
    for point, eigenvalue in zip(x, lam):
        lower, upper = 1e-3 * host.host.rh, .1 * np.linalg.norm(point) * 1000
        def ratio(r):
            return galaxy.hessian_ball_bound(point, r / 1000) / 1e6 * r ** 3 / (host.G * host.enclosed(r))
        radii.append(bracketed_root(lambda r: ratio(r) - .01, lower, upper))
        if eigenvalue <= 0:
            raise ValueError('No extensive frozen-frame eigenvalue')
        scales.append(bracketed_root(lambda r: host.G * host.enclosed(r) / r ** 3 - eigenvalue, lower, upper))
    radii, scales = np.asarray(radii), np.asarray(scales)
    def tdyn(r):
        return np.sqrt(r ** 3 / (host.G * host.enclosed(r))) / (K * 1000)
    rh = host.host.rh
    return {'tidal_norm_kms2_pc2': norm, 'tidal_variation_per_Myr': variation,
            'tidal_amplitude_variation_per_Myr': amplitude,
            'euler_over_extensive_eigenvalue': np.linalg.norm(omega_dot, axis=1) / (lam * 1e6),
            'r_sufficient_1pct_pc': radii, 'frozen_frame_scale_pc': scales,
            'mass_beyond_frozen_scale_fraction': 1 - host.enclosed(scales) / host.M,
            'eta_linear_rh': norm * rh ** 3 / (host.G * host.enclosed(rh)),
            'eta_linear_3rh': norm * (3 * rh) ** 3 / (host.G * host.enclosed(3 * rh)),
            'eta_bound_rh': galaxy.hessian_ball_bound(x, rh / 1000) / 1e6 * rh ** 3 / (host.G * host.enclosed(rh)),
            'eta_bound_3rh': galaxy.hessian_ball_bound(x, 3 * rh / 1000) / 1e6 * (3 * rh) ** 3 / (host.G * host.enclosed(3 * rh)),
            'timescale_ratio_rh': tdyn(rh) * variation,
            'amplitude_timescale_ratio_rh': tdyn(rh) * amplitude,
            'timescale_ratio_3rh': tdyn(3 * rh) * variation,
            'timescale_ratio_frozen_scale': tdyn(scales) * variation}


def run_case(source, name, distance, p, galaxy, out, refinement_level=0):
    cfg, limits = dict(p['integration']), p['acceptance']
    cfg['refinement_rtol']=max(cfg['refinement_rtol'],3e-14)
    times = np.arange(cfg['time_Myr'][0], cfg['time_Myr'][1] + cfg['output_step_Myr'] / 2, cfg['output_step_Myr'])
    y0 = initial_state(source, distance)
    sol = integrate(galaxy, y0, cfg['time_Myr'], cfg['max_step_Myr'], cfg['rtol'], cfg['atol'])
    states = sol.sol(times).T
    energy, lz = conserved(galaxy, states)
    case = {'case': name, 'distance_kpc': distance, 'initial_state_kpc_kms': y0.tolist(),
            'energy_relative_drift': float(np.max(np.abs((energy - energy[0]) / energy[0]))),
            'Lz_relative_drift': float(np.max(np.abs((lz - lz[0]) / lz[0]))),
            'maximum_abs_Z_kpc': float(np.max(np.abs(states[:, 2])))}
    events, extrema = [], {'pericenter': [], 'apocenter': []}
    for t, y in zip(sol.t_events[0], sol.y_events[0]):
        kind = 'pericenter' if np.sum(y[3:] ** 2) + y[:3] @ galaxy.fields(y[:3])[1] > 0 else 'apocenter'
        radius = float(np.linalg.norm(y[:3])); extrema[kind].append(radius)
        events.append({'event': kind, 't_Myr': float(t), 'R_kpc': radius, 'Z_kpc': float(y[2])})
    for t, y in zip(sol.t_events[1], sol.y_events[1]):
        events.append({'event': 'plane_up' if y[5] > 0 else 'plane_down', 't_Myr': float(t),
                       'R_kpc': float(np.linalg.norm(y[:3])), 'Z_kpc': float(y[2])})
    for label, values in extrema.items():
        case[label + '_mean_kpc'] = float(np.mean(values)) if values else None
        case[label + '_count'] = len(values)
    ref = integrate(galaxy, y0, cfg['time_Myr'], cfg['refinement_max_step_Myr'], cfg['refinement_rtol'], 1e-13)
    refined = ref.sol(times).T
    case['refinement_max_position_kpc'] = float(np.max(np.linalg.norm(refined[:, :3] - states[:, :3], axis=1)))
    case['refinement_max_velocity_kms'] = float(np.max(np.linalg.norm(refined[:, 3:] - states[:, 3:], axis=1)))
    reverse = integrate(galaxy, states[-1], cfg['time_Myr'][::-1], cfg['max_step_Myr'], cfg['rtol'], cfg['atol'])
    case['reversal_position_error_kpc'] = float(np.linalg.norm(reverse.y[:3, -1] - y0[:3]))
    keys = [('energy_relative_drift', 'relative_energy_drift'), ('Lz_relative_drift', 'relative_Lz_drift'),
            ('refinement_max_position_kpc', 'refined_position_kpc'), ('refinement_max_velocity_kms', 'refined_velocity_kms'),
            ('reversal_position_error_kpc', 'time_reversal_position_kpc')]
    case['all_orbit_numerical_checks_passed'] = all(case[k] < limits[lim] for k, lim in keys)
    case['integration_settings'] = dict(cfg)
    case['refinement_level'] = refinement_level
    if not case['all_orbit_numerical_checks_passed'] and refinement_level < 3:
        tightened = json.loads(json.dumps(p))
        for key in ['max_step_Myr', 'refinement_max_step_Myr']:
            tightened['integration'][key] /= 2
        for key in ['rtol', 'atol', 'refinement_rtol']:
            tightened['integration'][key] /= 10
        accepted, accepted_states, accepted_ref = run_case(source, name, distance, tightened, galaxy, out, refinement_level+1)
        accepted.setdefault('preceding_numerical_attempts', []).insert(0, case)
        save_json(out / f'{name}.json', accepted)
        return accepted, accepted_states, accepted_ref
    np.savez_compressed(out / f'{name}-track.npz', t_Myr=times, state_kpc_kms=states)
    save_json(out / f'{name}-events.json', events)
    save_json(out / f'{name}.json', case)
    return case, states, ref


def load_host(cluster, p):
    path = ROOT / p['host_path_pattern'].format(cluster=cluster)
    summary = path.with_name('summary.json')
    if not path.exists() or not summary.exists():
        return None, None, 'missing_reconstructed_host'
    hs = json.loads(summary.read_text())
    if hs.get('all_numerical_checks_passed') is not True:
        return None, hs, 'host_numerical_checks_not_passed'
    return EnclosedHost(Host(dict(np.load(path, allow_pickle=False)))), hs, 'available_conditional_host'


def run_cluster(cluster, source, p, catalog_only=False):
    out = OUT / cluster; out.mkdir(parents=True, exist_ok=True)
    galaxy, cfg = Galaxy(p['galaxy']), p['integration']
    host, hs, host_status = load_host(cluster, p)
    if catalog_only:
        host, host_status = None, 'catalog_only_requested'
    distances = [('catalog_reference', source['distance_kpc']),
                 ('catalog_distance_minus_sigma', source['distance_kpc'] - source['distance_error_kpc']),
                 ('catalog_distance_plus_sigma', source['distance_kpc'] + source['distance_error_kpc'])]
    if host is not None:
        distances.insert(1, ('host_distance', hs['summary']['distance_kpc']))
    records, numerical_errors = [], []
    for name, distance in distances:
        print(f'{cluster}: {name}', flush=True)
        case, states, refined = run_case(source, name, distance, p, galaxy, out)
        if not case['all_orbit_numerical_checks_passed']:
            numerical_errors.append(name + ': orbit convergence failed')
        if host is not None:
            try:
                tidal = host_diagnostics(galaxy, host, states, p['tensor_derivative_step_Myr'])
                for key, values in tidal.items():
                    case[key + '_min'], case[key + '_max'] = float(np.min(values)), float(np.max(values))
                np.savez_compressed(out / f'{name}-environment.npz', **tidal)
                case['host_diagnostics_status'] = 'calculated'
                if name == 'host_distance':
                    previous_radius=case['r_sufficient_1pct_pc_min']
                    history=[]
                    for level in range(1,5):
                        step=cfg['output_step_Myr']/2**level
                        times=np.arange(cfg['time_Myr'][0],cfg['time_Myr'][1]+step/2,step)
                        fine=host_diagnostics(galaxy,host,refined.sol(times).T,p['tensor_derivative_step_Myr']/2**level)
                        min_radius=float(np.min(fine['r_sufficient_1pct_pc']))
                        difference=abs(min_radius/previous_radius-1)
                        history.append({'output_step_Myr':step,'minimum_sufficient_radius_pc':min_radius,'relative_change':difference})
                        previous_radius=min_radius
                        if difference < p['acceptance']['bound_radius_relative_refinement']:break
                    case['diagnostic_refinement_history']=history
                    case['refined_minimum_sufficient_radius_pc']=min_radius
                    case['sufficient_radius_cadence_relative_difference']=difference
                    case['proposed_fixed_isolated_control_pc']=float(np.floor(min_radius))
                    case['timescale_ratio_rh_refinement_absolute_difference']=abs(float(np.max(fine['timescale_ratio_rh']))-case['timescale_ratio_rh_max'])
                    case['refined_timescale_ratio_rh_max']=float(np.max(fine['timescale_ratio_rh']))
                    if difference >= p['acceptance']['bound_radius_relative_refinement']:
                        numerical_errors.append(name + ': diagnostic cadence convergence failed')
            except ValueError as exc:
                case['host_diagnostics_status'] = 'not_supported_in_declared_domain'
                case['host_diagnostics_reason'] = str(exc)
                numerical_errors.append(name + ': ' + str(exc))
        save_json(out / f'{name}.json', case)
        records.append(case)
    catalog = records[0]
    comparison = {label + '_fractional_difference': abs(catalog[label + '_mean_kpc'] / source['catalog_mean_' + label + '_kpc'] - 1)
                  for label in ['pericenter', 'apocenter'] if catalog[label + '_mean_kpc'] is not None}
    comparison['within_5_percent'] = len(comparison) == 2 and max(comparison.values()) < .05
    static = []
    if host is not None:
        for vc in p['static_scenarios']['vc_kms']:
            try:
                row = static_diagnostics(host, 1000 * source['R_GC_kpc'], vc, p['static_scenarios']['tolerances'])
                coarse = static_diagnostics(EnclosedHost(host.host, stride=2), 1000 * source['R_GC_kpc'], vc, p['static_scenarios']['tolerances'])
                keys = ['rJ_extended_pc'] + [f'r_valid_{eps:g}_pc' for eps in p['static_scenarios']['tolerances']]
                row['max_interpolation_refinement_relative'] = max(abs(coarse[k] / row[k] - 1) for k in keys)
                errors = []
                for key in keys:
                    radius = min(row[key], host.rt)
                    mass = quad(lambda r: 4 * np.pi * r * r * host.host.rho_spline(r)[host.host.physical].sum(), 0, radius,
                                epsabs=1e-5, epsrel=1e-8, limit=500)[0]
                    errors.append(abs(mass / host.enclosed(radius) - 1))
                row['max_density_integral_relative_error'] = max(errors)
                row['all_static_numerical_checks_passed'] = (row['max_interpolation_refinement_relative'] < .0002
                    and row['max_density_integral_relative_error'] < .0002
                    and max(v for k, v in row.items() if 'residual' in k) < 1e-9)
                if not row['all_static_numerical_checks_passed']:
                    numerical_errors.append(f'static vc={vc}: numerical checks failed')
                static.append(row)
            except ValueError as exc:
                static.append({'vc_kms': vc, 'status': 'not_supported', 'reason': str(exc)})
                numerical_errors.append(f'static vc={vc}: {exc}')
    result = {'revision': p['revision'], 'cluster': cluster, 'host_status': host_status,
              'galactic_orbits_integrated': True, 'remnant_trajectories_performed': False,
              'host_re_equilibrated': False, 'cases': records, 'static_controls': static,
              'catalog_orbit_comparison': comparison, 'numerical_errors': numerical_errors,
              'all_numerical_checks_passed': not numerical_errors,
              'source_sha256': digest(DATA / f'{cluster}.json'), 'protocol_sha256': digest(DATA / 'protocol.json'),
              'host_sha256': digest(ROOT / p['host_path_pattern'].format(cluster=cluster)) if host is not None else None,
              'environment': {'python': sys.version.split()[0], 'numpy': np.__version__, 'scipy': scipy.__version__}}
    save_json(out / 'summary.json', result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--clusters', nargs='+', choices=CLUSTERS, default=CLUSTERS)
    parser.add_argument('--catalog-only', action='store_true')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    DATA.mkdir(parents=True, exist_ok=True); OUT.mkdir(parents=True, exist_ok=True)
    sources, p = prepare_sources(), protocol()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures=[pool.submit(run_cluster,cluster,sources[cluster],p,args.catalog_only) for cluster in args.clusters]
        for future in futures:future.result()
    summaries = [json.loads((OUT / c / 'summary.json').read_text()) for c in CLUSTERS if (OUT / c / 'summary.json').exists()]
    save_json(OUT / 'summary.json', {'revision': p['revision'], 'clusters': summaries,
                                  'all_six_clusters_present': len(summaries) == 6,
                                  'remnant_trajectories_performed': False})
    record_provenance()


def record_provenance():
    paths = [Path(__file__), ROOT / 'analysis/orbit/galaxy.py', ROOT / 'analysis/orbit/orbits.py',
             ROOT / 'analysis/host_reconstruction/host.py', ROOT / 'analysis/boundary/model.py']
    paths += sorted(DATA.glob('*'))
    save_json(OUT / 'provenance.json', {'input_sha256': {str(x.relative_to(ROOT)): digest(x) for x in paths if x.is_file()},
              'output_sha256': {str(x.relative_to(ROOT)): digest(x) for x in sorted(OUT.rglob('*')) if x.is_file() and x.name != 'provenance.json'}})


if __name__ == '__main__':
    main()
