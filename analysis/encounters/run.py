"""Finite-angle local BH rates, matched moments and sampler checks; no orbits."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'analysis/host_reconstruction'))
from host import Host
from physics import (DiscreteBath, PC_PER_KMS_MYR, closest_approach_annulus,
                     impact_moments, velocity_ball)

OUT = ROOT / 'results/finite-angle-encounters'
PROTOCOL = Path(__file__).with_name('protocol.json')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def density_scale(host, selected):
    rho0 = float(host.rho_spline(0)[selected].sum())
    if rho0 <= 0:
        raise ValueError('A BH density scale is undefined for an absent population')
    return brentq(lambda r: float(host.rho_spline(r)[selected].sum())-rho0/2,
                  0, host.rt, xtol=1e-13)


def local_scale(host, selected, radius, core):
    rho = float(host.rho_spline(radius)[selected].sum())
    if radius == 0:
        return core
    h = min(core*1e-4, radius/4)
    derivative = float((host.rho_spline(radius+h)[selected].sum()
                        - host.rho_spline(radius-h)[selected].sum())/(2*h))
    return min(core, rho/abs(derivative)) if derivative else core


def elastic_floor(speed, mass, remnant, G, bmax, prescription):
    c = prescription['speed_of_light_kms']
    total = mass + remnant
    eta = mass*remnant/total**2
    # The zero-speed incoming flux vanishes. Avoid infinite helper inputs.
    safe_speed = np.maximum(speed, np.finfo(float).tiny**.25)
    rg = G*total/c**2
    rcap = rg*((85*np.pi/(6*np.sqrt(2)))*eta*(c/safe_speed)**2)**(2/7)
    floor = np.maximum(prescription['weak_field_radius_GM_over_c2']*rg,
                       prescription['parabolic_capture_radius_multiplier']*rcap)
    return closest_approach_annulus(speed, mass, remnant, G, floor, bmax)[0]


def large_deflection_scale_trace(host, selected, radius, velocity, length, bmax, p, order):
    """Integrate b90>L on its exact relative-velocity ball per species.

    Clipping a coarse quadrature by an indicator produces oscillatory errors.
    Centering a separate ball at the test velocity resolves this boundary as
    an integration endpoint instead; the DF remains in the cluster frame.
    """
    total = 0.
    for index in np.flatnonzero(selected):
        mass = host.mj[index]
        cutoff = np.sqrt(host.G*(mass+p['remnant_mass_Msun'])/length)
        w, weights = velocity_ball(cutoff, order, radial_scale=np.sqrt(host.s2j[index]))
        df = host.number_df(radius, velocity[None, :]-w)[index]
        speed = np.linalg.norm(w, axis=1)
        lower = elastic_floor(speed, mass, p['remnant_mass_Msun'], host.G,
                              bmax, p['closest_approach'])
        moments = impact_moments(w, mass, p['remnant_mass_Msun'], host.G, lower, bmax)
        total += float(np.sum(df*weights*np.trace(moments['D'],axis1=-2,axis2=-1))/PC_PER_KMS_MYR)
    return total


def evaluate(host, selected, radius, velocity, length, eta, p, order):
    sigma_scale = float(np.sqrt(host.s2j[selected]).min())
    nodes, weights = velocity_ball(float(np.sqrt(2*host.psi(radius))), order,
                                    radial_scale=sigma_scale)
    numbers = host.number_df(radius, nodes)[selected]*weights[None, :]
    masses = host.mj[selected]
    bmax = eta*length
    lower = lambda speed, mass: elastic_floor(speed, mass, p['remnant_mass_Msun'],
                                               host.G, bmax, p['closest_approach'])
    bath = DiscreteBath(nodes, numbers, masses, p['remnant_mass_Msun'], host.G,
                        lower, bmax)
    ans = bath.moments(velocity)
    relative = velocity[None, :] - nodes
    speed = np.linalg.norm(relative, axis=1)[None, :]
    mass = masses[:, None]
    lo = lower(speed, mass)
    # The diagnostic is not an established error criterion for residence.
    large_trace = large_deflection_scale_trace(host, selected, radius, velocity,
                                               length, bmax, p, order)
    trace = float(np.trace(ans['D']))
    parent_rate = float(np.sum(numbers*speed*np.pi*bmax*bmax)/PC_PER_KMS_MYR)
    excluded_rate = float(np.sum(numbers*speed*np.pi*lo**2)/PC_PER_KMS_MYR)
    exact_number = float(np.sum(host.rho_spline(radius)[selected]/masses))
    ans.update(trace=trace, excluded_rate=excluded_rate, parent_rate=parent_rate,
               nonlocal_trace_fraction=large_trace/trace,
               number_density_relative_error=abs(float(numbers.sum())/exact_number-1),
               bmax_pc=bmax)
    return ans, bath


def difference(old, new, sigma):
    scale_a = max(np.linalg.norm(new['A']), new['trace']/sigma, np.finfo(float).tiny)
    return {
        'A_refinement': float(np.linalg.norm(new['A']-old['A'])/scale_a),
        'D_refinement': float(np.linalg.norm(new['D']-old['D'])/np.linalg.norm(new['D'])),
        'rate_refinement': abs(new['rate']/old['rate']-1),
        'excluded_rate_refinement': abs(new['excluded_rate']/old['excluded_rate']-1),
        'diagnostic_fraction_refinement': abs(new['nonlocal_trace_fraction']-old['nonlocal_trace_fraction']),
    }


def accepted(ans, diff, velocity, sigma, p):
    a = p['acceptance']
    drift_zero = np.linalg.norm(velocity) != 0 or np.linalg.norm(ans['A']) <= a['zero_drift_relative_scale']*ans['trace']/sigma
    return bool(max(diff[k] for k in ('A_refinement','D_refinement','rate_refinement','excluded_rate_refinement')) < a['moment_relative_refinement']
                and diff['diagnostic_fraction_refinement'] < a['diagnostic_fraction_absolute_refinement']
                and ans['number_density_relative_error'] < a['number_density_relative_error']
                and np.linalg.eigvalsh(ans['D'])[0] >= 0 and drift_zero)


def sample_check(bath, ans, velocity, n, seed):
    increments = bath.sample_events(velocity, np.random.default_rng(seed), n)
    measured = np.column_stack((increments, increments**2,
                                increments[:, 0]*increments[:, 1],
                                increments[:, 0]*increments[:, 2],
                                increments[:, 1]*increments[:, 2]))
    expected = np.r_[ans['A'], np.diag(ans['D']), ans['D'][0,1], ans['D'][0,2], ans['D'][1,2]]/ans['rate']
    mean = measured.mean(axis=0)
    se = measured.std(axis=0, ddof=1)/np.sqrt(n)
    if np.any((se == 0) & ~np.isclose(mean, expected, rtol=1e-12, atol=1e-15)):
        raise AssertionError('Degenerate sampler differs from its exact moment')
    z = np.divide(abs(mean-expected), se, out=np.zeros_like(se), where=se>0)
    return {'events': n, 'expected_conditional_moments': expected.tolist(),
            'sampled_conditional_moments': mean.tolist(), 'standard_errors': se.tolist(),
            'standardized_errors': z.tolist(), 'maximum_standardized_error': float(z.max()),
            'columns': ['dx','dy','dz','dx2','dy2','dz2','dxdy','dxdz','dydz']}


def run(clusters=None):
    p = json.loads(PROTOCOL.read_text())
    OUT.mkdir(parents=True, exist_ok=True)
    for cluster in clusters or p['clusters']:
        if cluster not in p['clusters']:
            raise ValueError('Cluster is not in the declared protocol')
        source = ROOT / f'results/cluster-comparison/{cluster}/reference/host.npz'
        host_record = json.loads(source.with_name('summary.json').read_text())
        if not host_record['all_numerical_checks_passed']:
            raise ValueError('Reference host has not passed its numerical checks')
        host = Host(dict(np.load(source, allow_pickle=False)))
        selected = host.types == 'BH'
        core = density_scale(host, selected)
        rows, tensors, sampling = [], [], None
        for fraction in p['r_over_density_radius']:
            radius = core*fraction
            density = host.rho_spline(radius)[selected]/host.mj[selected]
            v2 = (host.pr_spline(radius)[selected]+host.pt_spline(radius)[selected])/host.mj[selected]
            sigma = float(np.sqrt(v2.sum()/(3*density.sum())))
            length = local_scale(host, selected, radius, core)
            for speed in p['v_over_BH_sigma_1d']:
                for direction in (['radial'] if speed == 0 else p['orientations']):
                    axis = np.array([0.,0.,1.]) if direction == 'radial' else np.array([1.,0.,0.])
                    velocity = speed*sigma*axis
                    for eta in p['bmax_over_local_scale']:
                        previous, history = None, []
                        for order in p['quadrature_orders']:
                            answer, bath = evaluate(host, selected, radius, velocity, length, eta, p, order)
                            if previous is not None:
                                refinement = difference(previous, answer, sigma)
                                history.append({'coarse_order': prior_order, 'fine_order': order, **refinement})
                                if accepted(answer, refinement, velocity, sigma, p):
                                    break
                            previous, prior_order = answer, order
                        passed = accepted(answer, refinement, velocity, sigma, p)
                        row = dict(cluster=cluster, state=len(rows), r_over_density_radius=fraction,
                                   radius_pc=radius, density_radius_pc=core, local_scale_pc=length,
                                   sigma_BH_1d_kms=sigma, speed_over_sigma=speed, orientation=direction,
                                   bmax_over_local_scale=eta, bmax_pc=answer['bmax_pc'],
                                   allowed_rate_per_Myr=answer['rate'], excluded_rate_per_Myr=answer['excluded_rate'],
                                   parent_rate_per_Myr=answer['parent_rate'],
                                   excluded_flux_fraction=answer['excluded_rate']/answer['parent_rate'],
                                   trace_D_kms2_per_Myr=answer['trace'],
                                   drift_along_velocity_kms_per_Myr=float(answer['A']@axis),
                                   rms_increment_over_sigma=float(np.sqrt(answer['trace']/answer['rate'])/sigma),
                                   rate_times_local_crossing=answer['rate']*length/sigma*PC_PER_KMS_MYR,
                                   b90_above_L_trace_fraction=answer['nonlocal_trace_fraction'],
                                   number_density_relative_error=answer['number_density_relative_error'],
                                   minimum_D_eigenvalue=float(np.linalg.eigvalsh(answer['D'])[0]),
                                   quadrature_order=order, accepted=passed, **refinement)
                        rows.append(row)
                        tensors.append({'state': row['state'], 'A': answer['A'].tolist(), 'D': answer['D'].tolist(),
                                        'refinement_history': history})
                        if fraction == 1 and speed == 1 and direction == 'radial' and eta == .1:
                            sampling = sample_check(bath, answer, velocity, p['sampling']['events_per_host'],
                                                    np.random.SeedSequence([p['sampling']['seed'],p['clusters'].index(cluster)]))
                        print(cluster, row['state'], order, 'accepted' if passed else 'UNRESOLVED', flush=True)
        result = {'cluster': cluster, 'density_radius_pc': core, 'states': rows, 'tensors': tensors,
                  'sampling': sampling, 'all_numerical_checks_passed': all(r['accepted'] for r in rows)
                  and sampling['maximum_standardized_error'] < p['sampling']['standard_error_limit'],
                  'protocol_sha256': sha(PROTOCOL), 'host_sha256': sha(source),
                  'implementation_sha256': {f.name: sha(f) for f in [Path(__file__),Path(__file__).with_name('physics.py')]}}
        (OUT/f'{cluster}.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    if all((OUT/f'{c}.json').exists() for c in p['clusters']):
        results = [json.loads((OUT/f'{c}.json').read_text()) for c in p['clusters']]
        if any(r['protocol_sha256'] != sha(PROTOCOL) for r in results):
            raise ValueError('Stale per-host protocol; rerun all hosts')
        current_code = {f.name: sha(f) for f in [Path(__file__),Path(__file__).with_name('physics.py')]}
        for result in results:
            source = ROOT/f"results/cluster-comparison/{result['cluster']}/reference/host.npz"
            if result['implementation_sha256'] != current_code or result['host_sha256'] != sha(source):
                raise ValueError('Mixed or stale host/implementation; rerun affected hosts')
        rows = [s for r in results for s in r['states']]
        with (OUT/'states.csv').open('w',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
        summary = {'revision':p['revision'], 'clusters':len(results), 'states':len(rows),
                   'accepted_states':sum(r['accepted'] for r in rows),
                   'all_numerical_checks_passed':all(r['all_numerical_checks_passed'] for r in results),
                   'sampled_elastic_events':sum(r['sampling']['events'] for r in results),
                   'max_sample_standardized_error':max(r['sampling']['maximum_standardized_error'] for r in results),
                   'transport_simulations_performed':False, 'residence_response_evaluated':False,
                   'scope':p['interpretation']}
        (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--cluster',action='append')
    run(parser.parse_args().cluster)
