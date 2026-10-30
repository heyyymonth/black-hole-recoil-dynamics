"""Collisionless remnant controls and path-integrated BH encounter diagnostics."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'analysis/encounters'),
               str(ROOT/'analysis/host_reconstruction')]
from host import Host
from physics import DiscreteBath, PC_PER_KMS_MYR, velocity_ball
from analysis.encounters.run import density_scale, local_scale, elastic_floor
from analysis.residence.orbits import integrate_orbit, host_acceleration
from analysis.residence.radial_reference import radial_residence

OUT = ROOT/'results/collisionless-residence'
PROTOCOL = Path(__file__).with_name('protocol.json')
CODE = [Path(__file__), Path(__file__).with_name('orbits.py'),
        Path(__file__).with_name('radial_reference.py'),
        ROOT/'analysis/host_reconstruction/host.py',
        ROOT/'analysis/encounters/physics.py', ROOT/'analysis/encounters/run.py']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def local_budget(host, selected, radius, velocity, core, p, ep):
    """Rotate to the spherical DF frame; refine the same finite encounter law."""
    density = host.rho_spline(radius)[selected]
    exact_n = float(np.sum(density/host.mj[selected]))
    bmax = p['bmax_over_local_scale']*local_scale(host, selected, radius, core)
    sigma = float(np.sqrt(np.sum((host.pr_spline(radius)[selected]
                                 + host.pt_spline(radius)[selected])/host.mj[selected])
                          / (3*exact_n)))
    lower = lambda speed, mass: elastic_floor(speed, mass, p['remnant_mass_Msun'],
                                              host.G, bmax, ep['closest_approach'])
    previous, refinements = None, []
    for order in p['velocity_quadrature_orders']:
        # Center on relative velocity to resolve the integrable close-passage
        # flux near w=0. This ball contains the entire shifted DF support;
        # outside its original field-velocity support the DF is exactly zero.
        relative_nodes, weights = velocity_ball(np.sqrt(2*host.psi(radius))
                                                +np.linalg.norm(velocity), order,
                                                radial_scale=np.sqrt(host.s2j[selected]).min())
        nodes = velocity-relative_nodes
        numbers = host.number_df(radius, nodes)[selected]*weights
        bath = DiscreteBath(nodes, numbers, host.mj[selected], p['remnant_mass_Msun'],
                            host.G, lower, bmax)
        moment = bath.moments(velocity)
        relative = np.linalg.norm(velocity-nodes, axis=1)[None, :]
        lo = lower(relative, host.mj[selected, None])
        excluded = float(np.sum(numbers*relative*np.pi*lo**2)/PC_PER_KMS_MYR)
        parent = float(np.sum(numbers*relative*np.pi*bmax**2)/PC_PER_KMS_MYR)
        trace = float(np.trace(moment['D']))
        values = np.array([moment['rate'], excluded, np.linalg.norm(moment['A']), trace])
        number_error = abs(float(numbers.sum())/exact_n-1)
        if previous is not None:
            old_values, old_moment = previous
            # Same scaled drift convention as the local encounter calculation.
            changes = {
                'rate': abs(values[0]/old_values[0]-1),
                'excluded': abs(values[1]/old_values[1]-1),
                'A': float(np.linalg.norm(moment['A']-old_moment['A'])
                           / max(np.linalg.norm(moment['A']),trace/sigma)),
                'D': float(np.linalg.norm(moment['D']-old_moment['D'])
                           / np.linalg.norm(moment['D']))}
            refinements.append({'order':order, **changes})
            accepted = (max(changes.values()) < p['acceptance']['local_moment_relative_refinement']
                        and number_error < p['acceptance']['local_number_density_relative_error']
                        and np.linalg.eigvalsh(moment['D'])[0] >= 0)
            if accepted:
                break
        previous = values, moment
    return values, {'order':order, 'accepted':bool(accepted),
                    'number_density_relative_error':number_error,
                    'parent_rate_per_Myr':parent,
                    'rate_partition_relative_error':abs((values[0]+values[1])/parent-1),
                    'minimum_D_eigenvalue':float(np.linalg.eigvalsh(moment['D'])[0]),
                    'refinement_history':refinements}


def integrate_budgets(host, selected, orbit, core, vdyn, tdyn, launch_speed, p, ep):
    boundaries = sorted(set([0., orbit.end_time, *orbit.turning_times,
                             *(e['time'] for e in orbit.crossings)]))
    boundaries = [t for t in boundaries if 0 <= t <= orbit.end_time]
    rho0 = float(host.rho_spline(0)[selected].sum())
    previous, history, node_records = None, [], []
    for order in p['time_quadrature_orders_per_phase']:
        nodes, weights = leggauss(order)
        total = np.zeros(5)
        records = []
        for left,right in zip(boundaries[:-1],boundaries[1:]):
            if right-left < 1e-12:
                continue
            times = left+(nodes+1)*(right-left)/2
            weights_t = weights*(right-left)/2*tdyn
            for time,weight in zip(times,weights_t):
                state = orbit.solution(float(time))
                x, v = state[:3]*core, state[3:]*vdyn
                radius = float(np.linalg.norm(x))
                radial = float(x@v/radius) if radius else float(np.linalg.norm(v))
                tangential = float(np.linalg.norm(np.cross(x,v))/radius) if radius else 0.
                values, record = local_budget(host,selected,radius,
                                              np.array([tangential,0.,radial]),core,p,ep)
                rho = float(host.rho_spline(radius)[selected].sum())
                total += weight*np.r_[values,rho/rho0]
                records.append({'time_over_tdyn':float(time), 'weight_Myr':float(weight),
                                'radius_over_density_radius':radius/core,
                                'local_integrands':values.tolist(),
                                'rho_BH_over_central':rho/rho0, **record})
        normalized = total/np.array([1.,1.,launch_speed,launch_speed**2,
                                      p['horizon_over_t_dyn']*tdyn])
        node_records.append({'time_order':order,'nodes':records})
        if previous is not None:
            delta = abs(normalized-previous)
            relative = delta[:4]/np.maximum(abs(normalized[:4]),np.finfo(float).tiny)
            accepted = (max(relative) < p['acceptance']['budget_time_relative_refinement']
                        and delta[4] < p['acceptance']['density_exposure_time_absolute_refinement']
                        and all(x['accepted'] for x in records))
            history.append({'time_order':order, 'budget_relative_refinement':relative.tolist(),
                            'density_exposure_absolute_refinement':float(delta[4]),
                            'accepted':bool(accepted)})
            if accepted:
                break
        previous = normalized
    keys = ['allowed_rate_integral','excluded_flux_integral','drift_norm_integral_over_launch_speed',
            'diffusion_trace_integral_over_launch_speed_squared','normalized_BH_density_exposure']
    return dict(zip(keys,normalized.tolist())), history, node_records


def run_cluster(cluster, p, ep):
    source = ROOT/f'results/cluster-comparison/{cluster}/reference/host.npz'
    if not json.loads(source.with_name('summary.json').read_text())['all_numerical_checks_passed']:
        raise ValueError('Reference host must pass its checks')
    host = Host(dict(np.load(source,allow_pickle=False)))
    selected = host.types == 'BH'
    core = density_scale(host,selected)
    vdyn = float(np.sqrt(-core*host.psi_spline(core,1)))
    tdyn = core/vdyn*PC_PER_KMS_MYR
    r0 = p['launch']['radius_over_density_radius']*core
    apo = p['launch']['radial_apocenter_over_density_radius']*core
    speed = float(np.sqrt(2*(host.psi(r0)-host.psi(apo))))
    environment_file = ROOT/'results/cluster-comparison-environment/comparison.csv'
    with environment_file.open() as f:
        env = next(r for r in csv.DictReader(f) if r['cluster']==cluster)
    assert env['all_numerical_checks_passed']=='True'
    stop = min(host.rt,float(env['refined_minimum_1pct_radius_pc']))
    if apo >= stop:
        raise ValueError('Diagnostic launch exceeds the declared control boundary')
    acceleration = host_acceleration(host,core,vdyn)
    records = []
    for direction in p['launch']['orientations']:
        x0 = np.array([p['launch']['radius_over_density_radius'],0.,0.])
        axis = np.array([1.,0.,0.]) if direction=='radial_outward' else np.array([0.,1.,0.])
        v0 = speed/vdyn*axis
        orbits = [integrate_orbit(acceleration,x0,v0,p['horizon_over_t_dyn'],
                                 region_radius=1.,stop_radius=stop/core,**controls)
                  for controls in p['orbital_refinement']]
        orbit = orbits[-1]
        times = np.linspace(0,orbit.end_time,p['track_samples'])
        states = orbit.solution(times).T
        radii = np.linalg.norm(states[:,:3],axis=1)
        energy = np.sum(states[:,3:]**2,axis=1)/2-host.psi(core*radii)/vdyn**2
        angular = np.cross(states[:,:3],states[:,3:])
        orbit_checks = {
            'energy_drift_over_vdyn_squared':float(np.max(abs(energy-energy[0]))),
            'angular_momentum_drift_over_rho_vdyn':float(np.max(np.linalg.norm(angular-angular[0],axis=1))),
            'orbit_state_refinement_dimensionless':float(np.max(np.linalg.norm(states-orbits[0].solution(times).T,axis=1))),
            'residence_refinement_over_horizon':abs(orbit.residence-orbits[0].residence)/p['horizon_over_t_dyn']}
        psi = lambda r: float(host.psi(core*r))/vdyn**2
        independent = [radial_residence(psi,x0,v0,orbit.end_time,region_radius=1.,
                                        outer_radius=3.,order=order)
                       for order in p['independent_radial_orders']]
        orbit_checks['independent_residence_difference_over_horizon'] = max(
            abs(orbit.residence-ref['residence'])/p['horizon_over_t_dyn'] for ref in independent)
        values, time_refinement, local_checks = integrate_budgets(
            host,selected,orbit,core,vdyn,tdyn,speed,p,ep)
        passed = (all(value < p['acceptance'][key] for key,value in orbit_checks.items())
                  and time_refinement[-1]['accepted'])
        row = {'cluster':cluster,'orientation':direction,'density_radius_pc':core,
               't_dyn_Myr':tdyn,'launch_speed_kms':speed,
               'horizon_Myr':p['horizon_over_t_dyn']*tdyn,'stop_radius_pc':stop,
               'end_time_Myr':orbit.end_time*tdyn,'stop_cause':orbit.stop_cause,
               'residence_fraction':orbit.residence/p['horizon_over_t_dyn'],
               'crossing_count':len(orbit.crossings),
               **values, 'accepted':bool(passed)}
        records.append({'summary':row,'orbit_checks':orbit_checks,
                        'crossings':orbit.crossings,'turning_times':orbit.turning_times,
                        'independent_radial_references':independent,
                        'time_refinement':time_refinement,'local_checks':local_checks,
                        'track_columns':['t_over_tdyn','x_over_rho','y_over_rho','z_over_rho',
                                         'vx_over_vdyn','vy_over_vdyn','vz_over_vdyn'],
                        'track':np.column_stack([times,states]).tolist()})
        print(cluster,direction,'accepted' if passed else 'UNRESOLVED',json.dumps(row),flush=True)
    result = {'cluster':cluster,'records':records,'all_numerical_checks_passed':all(r['summary']['accepted'] for r in records),
              'protocol_sha256':digest(PROTOCOL),'host_sha256':digest(source),
              'encounter_protocol_sha256':digest(ROOT/p['encounter_protocol']),
              'environment_sha256':digest(environment_file),
              'implementation_sha256':{str(f.relative_to(ROOT)):digest(f) for f in CODE}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/f'{cluster}.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')


def aggregate(p):
    results = [json.loads((OUT/f'{c}.json').read_text()) for c in p['clusters']]
    for result in results:
        assert result['protocol_sha256']==digest(PROTOCOL)
        assert result['encounter_protocol_sha256']==digest(ROOT/p['encounter_protocol'])
        assert result['environment_sha256']==digest(ROOT/'results/cluster-comparison-environment/comparison.csv')
        assert result['host_sha256']==digest(ROOT/f"results/cluster-comparison/{result['cluster']}/reference/host.npz")
        assert result['implementation_sha256']=={str(f.relative_to(ROOT)):digest(f) for f in CODE}
    rows = [record['summary'] for result in results for record in result['records']]
    with (OUT/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    summary={'revision':p['revision'],'scope':p['scope'],'trajectories':len(rows),
             'accepted_trajectories':sum(row['accepted'] for row in rows),
             'all_numerical_checks_passed':all(row['accepted'] for row in rows),
             'collisionless_trajectories_performed':True,'collisional_trajectories_performed':False,
             'residence_response_evaluated':False,'interpretation':p['interpretation']}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cluster',action='append');parser.add_argument('--aggregate',action='store_true')
    args=parser.parse_args();protocol=json.loads(PROTOCOL.read_text())
    encounter_protocol=json.loads((ROOT/protocol['encounter_protocol']).read_text())
    if not args.aggregate:
        for cluster in args.cluster or protocol['clusters']:
            if cluster not in protocol['clusters']:
                raise ValueError('Cluster not declared in protocol')
            run_cluster(cluster,protocol,encounter_protocol)
    if all((OUT/f'{c}.json').exists() for c in protocol['clusters']):
        aggregate(protocol)
