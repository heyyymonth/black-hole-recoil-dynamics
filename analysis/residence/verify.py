"""Verify recorded collisionless controls without integrating trajectories.

Only the standard library and NumPy are needed. This checks recorded numerical
identities, thresholds and hashes; it does not independently rerun the host,
ODE, velocity integrals or the physical model.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json

import numpy as np
from numpy.polynomial.legendre import leggauss


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/collisionless-residence'
PROTOCOL = ROOT/'analysis/residence/protocol.json'
TIME_CONVERSION = 0.9777922216807892
BUDGET_KEYS = (
    'allowed_rate_integral', 'excluded_flux_integral',
    'drift_norm_integral_over_launch_speed',
    'diffusion_trace_integral_over_launch_speed_squared',
    'normalized_BH_density_exposure',
)
ORBIT_CHECK_KEYS = (
    'energy_drift_over_vdyn_squared',
    'angular_momentum_drift_over_rho_vdyn',
    'orbit_state_refinement_dimensionless',
    'residence_refinement_over_horizon',
    'independent_residence_difference_over_horizon',
)
IMPLEMENTATION = (
    'analysis/residence/run.py', 'analysis/residence/orbits.py',
    'analysis/residence/radial_reference.py',
    'analysis/host_reconstruction/host.py',
    'analysis/encounters/physics.py', 'analysis/encounters/run.py',
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_tree(value):
    if isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)
    elif isinstance(value, float):
        require(np.isfinite(value), 'nonfinite number in recorded evidence')


def close(actual, expected, name, *, rtol=2e-12, atol=2e-13):
    require(np.allclose(actual, expected, rtol=rtol, atol=atol), name)


def check_local_node(node, protocol):
    a = protocol['acceptance']
    values = np.asarray(node['local_integrands'], dtype=float)
    require(values.shape == (4,) and np.all(values >= 0), 'invalid local budget values')
    require(0 <= node['rho_BH_over_central'] <= 1+1e-6, 'invalid normalized BH density')
    require(node['radius_over_density_radius'] >= 0, 'invalid local radius')
    require(node['number_density_relative_error'] >= 0, 'negative density error')
    parent = node['parent_rate_per_Myr']
    require(parent > 0, 'positive parent encounter rate required')
    partition = abs((values[0]+values[1])/parent-1)
    close(node['rate_partition_relative_error'], partition, 'local rate partition identity')
    require(partition < 1e-12, 'local incoming flux partition does not close')
    orders = protocol['velocity_quadrature_orders']
    require(node['order'] in orders[1:], 'unregistered or unrefined velocity quadrature')
    history = node['refinement_history']
    require([row['order'] for row in history] == orders[1:orders.index(node['order'])+1],
            'incomplete local velocity-refinement history')
    for refinement in history:
        require(set(refinement) == {'order', 'rate', 'excluded', 'A', 'D'},
                'unexpected local refinement fields')
        require(all(refinement[key] >= 0 for key in ('rate','excluded','A','D')),
                'negative local refinement error')
    last = history[-1]
    passed = (max(last[key] for key in ('rate','excluded','A','D'))
              < a['local_moment_relative_refinement']
              and node['number_density_relative_error']
              < a['local_number_density_relative_error']
              and node['minimum_D_eigenvalue'] >= 0)
    require(node['accepted'] is bool(passed), 'local acceptance flag contradicts thresholds')
    return np.r_[values, node['rho_BH_over_central']], parent, partition


def check_trajectory(record, protocol):
    row, limits = record['summary'], protocol['acceptance']
    horizon = protocol['horizon_over_t_dyn']
    tdyn, core, speed = row['t_dyn_Myr'], row['density_radius_pc'], row['launch_speed_kms']
    require(min(tdyn, core, speed, row['stop_radius_pc']) > 0, 'positive trajectory scales required')
    close(row['horizon_Myr'], horizon*tdyn, 'physical horizon conversion')
    end = row['end_time_Myr']/tdyn
    require(0 < end <= horizon*(1+1e-13), 'invalid terminal time')
    require(row['stop_cause'] in ('horizon','validity_boundary'), 'unregistered stopping cause')
    if row['stop_cause'] == 'horizon':
        close(end, horizon, 'horizon-stopped endpoint')
    require(0 <= row['residence_fraction'] <= end/horizon+1e-13,
            'residence exceeds modeled duration')
    require(set(record['orbit_checks']) == set(ORBIT_CHECK_KEYS), 'missing orbital check')
    checks = record['orbit_checks']
    for name, value in checks.items():
        require(0 <= value < limits[name], f'unaccepted orbital check: {name}')

    track = np.asarray(record['track'], dtype=float)
    require(record['track_columns'] == ['t_over_tdyn','x_over_rho','y_over_rho','z_over_rho',
                                       'vx_over_vdyn','vy_over_vdyn','vz_over_vdyn'],
            'trajectory units/columns changed')
    require(track.shape == (protocol['track_samples'],7), 'unexpected trajectory track size')
    close(track[:,0], np.linspace(0,end,len(track)), 'trajectory sampling times')
    close(track[0,1:4], [protocol['launch']['radius_over_density_radius'],0,0],
          'registered launch position')
    vdyn = core/tdyn*TIME_CONVERSION
    axis = [1,0,0] if row['orientation']=='radial_outward' else [0,1,0]
    close(track[0,4:], speed/vdyn*np.asarray(axis), 'registered launch velocity')
    angular = np.cross(track[:,1:4], track[:,4:])
    close(np.max(np.linalg.norm(angular-angular[0],axis=1)),
          checks['angular_momentum_drift_over_rho_vdyn'], 'recorded angular-momentum drift')
    radii = np.linalg.norm(track[:,1:4],axis=1)
    require(np.max(radii) <= row['stop_radius_pc']/core+1e-8,
            'trajectory continues beyond its validity boundary')
    if row['stop_cause']=='validity_boundary':
        close(radii[-1], row['stop_radius_pc']/core, 'validity-stop endpoint', atol=1e-8)

    crossings = record['crossings']
    require(row['crossing_count']==len(crossings), 'crossing count mismatch')
    require([event['time'] for event in crossings] == sorted(set(event['time'] for event in crossings)),
            'crossing times must be unique and ordered')
    inside, previous, residence = True, 0., 0.
    for event in crossings:
        time = event['time']
        require(previous < time < end, 'crossing lies outside trajectory')
        require(event['direction']==('exit' if inside else 'entry'), 'inconsistent crossing sequence')
        if inside:
            residence += time-previous
        inside, previous = not inside, time
    if inside:
        residence += end-previous
    close(residence/horizon, row['residence_fraction'], 'crossing-based residence')
    turns = record['turning_times']
    require(turns==sorted(set(turns)) and all(0 < time < end for time in turns),
            'turning times must be unique, ordered and interior')
    references = record['independent_radial_references']
    require(len(references)==len(protocol['independent_radial_orders']),
            'missing independent radial refinements')
    require([ref['quadrature_order'] for ref in references]
            == protocol['independent_radial_orders'], 'unregistered independent radial orders')
    require(all(0 <= ref['residence'] <= end+1e-10 for ref in references),
            'invalid independent reference residence')
    for ref in references:
        period, phase, region_time = (ref['radial_period'],ref['initial_phase'],
                                      ref['region_outward_time'])
        require(0 <= ref['pericenter'] < ref['apocenter'], 'invalid radial turning interval')
        require(period > 0 and 0 <= phase <= period and 0 <= region_time <= period/2,
                'invalid radial reference phase')
        close(ref['half_period']*2,period,'radial half-period relation')
        close(ref['angular_momentum'],np.linalg.norm(angular[0]),'reference angular momentum')
        if row['orientation']=='radial_outward':
            close(ref['pericenter'],0.,'radial pericenter')
            close(ref['apocenter'],protocol['launch']['radial_apocenter_over_density_radius'],
                  'radial apocenter from launch energy', atol=1e-10)
        def occupation(phase_time):
            cycles = np.floor(phase_time/period)
            remainder = phase_time-cycles*period
            return (cycles*2*region_time+min(remainder,region_time)
                    +max(remainder-(period-region_time),0.))
        close(ref['residence'],occupation(phase+end)-occupation(phase),
              'independent finite-horizon phase accounting',atol=2e-11)
        close(ref['residence_fraction'],ref['residence']/end,'independent residence normalization')
    error = max(abs(residence-ref['residence'])/horizon for ref in references)
    close(error, checks['independent_residence_difference_over_horizon'],
          'independent residence discrepancy')

    boundaries = sorted(set([0.,end,*turns,*(event['time'] for event in crossings)]))
    levels = record['local_checks']
    expected_orders = protocol['time_quadrature_orders_per_phase']
    require(len(levels) >= 2 and [level['time_order'] for level in levels]
            == expected_orders[:len(levels)], 'incomplete time-quadrature refinement')
    require(len(record['time_refinement'])==len(levels)-1, 'time-refinement record length')
    previous = None
    largest_partition = 0.
    node_count = 0
    for index, level in enumerate(levels):
        nodes, weights = leggauss(level['time_order'])
        expected_times, expected_weights = [], []
        for left,right in zip(boundaries[:-1],boundaries[1:]):
            if right-left < 1e-12:
                continue
            expected_times.extend(left+(nodes+1)*(right-left)/2)
            expected_weights.extend(weights*(right-left)/2*tdyn)
        records = level['nodes']
        require(len(records)==len(expected_times), 'time-quadrature node count')
        actual_times = [node['time_over_tdyn'] for node in records]
        actual_weights = np.asarray([node['weight_Myr'] for node in records])
        require(np.all(actual_weights>0), 'time-quadrature weights must be positive')
        close(actual_times, expected_times, 'registered time-quadrature nodes')
        close(actual_weights, expected_weights, 'registered time-quadrature weights')
        close(actual_weights.sum(), end*tdyn, 'weights must cover entire modeled duration')
        integrands, parents = [], []
        for node in records:
            values, parent, partition = check_local_node(node,protocol)
            integrands.append(values); parents.append(parent)
            largest_partition = max(largest_partition,partition)
        values = np.asarray(integrands)
        total = actual_weights @ values
        parent_integral = actual_weights @ np.asarray(parents)
        integrated_partition = abs((total[0]+total[1])/parent_integral-1)
        largest_partition = max(largest_partition,integrated_partition)
        require(integrated_partition < 1e-12, 'integrated incoming flux partition does not close')
        normalized = total/np.array([1,1,speed,speed**2,horizon*tdyn])
        node_count += len(records)
        if previous is not None:
            refinement = record['time_refinement'][index-1]
            require(refinement['time_order']==level['time_order'], 'time-refinement order mismatch')
            delta = abs(normalized-previous)
            relative = delta[:4]/np.maximum(abs(normalized[:4]),np.finfo(float).tiny)
            close(refinement['budget_relative_refinement'], relative,
                  'weighted budget time refinement', atol=2e-12)
            close(refinement['density_exposure_absolute_refinement'], delta[4],
                  'density exposure time refinement', atol=2e-12)
            passed = (max(relative) < limits['budget_time_relative_refinement']
                      and delta[4] < limits['density_exposure_time_absolute_refinement']
                      and all(node['accepted'] for node in records))
            require(refinement['accepted'] is bool(passed), 'time acceptance flag contradicts thresholds')
            if index < len(levels)-1:
                require(not passed, 'calculation continued after declared convergence stop')
        previous = normalized
    require(record['time_refinement'][-1]['accepted'] is True, 'unaccepted final time refinement')
    close([row[key] for key in BUDGET_KEYS], normalized, 'reported path budget integrals')
    require(row['accepted'] is True, 'trajectory not accepted')
    return node_count, largest_partition


def verify(record=False):
    protocol = json.loads(PROTOCOL.read_text())
    summary = json.loads((OUT/'summary.json').read_text())
    require(len(protocol['clusters'])==len(set(protocol['clusters']))==6, 'six distinct hosts required')
    require(protocol['launch']['orientations']==['radial_outward','tangential'], 'registered launch directions')
    encounter = ROOT/protocol['encounter_protocol']
    ep = json.loads(encounter.read_text())
    require(ep['remnant_mass_Msun']==protocol['remnant_mass_Msun'], 'encounter mass mismatch')
    require(protocol['bmax_over_local_scale'] in ep['bmax_over_local_scale'], 'unregistered impact cutoff')
    environment = ROOT/'results/cluster-comparison-environment/comparison.csv'
    with environment.open() as handle:
        environments = {item['cluster']:item for item in csv.DictReader(handle)}
    input_files = [PROTOCOL, encounter, environment,
                   *(ROOT/'analysis/residence').glob('*.py'),
                   *(ROOT/name for name in IMPLEMENTATION)]
    output_files = [OUT/'comparison.csv',OUT/'summary.json']
    implementation_hashes = {name:digest(ROOT/name) for name in IMPLEMENTATION}
    all_rows, total_nodes, largest_partition = [], 0, 0.
    for cluster in protocol['clusters']:
        host = ROOT/f'results/cluster-comparison/{cluster}/reference/host.npz'
        source = OUT/f'{cluster}.json'
        result = json.loads(source.read_text())
        finite_tree(result)
        input_files.extend([host,host.with_name('summary.json')]); output_files.append(source)
        require(result['cluster']==cluster, 'cluster identifier mismatch')
        require(result['host_sha256']==digest(host), 'stale host input')
        require(result['protocol_sha256']==digest(PROTOCOL), 'stale trajectory protocol')
        require(result['encounter_protocol_sha256']==digest(encounter), 'stale encounter protocol')
        require(result['environment_sha256']==digest(environment), 'stale environmental input')
        require(result['implementation_sha256']==implementation_hashes, 'stale implementation')
        require(json.loads(host.with_name('summary.json').read_text())['all_numerical_checks_passed'] is True,
                'unaccepted host reconstruction')
        require(environments[cluster]['all_numerical_checks_passed']=='True', 'unaccepted environment calculation')
        with np.load(host,allow_pickle=False) as data:
            expected_stop = min(float(data['rt']),float(environments[cluster]['refined_minimum_1pct_radius_pc']))
        require(len(result['records'])==2, 'two trajectories required per host')
        require([item['summary']['orientation'] for item in result['records']]
                == protocol['launch']['orientations'], 'missing or duplicate launch orientation')
        for item in result['records']:
            row = item['summary']
            require(row['cluster']==cluster, 'row cluster mismatch')
            close(row['stop_radius_pc'],expected_stop,'declared validity boundary')
            count, partition = check_trajectory(item,protocol)
            total_nodes += count; largest_partition = max(largest_partition,partition)
            all_rows.append(row)
        require(result['all_numerical_checks_passed'] is True, 'unaccepted host trajectory calculation')
    require(summary['revision']==protocol['revision'] and summary['scope']==protocol['scope'],
            'aggregate protocol metadata mismatch')
    require(summary['interpretation']==protocol['interpretation'], 'aggregate interpretation changed')
    require(summary['trajectories']==summary['accepted_trajectories']==len(all_rows)==12,
            'twelve accepted trajectories required')
    require(summary['all_numerical_checks_passed'] is True, 'aggregate not accepted')
    require(summary['collisionless_trajectories_performed'] is True
            and summary['collisional_trajectories_performed'] is False
            and summary['residence_response_evaluated'] is False,
            'collisionless controls must remain separate from collisional results')
    with (OUT/'comparison.csv').open() as handle:
        csv_rows = list(csv.DictReader(handle))
    require(len(csv_rows)==len(all_rows), 'comparison CSV row count')
    for csv_row,row in zip(csv_rows,all_rows,strict=True):
        require(csv_row=={key:str(value) for key,value in row.items()}, 'CSV/JSON row mismatch')
    provenance = {
        'revision':protocol['revision'], 'scope':protocol['interpretation'],
        'checked_trajectories':len(all_rows), 'checked_time_nodes':total_nodes,
        'maximum_local_or_integrated_rate_partition_error':largest_partition,
        'input_sha256':{str(path.relative_to(ROOT)):digest(path) for path in sorted(set(input_files))},
        'output_sha256':{str(path.relative_to(ROOT)):digest(path) for path in sorted(set(output_files))},
    }
    path = OUT/'provenance.json'
    if record:
        path.write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
    else:
        require(json.loads(path.read_text())==provenance, 'stale collisionless provenance')
    print(f'Verified 12 collisionless trajectories and {total_nodes} recorded time nodes; '
          'no collisional residence or response result.')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record',action='store_true')
    verify(parser.parse_args().record)
