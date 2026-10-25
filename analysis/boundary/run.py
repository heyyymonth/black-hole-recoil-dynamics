"""Execute the registered static-boundary calculation, without remnant paths."""
from pathlib import Path
import csv, hashlib, json, sys
import numpy as np
import scipy
from scipy.integrate import quad
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'analysis/host_reconstruction'))
from host import Host
from model import EnclosedHost, diagnostics, max_force_ratio


def main():
    protocol = json.loads((ROOT/'analysis/boundary/protocol.json').read_text())
    host = Host(dict(np.load(ROOT/protocol['host'], allow_pickle=False)))
    model = EnclosedHost(host)
    records = []
    density_errors = []
    refinement_errors = []
    for scenario in protocol['scenarios']:
        args = scenario['R_pc'], scenario['vc_kms'], protocol['neglected_external_force_tolerances']
        result = diagnostics(model, *args)
        coarse = diagnostics(EnclosedHost(host, stride=2), *args)
        for key in ['rJ_extended_pc']+[f'r_valid_{e:g}_pc' for e in args[-1]]:
            r = result[key]
            integrated = quad(lambda x: 4*np.pi*x*x*host.rho_spline(x)[host.physical].sum(), 0, r, epsabs=1e-5, epsrel=1e-8, limit=500)[0]
            density_errors.append(abs(integrated/model.enclosed(r)-1))
            refinement_errors.append(abs(coarse[key]/r-1))
        records.append({'scenario': scenario['name'], **result})
    checks = {
        'max_root_relative_residual': max(v for row in records for k,v in row.items() if 'residual' in k),
        'max_density_integral_relative_error': max(density_errors),
        'max_mass_interpolation_refinement_relative': max(refinement_errors),
    }
    a = protocol['acceptance']
    passed = checks['max_root_relative_residual']<a['root_relative_residual'] and checks['max_density_integral_relative_error']<a['density_integral_relative'] and checks['max_mass_interpolation_refinement_relative']<a['mass_interpolation_refinement_relative']
    if not passed:
        raise AssertionError(checks)
    out = ROOT/'results/outer-boundary';out.mkdir(exist_ok=True)
    with (out/'scenarios.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=records[0], lineterminator='\n');writer.writeheader();writer.writerows(records)
    R, vc = records[0]['R_pc'], records[0]['vc_kms']
    with (out/'radial-profile.csv').open('w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n');writer.writerow(['r_pc','enclosed_mass_fraction','max_external_to_internal_force','linear_tide_ratio'])
        for r in np.geomspace(.01*host.rh, min(host.rt, .2*R), 320):
            mass = float(model.enclosed(r))
            writer.writerow([r,mass/host.M,float(max_force_ratio(r,model.enclosed,host.G,R,vc)),vc**2/R**2*r**3/(host.G*mass)])
    summary = {'revision':protocol['revision'],'transport_simulations_performed':False,'galactic_orbits_integrated':False,
               'host_mass_msun':host.M,'host_rh_pc':host.rh,'host_formal_rt_pc':host.rt,'G_pc_kms2_msun':host.G,
               'baseline':records[0], 'checks':checks,'all_numerical_checks_passed':passed,
               'numerical_environment':{'python':sys.version.split()[0],'numpy':np.__version__,'scipy':scipy.__version__},
               'protocol_sha256':hashlib.sha256((ROOT/'analysis/boundary/protocol.json').read_bytes()).hexdigest(),
               'interpretation':protocol['interpretation']}
    (out/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))
if __name__=='__main__':main()
