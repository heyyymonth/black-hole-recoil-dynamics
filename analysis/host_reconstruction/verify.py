"""Read-only calculation integrity checks; --record pins a reviewed result set."""
from pathlib import Path
import argparse,csv,hashlib,json,math
ROOT=Path(__file__).resolve().parents[2]
ANALYSIS=['analysis/host_reconstruction/'+n for n in ['run.py','host.py','verify.py','test_host.py','protocol.json','requirements-lock.txt','environment.json']]+['analysis/transport/'+n for n in ['coefficients.py','run.py','test_coefficients.py','protocol.json']]
OUTPUTS=['results/host-reconstruction/'+n for n in ['summary.json','host.npz','species.csv','df-checks.csv','projection-checks.csv','published-profile-comparison.csv']]+['results/local-transport-kernels/'+n for n in ['summary.json','states.csv','species.csv','tensors.json']]

def load(path):
    def reject(s):raise ValueError('Nonfinite JSON: '+s)
    return json.loads((ROOT/path).read_text(),parse_constant=reject)

def rows(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))

def hashes(paths):return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}

def validate():
    hp=load('analysis/host_reconstruction/protocol.json');kp=load('analysis/transport/protocol.json')
    h=load('results/host-reconstruction/summary.json');k=load('results/local-transport-kernels/summary.json')
    assert h['protocol_sha256']==hashes(['analysis/host_reconstruction/protocol.json'])['analysis/host_reconstruction/protocol.json']
    assert k['protocol_sha256']==hashes(['analysis/transport/protocol.json'])['analysis/transport/protocol.json']
    assert k['host_sha256']==hashes(['results/host-reconstruction/host.npz'])['results/host-reconstruction/host.npz']
    assert h['all_numerical_checks_passed'] and k['all_numerical_checks_passed']
    assert all(h['checks_passed'].values()) and all(k['checks'].values())
    for key,value in h['checks'].items():assert 0<=value<(.002 if key.startswith('projection') else .001),(key,value)
    assert k['max_A_refinement_error']<kp['acceptance']['relative_vector_and_tensor_refinement']
    assert k['max_D_refinement_relative_error']<kp['acceptance']['relative_vector_and_tensor_refinement']
    states=rows('results/local-transport-kernels/states.csv');species=rows('results/local-transport-kernels/species.csv');host_species=rows('results/host-reconstruction/species.csv');tensors=load('results/local-transport-kernels/tensors.json')
    assert len(states)==k['states']==20 and len(tensors)==len(states)
    assert math.isclose(sum(float(r['M_Msun']) for r in host_species),h['summary']['M_Msun'],rel_tol=1e-12)
    for state,tensor in zip(states,tensors):
        assert int(state['state'])==tensor['state']
        terms=[r for r in species if r['state']==state['state']]
        total=float(state['trace_D_over_lnLambda_kms2_Myr'])
        assert math.isclose(sum(float(r['trace_D_over_lnLambda_kms2_Myr']) for r in terms),total,rel_tol=1e-12)
        assert math.isclose(sum(float(r['trace_D_fraction']) for r in terms),1,rel_tol=1e-12)
        assert math.isclose(sum(tensor['D_over_lnLambda'][i][i] for i in range(3)),total,rel_tol=1e-12)
        assert float(state['minimum_D_eigenvalue'])>=0
        assert float(state['trace_identity_relative_error'])<kp['acceptance']['trace_identity_relative']
    return {'host':'host-reconstruction-1','coefficients':'local-transport-kernels-1','scope':'one joint posterior host and local weak-scattering kernels; no remnant trajectories','source_files':hp['sources'],'input_sha256':hashes(ANALYSIS),'output_sha256':hashes(OUTPUTS)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    record=validate();path=ROOT/'results/calculations.json'
    if args.record:path.write_text(json.dumps(record,indent=2)+'\n')
    else:assert load('results/calculations.json')==record,'Calculation provenance changed'
    print('Verified fixed host and 20 local collision states; no trajectory execution.')
if __name__=='__main__':main()
