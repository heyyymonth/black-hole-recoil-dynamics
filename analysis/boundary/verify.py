"""Read-only boundary evidence checks; --record binds a reviewed calculation."""
from pathlib import Path
import argparse,csv,hashlib,json,math
ROOT=Path(__file__).resolve().parents[2]
INPUTS=['analysis/boundary/'+n for n in ['protocol.json','sources.json','model.py','run.py','test_boundary.py','verify.py']]+['results/host-reconstruction/host.npz','results/calculations.json','data/observed-clusters/derived/catalog_scalars.csv']
OUTPUTS=['results/outer-boundary/'+n for n in ['summary.json','scenarios.csv','radial-profile.csv']]

def hashes(paths):return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
def read(path):return json.loads((ROOT/path).read_text(),parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
def rows(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))

def validate():
    p=read('analysis/boundary/protocol.json');s=read(OUTPUTS[0]);data=rows(OUTPUTS[1]);profile=rows(OUTPUTS[2])
    assert s['revision']==p['revision']=='outer-boundary-1'
    assert s['protocol_sha256']==hashes(['analysis/boundary/protocol.json'])['analysis/boundary/protocol.json']
    assert not s['transport_simulations_performed'] and not s['galactic_orbits_integrated']
    assert s['all_numerical_checks_passed'] and len(data)==len(p['scenarios'])==4
    assert s['checks']['max_root_relative_residual']<p['acceptance']['root_relative_residual']
    assert s['checks']['max_density_integral_relative_error']<p['acceptance']['density_integral_relative']
    assert s['checks']['max_mass_interpolation_refinement_relative']<p['acceptance']['mass_interpolation_refinement_relative']
    sources=read(p['sources']);catalog=next(r for r in rows(p['catalog']) if r['cluster']=='NGC3201')
    assert float(catalog['tidal_catalog_pc'])==sources['catalog_snapshot']['tidal_catalog_pc']
    raw=ROOT/sources['catalog_snapshot']['path']
    if raw.exists():assert hashlib.sha256(raw.read_bytes()).hexdigest()==sources['catalog_snapshot']['sha256']
    for scenario,row in zip(p['scenarios'],data):
        assert row['scenario']==scenario['name']
        for key in ['R_pc','vc_kms']:assert float(row[key])==scenario[key]
        numbers={k:float(v) for k,v in row.items() if k!='scenario'}
        assert all(math.isfinite(v) for v in numbers.values())
        assert 0<numbers['mass_outside_rJ_fraction']<1
        assert numbers['r_valid_0.01_pc']<numbers['r_valid_0.03_pc']<numbers['r_valid_0.1_pc']<numbers['rJ_extended_pc']
        assert numbers['exact_inward_saddle_pc']<numbers['rJ_extended_pc']<numbers['exact_outward_saddle_pc']
        if row['scenario']=='baseline':
            assert all(math.isclose(v,s['baseline'][k],rel_tol=1e-14,abs_tol=1e-20) for k,v in numbers.items())
    assert len(profile)==320
    eta=[float(r['max_external_to_internal_force']) for r in profile]
    assert all(0<a<b for a,b in zip(eta,eta[1:]))
    return {'revision':p['revision'],'scope':p['scope'],'input_sha256':hashes(INPUTS),'output_sha256':hashes(OUTPUTS)}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    record=validate();path=ROOT/'results/outer-boundary/provenance.json'
    if args.record:path.write_text(json.dumps(record,indent=2)+'\n')
    else:assert read(str(path.relative_to(ROOT)))==record,'Boundary evidence changed'
    print('Verified four static boundary scenarios; no orbit or remnant trajectory execution.')
if __name__=='__main__':main()
