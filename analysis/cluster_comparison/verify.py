"""Verify recorded cross-cluster results and source lineage without rerunning fits."""
from pathlib import Path
import argparse,csv,hashlib,json,math
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/cluster-comparison'
CLUSTERS=['NGC0104','NGC1851','NGC3201','NGC5139','NGC6121','NGC6397']

def load(p):
    return json.loads(p.read_text(),parse_constant=lambda x:(_ for _ in ()).throw(ValueError(x)))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def hashes(paths):return {str(p.relative_to(ROOT)):sha(p) for p in sorted(set(paths))}

def validate():
    summary=load(OUT/'summary.json');runs=summary['runs']
    assert summary['planned']==54 and summary['attempted']==54 and len(runs)==54
    assert not summary['remnant_trajectories_performed'] and not summary['pilot_is_converged_posterior']
    inventory=load(OUT/'posterior-inventory.json')['clusters'];assert [r['cluster'] for r in inventory]==CLUSTERS
    for c in CLUSTERS:
        selected=load(OUT/c/'selection.json')['rows'];actual=[r for r in runs if r['cluster']==c]
        assert len(selected)==len(actual)==9
        assert sum(r['kind']=='reference' for r in actual)==1
        for chosen,r in zip(selected,actual):
            for key in ['cluster','kind','label','row_index','theta','normalized_source_weight']:assert chosen[key]==r[key]
            assert r==load(OUT/c/r['label']/'result.json')
            assert r['host_status'] in ['passed','failed']
            if r['host_status']=='failed':
                assert r.get('failure_log');continue
            h=r['host_summary'];assert h['all_numerical_checks_passed'] and all(h['checks_passed'].values())
            assert h['theta']==r['theta'] and h['joint_row_index']==r['row_index']
            assert h['protocol_sha256']==sha(ROOT/r['protocol_path'])
            assert r['core_profile_status']=='passed'
            for p in r['core_profile']:
                assert p['status']=='passed' and p['refinement_relative_error']<.005
                assert p['trace_D_per_lnLambda']>0 and 0<=p['BH_trace_D_fraction']<=1+1e-12
                assert 0<=p['BH_density_fraction']<=1+1e-12
            m=r['BH_metrics']
            if m['BH_absent']:
                assert m['BH_mass_Msun']==m['BH_expected_number']==0
                assert m['probe_to_BH_reservoir_mass_ratio'] is None
                assert all(p['BH_trace_D_fraction']==0 for p in r['core_profile'])
            else:assert math.isclose(m['probe_to_BH_reservoir_mass_ratio']*m['BH_mass_Msun'],60,rel_tol=1e-12)
            if r['kind']=='reference':
                assert r['transport_status']=='passed'
                k=r['transport_summary'];assert k['all_numerical_checks_passed'] and all(k['checks'].values()) and k['states']==20
                assert k['host_sha256']==sha(OUT/c/'reference/host.npz')==r['host_arrays_sha256']
                with (OUT/c/'reference/transport/species.csv').open() as f:species=list(csv.DictReader(f))
                with (OUT/c/'reference/transport/states.csv').open() as f:states=list(csv.DictReader(f))
                for row in states:
                    terms=[x for x in species if x['state']==row['state']]
                    assert math.isclose(sum(float(x['trace_D_fraction']) for x in terms),1,rel_tol=1e-12)
                assert math.isclose(float(states[0]['trace_D_over_lnLambda_kms2_Myr']),r['core_profile'][0]['trace_D_per_lnLambda'],rel_tol=.005)
    environment=load(ROOT/'results/cluster-comparison-environment/summary.json')
    assert environment['all_six_clusters_present'] and len(environment['clusters'])==6
    for r in environment['clusters']:
        assert r['all_numerical_checks_passed'],r['cluster']
        assert len(r['cases'])==4 and len(r['static_controls'])==3
        assert r['host_sha256']==sha(OUT/r['cluster']/'reference/host.npz')
        assert r['protocol_sha256']==sha(ROOT/'data/cluster-comparison-environment/protocol.json')
    return {'reference_hosts':6,'reference_local_states':120,'pilot_attempts':48,
            'accepted_pilot_hosts':sum(r['kind']=='pilot' and r['host_status']=='passed' for r in runs),
            'galactic_reference_orbits':24,'static_controls':18,'published_derived_realizations':sum(r['derived_rows'] for r in inventory)}

def record():
    report=validate()
    inputs=list((ROOT/'analysis/cluster_comparison').rglob('*.py'))+list((ROOT/'analysis/cluster_comparison').rglob('*.json'))
    inputs += [ROOT/p for p in ['analysis/host_reconstruction/run.py','analysis/host_reconstruction/host.py','analysis/transport/run.py','analysis/transport/coefficients.py','analysis/transport/protocol.json','analysis/host_reconstruction/requirements-lock.txt','data/observed-clusters/sources.json','analysis/orbit/galaxy.py','analysis/orbit/orbits.py','analysis/boundary/model.py']]
    inputs += [p for p in (ROOT/'data/cluster-comparison-environment').glob('*') if p.is_file()]
    inputs += [ROOT/'docs/six-cluster-data-quality.md']
    notebook=ROOT/'analysis/cluster_comparison/cluster_comparison.ipynb'
    nb=load(notebook)
    assert all(c.get('execution_count') is not None and not any(o.get('output_type')=='error' for o in c.get('outputs',[])) for c in nb['cells'] if c['cell_type']=='code')
    outputs=[p for folder in [OUT,ROOT/'results/cluster-comparison-environment'] for p in folder.rglob('*') if p.is_file() and p.name!='provenance.json']
    outputs += [notebook,ROOT/'docs/six-cluster-comparison.md',ROOT/'docs/six-cluster-environment.md']
    return {'revision':'six-cluster-comparison-1','validated_counts':report,'input_sha256':hashes(inputs),'output_sha256':hashes(outputs),
            'scratch_policy':'Pilot meshes are reproducible scratch. Their hashes, selected complete parameter rows and validation outputs are retained; scratch is not required for read-only verification.'}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--record',action='store_true');args=parser.parse_args()
    value=record();p=OUT/'provenance.json'
    if args.record:p.write_text(json.dumps(value,indent=2)+'\n')
    else:assert value==load(p),'Comparison evidence changed'
    print(json.dumps(value['validated_counts']))
if __name__=='__main__':main()
