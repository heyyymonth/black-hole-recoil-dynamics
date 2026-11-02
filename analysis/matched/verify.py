"""Read-only integrity/accounting checks; NumPy only, no physical integration."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analysis.matched.analyze import analyze,settings,name
OUT=ROOT/'results/matched-residence'
PROTOCOL=ROOT/'analysis/matched/protocol.json'


def require(condition,message):
    if not condition:raise AssertionError(message)
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def compare_tree(actual,expected,path='record'):
    if isinstance(expected,dict):
        require(set(actual)==set(expected),f'{path}: field mismatch')
        for k in expected:compare_tree(actual[k],expected[k],path+'.'+k)
    elif isinstance(expected,list):
        require(len(actual)==len(expected),f'{path}: list length mismatch')
        for i,(a,b) in enumerate(zip(actual,expected)):compare_tree(a,b,f'{path}[{i}]')
    elif isinstance(expected,float):
        require(np.isfinite(actual) and np.isfinite(expected) and np.isclose(actual,expected,rtol=2e-11,atol=2e-13),f'{path}: numeric mismatch')
    else:require(actual==expected,f'{path}: mismatch')
def check_hashes(records):
    for file,sha in records.items():
        require(not Path(file).is_absolute() and '..' not in Path(file).parts,'unsafe path')
        require(digest(ROOT/file)==sha,'stale evidence '+file)


def check_quadrature(row,p):
    orders=p['kernel']['velocity_quadrature_orders'];limit=p['kernel']['local_acceptance']
    require(row['order'] in orders[1:],'unrefined quadrature')
    require([h['order'] for h in row['refinement_history']]==orders[1:orders.index(row['order'])+1], 'incomplete refinement')
    last=row['refinement_history'][-1]
    require(all(np.isfinite(last[k]) and 0<=last[k]<limit['moment_relative_refinement'] for k in ('rate','excluded','A','D')),'unaccepted moment quadrature')
    require(0<=row['number_density_relative_error']<limit['number_density_relative_error'],'unaccepted density quadrature')
    require(0<=row['rate_partition_relative_error']<1e-12,'incoming flux partition')
    require(row['accepted'] and row['minimum_D_eigenvalue']>=0,'unaccepted vertex')


def verify_sampling_records(p):
    """Check retained diagnostics without importing the physical sampler.

    Draws are not retained, so this verifies the recorded Monte Carlo summaries,
    their acceptance arithmetic and their provenance; it does not remeasure
    their standard errors or replay the random draws.
    """
    design_path=OUT/'sampling-validation-design.json'
    design=json.loads(design_path.read_text())
    result=json.loads((OUT/'sampling-validation.json').read_text())
    require(result['design_sha256']==digest(design_path),'stale sampler design')
    require(design['revision']==result['revision']=='matched-sampling-validation-1','sampler validation revision')
    require(design['seed']==2026110101 and design['samples_per_state_per_grid']==40000,'sampler validation controls')
    require(design['grids']==p['transport']['grids'],'sampler validation grids')
    require(len(design['states'])==6,'sampler validation state inventory')
    labels=[row['label'] for row in design['states']]
    require(len(set(labels))==len(labels),'duplicate sampling states')
    for row in design['states']:
        state=np.asarray(row['state'],dtype=float)
        require(state.shape==(6,) and np.isfinite(state).all(),'invalid validation state')
        compare_tree(float(np.linalg.norm(state[:3])),row['radius'])
        compare_tree(float(np.linalg.norm(state[3:])),row['speed'])
    require(design['source_sha256']==result['source_sha256'],'sampler source inventory mismatch')
    check_hashes(design['source_sha256'])
    require(not result['sources_changed_during_run'],'sampler sources changed during run')
    expected_inventory=[(grid,label) for grid in design['grids'] for label in labels]
    require([(row['grid'],row['label']) for row in result['records']]==expected_inventory,'sampler result inventory')
    for row in result['records']:
        require(np.isfinite(row['rate_per_t_dyn']) and row['rate_per_t_dyn']>0,'invalid sampling rate')
        for key,shape in [('mean',(3,)),('raw_second',(3,3))]:
            check=row[key]
            arrays={field:np.asarray(check[field],dtype=float) for field in
                    ('expected','estimate','standard_error','residual','standardized_residual')}
            require(all(value.shape==shape and np.isfinite(value).all() for value in arrays.values()),'invalid sampler moment shape/value')
            se=arrays['standard_error'];residual=arrays['estimate']-arrays['expected']
            require(np.all(se>=0),'negative sampler standard error')
            np.testing.assert_allclose(arrays['residual'],residual,rtol=2e-11,atol=2e-13)
            standardized=np.divide(residual,se,out=np.zeros_like(residual),where=se>0)
            np.testing.assert_allclose(arrays['standardized_residual'],standardized,rtol=2e-11,atol=2e-13)
            passed=bool(np.all(abs(residual)<6*se+1e-13))
            require(check['passed']==passed and passed,'failed sampling moment check')
        require(row['passed']==all(row[key]['passed'] for key in ('mean','raw_second')),'sampling aggregate status')
    require(result['all_checks_passed'] and all(row['passed'] for row in result['records']),'sampling validation failed')


def verify():
    p=json.loads(PROTOCOL.read_text());protocol_hash=digest(PROTOCOL)
    for grid in p['transport']['grids']:
        path=OUT/f'kernel-{grid}.npz';k=json.loads(path.with_suffix('.json').read_text());a=dict(np.load(path,allow_pickle=False))
        require(k['protocol_sha256']==protocol_hash and k['host_sha256']==digest(ROOT/p['host_path']),'stale kernel inputs')
        require(k['encounter_protocol_sha256']==digest(ROOT/p['encounter_protocol']),'stale encounter cutoffs')
        check_hashes(k['implementation_sha256']);require(k['table_sha256']==digest(path),'table changed')
        axes=[np.array(p['kernel']['coarse_axes'][key]) for key in ('radius','speed','pitch')]
        if grid=='fine':axes=[np.sort(np.r_[x,(x[1:]+x[:-1])/2]) for x in axes]
        for key,axis in zip(('radius','speed','pitch'),axes):np.testing.assert_array_equal(a[key],axis)
        shape=tuple(map(len,axes));count=int(np.prod(shape))
        require(k['shape']==list(shape) and k['vertex_count']==count and len(k['vertices'])==count,'grid inventory')
        require(a['completed'].all() and a['accepted'].all() and k['completed'] and k['all_vertex_quadrature_passed'] and k['all_holdout_quadrature_passed'],'incomplete kernel')
        for key in ['rate','excluded','A','D']:
            require(np.isfinite(a[key]).all(),'nonfinite moment')
        require(np.all(a['rate']>0) and np.all(a['excluded']>=0),'negative rate')
        require(np.min(np.linalg.eigvalsh(a['D']))>=-1e-13,'nonpositive diffusion')
        require(np.allclose(a['D'],a['D'].swapaxes(-1,-2),rtol=0,atol=1e-13),'asymmetric diffusion')
        for i,row in enumerate(k['vertices']):
            check_quadrature(row,p);idx=np.unravel_index(i,shape)
            compare_tree(row['coordinates'],[float(axes[j][v]) for j,v in enumerate(idx)])
            require(a['order'][idx]==row['order'],'quadrature order mismatch')
        require(len(k['holdout_points'])==p['kernel']['holdout_points']+8,'holdout inventory')
        for row in k['holdout_points']:check_quadrature(row['direct_quadrature'],p)
        compare_tree(k['holdout_error_maxima'],{key:max(row['errors'][key] for row in k['holdout_points']) for key in k['holdout_error_maxima']})
    verify_sampling_records(p)
    n=p['transport']['paths_per_setting'];H=p['horizon_over_t_dyn']
    for index,setting in enumerate(settings(p)):
        grid,dt,orientation,method=setting;path=OUT/f'{name(*setting)}.npz'
        a=dict(np.load(path,allow_pickle=False));r=json.loads(path.with_suffix('.json').read_text())
        require(r['protocol_sha256']==protocol_hash and r['samples_sha256']==digest(path),'stale trajectory inputs/samples')
        require(r['host_sha256']==digest(ROOT/p['host_path']) and r['encounter_protocol_sha256']==digest(ROOT/p['encounter_protocol']),'stale physical inputs')
        require(r['kernel_sha256']==digest(OUT/f'kernel-{grid}.npz') and r['kernel_record_sha256']==digest(OUT/f'kernel-{grid}.json'),'stale kernel provenance')
        check_hashes(r['implementation_sha256'])
        require([r[k] for k in ('grid','dt','orientation','method')]==list(setting),'unregistered setting')
        require(r['seed']==[p['transport']['seed'],index] and r['paths']==n and r['horizon']==H,'sample controls')
        radius_edge=p['kernel']['coarse_axes']['radius'][-1]
        velocity_edge=p['kernel']['coarse_axes']['speed'][-1]
        require(r['spatial_stop']==min(radius_edge,r['inherited_force_domain_over_density_radius']),'unregistered spatial stop')
        require(r['velocity_stop']==velocity_edge,'unregistered velocity stop')
        require(a['states'].shape==(n,6),'path count/shape')
        for key,value in a.items():
            require(value.shape[0]==n,'unaligned path rows')
            if value.dtype.kind in 'fiu':require(np.isfinite(value).all(),'nonfinite path statistic')
        require(np.all((a['end_time']>=0)&(a['end_time']<=H+1e-12)),'invalid duration')
        require(np.all((a['residence']>=-1e-12)&(a['residence']<=a['end_time']+1e-10)),'occupation exceeds modeled time')
        require(np.all(a['exposure']>=0) and np.all(a['exposure']<=a['end_time']+1e-7),'density exposure outside its bounds')
        require(np.allclose(a['normalized_residence'],a['residence']/H,rtol=0,atol=1e-13),'residence normalization')
        require(np.allclose(a['normalized_exposure'],a['exposure']/H,rtol=0,atol=1e-13),'exposure normalization')
        causes=a['stop_cause'];allowed=['horizon','spatial_domain','velocity_domain','excluded_channel']
        require(np.isin(causes,allowed).all(),'unregistered stop')
        require(np.array_equal(a['alive'],causes=='horizon'),'survival accounting')
        require(np.all(a['end_time'][a['alive']]==H),'horizon timing')
        require(np.all(a['elastic_events']>=0) and np.all(a['elastic_events']==np.floor(a['elastic_events'])),'invalid event counts')
        require(np.all(a['evaluated_steps']>=0) and np.all(a['evaluated_steps']==np.floor(a['evaluated_steps'])),'invalid evaluated step counts')
        require(np.all(a['evaluated_steps']<=int(np.ceil(H/dt))),'too many evaluated steps')
        if method=='diffusion':require(not a['elastic_events'].any(),'Gaussian arm contains counted events')
        spatial=causes=='spatial_domain'
        require(np.allclose(np.linalg.norm(a['states'][spatial,:3],axis=1),r['spatial_stop'],rtol=0,atol=1e-10),'spatial endpoint')
        velocity=causes=='velocity_domain'
        require(np.all(np.linalg.norm(a['states'][velocity,3:],axis=1)>=r['velocity_stop']-1e-12),'velocity-stop endpoint')
        modeled=~(spatial|velocity)
        require(np.all(np.linalg.norm(a['states'][modeled,:3],axis=1)<r['spatial_stop']+1e-12),'modeled path outside spatial domain')
        require(np.all(np.linalg.norm(a['states'][modeled,3:],axis=1)<r['velocity_stop']+1e-12),'modeled path outside velocity domain')
    result=analyze();stored=json.loads((OUT/'summary.json').read_text());compare_tree(stored,result)
    print(f"Verified {result['ensembles']} matched ensembles, {result['total_trajectories']} retained paths and both kernel grids; no trajectories rerun.")
    return result


if __name__=='__main__':
    argparse.ArgumentParser(description=__doc__).parse_args();verify()
