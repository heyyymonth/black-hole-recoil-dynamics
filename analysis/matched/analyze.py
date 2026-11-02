"""Summarize recorded independent ensembles without running physical dynamics."""
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analysis.matched.statistics import summarize_samples,compare_samples
OUT=ROOT/'results/matched-residence'
PROTOCOL=ROOT/'analysis/matched/protocol.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def name(grid,dt,orientation,method):return f'{grid}-dt{dt:.3f}-{orientation}-{method}'
def settings(p):
    return list(itertools.product(p['transport']['grids'],p['transport']['dt_over_t_dyn'],
                                  p['launch']['orientations'],p['transport']['methods']))


def analyze():
    p=json.loads(PROTOCOL.read_text());confidence=p['reporting']['confidence_level']
    records={};data={};inputs={str(PROTOCOL.relative_to(ROOT)):digest(PROTOCOL)}
    for setting in settings(p):
        key=name(*setting);path=OUT/f'{key}.npz';meta=path.with_suffix('.json')
        sample=dict(np.load(path,allow_pickle=False));record=json.loads(meta.read_text())
        if record['samples_sha256']!=digest(path):raise ValueError('Sample digest mismatch')
        for f in [path,meta]:inputs[str(f.relative_to(ROOT))]=digest(f)
        records[key]={'grid':setting[0],'dt':setting[1],'orientation':setting[2],'method':setting[3],
            'statistics':summarize_samples(sample['normalized_residence'],sample['normalized_exposure'],
                sample['end_time'],sample['stop_cause'],sample['elastic_events'],
                p['horizon_over_t_dyn'],confidence)}
        if setting[3]=='jump':
            zero=sample['elastic_events']==0
            full_horizon=zero & (sample['stop_cause']=='horizon')
            records[key]['statistics']['zero_event_paths']={
                'count':int(zero.sum()),'fraction':float(zero.mean()),
                'horizon_count':int(full_horizon.sum()),
                'horizon_fraction':float(full_horizon.mean()),
                'interpretation':'Observed zero resolved finite-event paths; horizon subset excludes numerical and excluded-channel stops. Descriptive counts, not a confidence statement.'}
        data[setting]=sample
    comparisons=[]
    for grid,dt,orientation in itertools.product(p['transport']['grids'],p['transport']['dt_over_t_dyn'],p['launch']['orientations']):
        j=data[grid,dt,orientation,'jump'];d=data[grid,dt,orientation,'diffusion']
        comparisons.append({'grid':grid,'dt':dt,'orientation':orientation,
            'sign':'jump minus diffusion','residence':compare_samples(j['normalized_residence'],d['normalized_residence'],confidence),
            'exposure':compare_samples(j['normalized_exposure'],d['normalized_exposure'],confidence)})
    refinements=[]
    def refinement(kind,fine,coarse):
        refinements.append({'kind':kind,'first':name(*fine),'second':name(*coarse),
                            'sign':'finer minus coarser',
                            'residence':compare_samples(data[fine]['normalized_residence'],data[coarse]['normalized_residence'],confidence),
                            'exposure':compare_samples(data[fine]['normalized_exposure'],data[coarse]['normalized_exposure'],confidence)})
    coarse_grid,fine_grid=p['transport']['grids'];coarse_dt,fine_dt=p['transport']['dt_over_t_dyn']
    for grid,orientation,method in itertools.product(p['transport']['grids'],p['launch']['orientations'],p['transport']['methods']):
        refinement('time',(grid,fine_dt,orientation,method),(grid,coarse_dt,orientation,method))
    for dt,orientation,method in itertools.product(p['transport']['dt_over_t_dyn'],p['launch']['orientations'],p['transport']['methods']):
        refinement('grid',(fine_grid,dt,orientation,method),(coarse_grid,dt,orientation,method))
    kernels={}
    for grid in p['transport']['grids']:
        path=OUT/f'kernel-{grid}.json';k=json.loads(path.read_text())
        inputs[str(path.relative_to(ROOT))]=digest(path)
        kernels[grid]={key:k[key] for key in ['shape','vertex_count','accepted_vertices','all_vertex_quadrature_passed','all_holdout_quadrature_passed','holdout_error_maxima','scales']}
    for f in ['analysis/matched/analyze.py','analysis/matched/statistics.py']:
        inputs[f]=digest(ROOT/f)
    sampling_design_path=OUT/'sampling-validation-design.json'
    sampling_result_path=OUT/'sampling-validation.json'
    sampling_design=json.loads(sampling_design_path.read_text())
    sampling_result=json.loads(sampling_result_path.read_text())
    if sampling_result['design_sha256']!=digest(sampling_design_path):
        raise ValueError('Stale sampling validation design')
    for path in (sampling_design_path,sampling_result_path):
        inputs[str(path.relative_to(ROOT))]=digest(path)
    sampling_summary={
        'revision':sampling_design['revision'],
        'states_per_grid':len(sampling_design['states']),
        'grids':sampling_design['grids'],
        'samples_per_state_per_grid':sampling_design['samples_per_state_per_grid'],
        'total_sampled_increments':len(sampling_design['states'])*len(sampling_design['grids'])*sampling_design['samples_per_state_per_grid'],
        'all_checks_passed':sampling_result['all_checks_passed'],
        'interpretation':sampling_design['interpretation']}
    return {'revision':p['revision'],'cluster':p['cluster'],'collisional_trajectories_performed':True,
            'matched_encounter_comparison_performed':True,'cluster_residence_response_evaluated':False,
            'ensembles':len(data),'paths_per_ensemble':p['transport']['paths_per_setting'],
            'total_trajectories':len(data)*p['transport']['paths_per_setting'],
            'scope':p['scope'],'reporting':p['reporting'],'kernels':kernels,
            'sampling_validation':sampling_summary,
            'settings':records,'comparisons':comparisons,'refinements':refinements,
            'input_sha256':inputs,
            'interpretation':'Descriptive independent-ensemble intervals, not simultaneous coverage. Finite grid/time changes are sensitivity tests, not rigorous truncation bounds. Numerical-domain stops are censoring, not physical escape. The response predictor has not been evaluated.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=analyze();path=OUT/'summary.json'
    if args.verify:
        stored=json.loads(path.read_text())
        # Numeric checks are tolerance based for platform-independent reductions.
        from analysis.matched.verify import compare_tree
        compare_tree(stored,result)
        print('Verified matched ensemble summaries from retained samples; no trajectories rerun.')
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print('Summarized',result['total_trajectories'],'conditional trajectories.')
if __name__=='__main__':main()
