"""Recompute the independent response retest from retained arrays; no dynamics."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from analysis.retest.statistics import mean_interval,residual,contrast,classify,continuation_upper,tail_diagnostics
OUT=ROOT/'results/retest'
PROTOCOL=ROOT/'analysis/retest/protocol.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):return json.loads(Path(path).read_text())


def analyze():
    p=read(PROTOCOL);report=p['reporting'];confidence=report['normal_confidence'];family=report['normal_family_size']
    delta=report['boundary_family_alpha']/report['boundary_family_size']
    inputs={str(f.relative_to(ROOT)):digest(f) for f in [PROTOCOL,Path(__file__),
        ROOT/'analysis/retest/statistics.py',OUT/'execution-design.json',OUT/'prediction-register.json']}
    predictions=[];validations=[];bounds=[];samples={};records=[]
    def load(stage,setting,launch,a=None):
        suffix='weighted' if stage=='prediction' else f'a{a:.3f}'
        name=f"{stage}-{setting['name']}-{launch['name']}-{suffix}"
        path=OUT/(name+'.npz');meta=path.with_suffix('.json');r=read(meta)
        if not r['complete'] or digest(path)!=r['samples_sha256']:raise ValueError('Changed outcome '+name)
        for f in (path,meta):inputs[str(f.relative_to(ROOT))]=digest(f)
        with np.load(path,allow_pickle=False) as data:arrays=dict(data)
        samples[(stage,setting['name'],launch['name'],a)]=arrays;records.append(r)
        maximum=2*r['n_steps'] if stage=='prediction' else 2.
        bound=continuation_upper(arrays['numerical_missing_bound'],maximum,delta=delta)
        bound.update({'name':name,'budget':report['numerical_missing_residence_budget'],
                      'population_budget_met':bound['upper_bound']<=report['numerical_missing_residence_budget']})
        bounds.append(bound)
        return name,r,arrays,bound
    for setting in p['numerical_settings']:
        for launch in p['launches']:
            name,r,a,b=load('prediction',setting,launch)
            stats=mean_interval(a['response_sample'],confidence,family)
            predictions.append({'name':name,'setting':setting['name'],'launch':launch['name'],'role':setting['role'],
                'grid':setting['grid'],'dt':setting['dt'],'response':stats,
                'precision_target_met':stats['half_width']<=p['prediction']['halfwidth_target'],
                'tails':tail_diagnostics(a['response_sample']),'continuation':b,'stops':r['stops']})
    for setting in p['numerical_settings']:
        values=p['validation']['primary_a_values'] if setting['role']=='primary' else p['validation']['control_a_values']
        for launch in p['launches']:
            pred=next(r for r in predictions if r['setting']==setting['name'] and r['launch']==launch['name'])
            scores=samples[('prediction',setting['name'],launch['name'],None)]['response_sample']
            for value in values:
                name,r,a,b=load('validation',setting,launch,value)
                shift=mean_interval(a['paired_difference'],confidence,family)
                error=residual(a['paired_difference'],scores,value,confidence,family)
                allowance=b['upper_bound']+value*pred['continuation']['upper_bound']
                envelope=[error['interval'][0]-allowance,error['interval'][1]+allowance]
                validations.append({'name':name,'setting':setting['name'],'launch':launch['name'],'role':setting['role'],
                    'a':value,'observed_shift':shift,'prediction_error':error,
                    'scaled_remainder':{'estimate':error['estimate']/value,'interval':[x/value for x in error['interval']],
                        'scope':'Same residual event rescaled; finite-a comparison is not the exact derivative'},
                    'finite_model_classification':classify(error['interval'],report['prediction_absolute_error_target']) if value==1 else None,
                    'continuation':b,'residual_continuation_allowance':allowance,
                    'continuation_adjusted_residual_interval':envelope,
                    'continuation_adjusted_classification':classify(envelope,report['prediction_absolute_error_target']) if value==1 else None,
                    'continuation_scope':'Enlarged normal interval using simultaneous population allowances for bounded-reward continuations of the same stopped prefixes; no model specified beyond the numerical domain',
                    'tails':tail_diagnostics(a['paired_difference']),'stops':r['stops']})
    refinements=[];coarse,fine,primary=p['numerical_settings']
    for launch in p['launches']:
        for kind,first,second in [('grid',fine,coarse),('time',primary,fine)]:
            for observable in ('response','endpoint_shift'):
                stage='prediction' if observable=='response' else 'validation'
                value=None if stage=='prediction' else 1.
                key='response_sample' if stage=='prediction' else 'paired_difference'
                a,b=[samples[(stage,s['name'],launch['name'],value)][key] for s in (first,second)]
                stats=contrast(a,b,confidence,family)
                refinements.append({'kind':kind,'launch':launch['name'],'observable':observable,
                    'first':first['name'],'second':second['name'],'difference':stats,
                    'classification':classify(stats['interval'],report['numerical_sensitivity_budget']),
                    'scope':'Independent finite-resolution sensitivity, not a truncation-error proof'})
    count=len(predictions)+2*len(validations)+len(refinements)
    if count!=report['inferential_statements'] or len(bounds)!=report['boundary_family_size']:
        raise ValueError('Incorrect inferential family inventory')
    return {'revision':p['revision'],'scope':p['scope'],'complete':True,'reporting':report,
        'counts':{'prediction_ensembles':len(predictions),'validation_ensembles':len(validations),
            'prediction_quartets':sum(r['job']['paths'] for r in records if r['job']['stage']=='prediction'),
            'validation_pairs':sum(r['job']['paths'] for r in records if r['job']['stage']=='validation'),
            'independent_groups':sum(r['job']['paths'] for r in records),
            'branch_trajectories':sum(r['job']['paths']*len(r['branches']) for r in records),
            'normal_statements':count,'continuation_bounds':len(bounds)},
        'predictions':predictions,'validations':validations,'refinements':refinements,
        'continuation_bounds':bounds,'continuum_accuracy_classification':'unresolved',
        'interpretation':'Independent test of the numerically stopped finite-step response only. Sampling tails, resolution sensitivity, numerical continuation and physical calibration remain distinct.',
        'input_sha256':inputs}


if __name__=='__main__':
    result=analyze();(OUT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Analyzed independent response retest:',result['counts'])
