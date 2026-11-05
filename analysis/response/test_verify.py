"""Evidence-verifier rejection controls, without cluster simulations."""
import numpy as np
import pytest

from analysis.response.verify import (
    BRANCHES,PROTOCOL,count_selection,expected_jobs,read,summary,verify_sample_arrays,
    timestamp,verify_adapter_record,
)


def fixture():
    n=4;a={};branches=BRANCHES['weighted']
    for branch in branches:
        b={'states':np.zeros((n,6)),'residence':np.ones(n),
           'normalized_residence':np.ones(n),'exposure':np.ones(n),
           'normalized_exposure':np.ones(n),'end_time':np.ones(n),
           'stop_cause':np.full(n,'horizon',dtype='U24'),'alive':np.ones(n,dtype=bool),
           'elastic_events':np.zeros(n,dtype=np.int64),'evaluated_steps':np.full(n,2,dtype=np.int64)}
        a.update({f'{branch}__{key}':value for key,value in b.items()})
    a.update({'numerical_missing_bound':np.zeros(n),'response_sample':np.zeros(n),
        'replacement_steps':np.array([0,1,0,1]),'active_at_replacement':np.ones(n,dtype=bool),
        'replacement_states':np.zeros((n,6)),'jump_zero_weight':np.ones(n),
        'jump_positive_weight':np.zeros(n),'selected_poisson_mean':np.zeros(n),
        'positive_branch_count':np.zeros(n,dtype=np.int64)})
    r={'job':{'method':'weighted','paths':n},'branches':list(branches),'horizon':1.,'n_steps':2,
       'spatial_stop':4.,'velocity_stop':4.,'stops':{branch:{'horizon':n} for branch in branches},
       'summary':summary(a['response_sample']),'numerical_missing_bound':summary(a['numerical_missing_bound']),
       'array_inventory':{key:{'shape':list(value.shape),'dtype':str(value.dtype)} for key,value in a.items()}}
    return a,r


def test_valid_zero_response_accounting_passes():
    a,r=fixture();verify_sample_arrays(a,r)


def test_adjusting_metadata_cannot_hide_corrupted_response_score():
    a,r=fixture();a['response_sample'][0]=.2;r['summary']=summary(a['response_sample'])
    with pytest.raises(AssertionError):verify_sample_arrays(a,r)


def test_invalid_surviving_duration_is_rejected():
    a,r=fixture();a['jump_zero__end_time'][0]=.5
    with pytest.raises(AssertionError,match='wrong horizon'):verify_sample_arrays(a,r)


def test_numerical_stop_requires_missing_residence_allowance():
    a,r=fixture();name='jump_zero'
    for key in ['residence','normalized_residence','exposure','normalized_exposure','end_time']:
        a[f'{name}__{key}'][0]=.5
    a[f'{name}__alive'][0]=False;a[f'{name}__stop_cause'][0]='spatial_domain'
    a[f'{name}__states'][0,0]=4.
    r['stops'][name]={'horizon':3,'spatial_domain':1}
    a['response_sample'][0]=-1.;r['summary']=summary(a['response_sample'])
    with pytest.raises(AssertionError):verify_sample_arrays(a,r)
    a['numerical_missing_bound'][0]=1.;r['numerical_missing_bound']=summary(a['numerical_missing_bound'])
    verify_sample_arrays(a,r)


def test_zero_weight_branch_cannot_claim_positive_poisson_count():
    a,r=fixture();a['positive_branch_count'][0]=1
    with pytest.raises(AssertionError,match='conditional positive counts'):verify_sample_arrays(a,r)


def test_registered_inventory_has_disjoint_stages_and_variants():
    p=read(PROTOCOL);counts={launch['name']:8192 for launch in p['launches']}
    groups=[list(expected_jobs(p,stage,counts)) for stage in ['pilot','prediction','validation']]
    assert [len(group) for group in groups]==[6,9,21]
    seeds=[tuple(job['seed']) for group in groups for job in group]
    assert len(set(seeds))==36
    assert {seed[1] for seed in seeds}=={10,20,30}


def test_sample_count_rule_caps_without_making_precision_claim():
    p=read(PROTOCOL)
    low,_=count_selection(0.,p)
    high,_=count_selection(1.,p)
    assert low['selected_count']==8192 and not low['cap_hit']
    assert high['selected_count']==65536 and high['cap_hit']
    assert high['uncapped_variance_rule_count']>high['selected_count']


def test_storage_adapter_rejects_unregistered_usage():
    record={'storage_adapter':{'path':'results/response/storage-adapter.json','sha256':'abc'}}
    with pytest.raises(AssertionError,match='unregistered adapter usage'):
        verify_adapter_record(None,'prediction-example',record)


@pytest.mark.parametrize('failure',['stale_hash','late_manifest','missing_diagnostics'])
def test_storage_adapter_requires_prior_identified_manifest(failure):
    reference={'path':'results/response/storage-adapter.json','sha256':'abc'}
    adapter={'adapted':{'prediction-example'},'created':timestamp('2026-11-04T10:00:00+00:00'),
             'reference':reference}
    record={'storage_adapter':dict(reference),'disk_cache':{},'started_utc':'2026-11-04T11:00:00+00:00'}
    verify_adapter_record(adapter,'prediction-example',record)
    if failure=='stale_hash':record['storage_adapter']['sha256']='different'
    elif failure=='late_manifest':record['started_utc']='2026-11-04T09:00:00+00:00'
    else:record.pop('disk_cache')
    with pytest.raises(AssertionError):verify_adapter_record(adapter,'prediction-example',record)
