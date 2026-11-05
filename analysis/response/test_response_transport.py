"""Transition-law, stopping, and paired-response controls independent of hosts."""
import numpy as np
import pytest

from analysis.matched.transport import run_ensemble
from analysis.response.transport import (
    advance_step, initialize_batch, run_mixture, run_paired_mixture,
    run_paired_response, run_weighted_response, sample_positive_poisson,
)

ARRAYS=('states','residence','exposure','end_time','stop_cause','alive',
        'elastic_events','evaluated_steps','normalized_residence','normalized_exposure')


class ConstantLaw:
    def __init__(self,rate=0.,increment=(.3,-.1,.2),excluded=0.):
        self.rate=rate;self.increment=np.array(increment);self.excluded=excluded

    def moments(self,states):
        n=len(states)
        return {'rate':np.full(n,self.rate),'excluded':np.full(n,self.excluded),
                'A':np.tile(self.rate*self.increment,(n,1)),
                'D':np.tile(self.rate*np.outer(self.increment,self.increment),(n,1,1))}

    def sample_events(self,states,rng):
        return np.tile(self.increment,(len(states),1))


class StateLaw(ConstantLaw):
    def moments(self,states):
        value=super().moments(states)
        scale=1+np.linalg.norm(states[:,:3],axis=1)
        value['rate']*=scale;value['excluded']*=scale
        value['A']*=scale[:,None];value['D']*=scale[:,None,None]
        return value


def zero_force(x):
    return np.zeros_like(x)


def assert_same_arrays(a,b):
    for key in ARRAYS:
        np.testing.assert_array_equal(a[key],b[key],err_msg=key)


@pytest.mark.parametrize('a,method',[(0.,'diffusion'),(1.,'jump')])
def test_mixture_endpoints_bitwise_match_original_complete_solver(a,method):
    initial=np.random.default_rng(24).normal(scale=.25,size=(128,6))
    initial[0,0]=4.;initial[1,3]=5.
    law=StateLaw(rate=3.,excluded=.4)
    options=dict(spatial_stop=1.2,velocity_stop=2.,density=lambda x:np.exp(-np.sum(x*x,axis=1)))
    old=run_ensemble(law,lambda x:-x,initial,1.03,.1,1234,method,**options)
    new=run_mixture(law,lambda x:-x,initial,1.03,.1,1234,a,**options)
    assert_same_arrays(old,new)


def test_paired_baseline_and_mixture_marginals_match_standalone_bitwise():
    initial=np.zeros((256,6));initial[:,0]=.2
    law=StateLaw(rate=4.,excluded=.3)
    options=dict(spatial_stop=2.,velocity_stop=2.,density=lambda x:1+np.sum(x*x,axis=1))
    pair=run_paired_mixture(law,lambda x:-x,initial,1.,.1,1235,.3,**options)
    old=run_ensemble(law,lambda x:-x,initial,1.,.1,1235,'diffusion',**options)
    alone=run_mixture(law,lambda x:-x,initial,1.,.1,1235,.3,**options)
    assert_same_arrays(pair['baseline'],old)
    assert_same_arrays(pair['mixture'],alone)


def test_single_replacement_baseline_is_exact_old_diffusion():
    initial=np.zeros((128,6));initial[:,0]=.2
    law=StateLaw(rate=2.,excluded=.2)
    pair=run_paired_response(law,lambda x:-x,initial,1.03,.1,1236)
    old=run_ensemble(law,lambda x:-x,initial,1.03,.1,1236,'diffusion')
    assert_same_arrays(pair['baseline'],old)
    assert pair['n_steps']==11
    assert np.all((pair['replacement_steps']>=0)&(pair['replacement_steps']<11))


def test_complete_step_replacement_includes_within_step_rewards_and_stops():
    law=ConstantLaw(rate=1.,increment=(2.,0,0))
    initial=np.zeros((1,6))
    diffusion=initialize_batch(initial,1.,4.,4.)
    no_event=initialize_batch(initial,1.,4.,4.)
    normal=np.zeros((1,3));exponential=np.ones(1)*100
    count_rng=np.random.default_rng(1);event_rng=np.random.default_rng(2)
    common=dict(spatial_stop=4.,velocity_stop=4.,region_radius=.25,density=lambda x:1+x[:,0])
    advance_step(law,zero_force,diffusion,0.,1.,normal,exponential,np.array([False]),count_rng,event_rng,**common)
    advance_step(law,zero_force,no_event,0.,1.,normal,exponential,np.array([True]),count_rng,event_rng,
                 forced_counts=np.zeros(1,dtype=int),**common)
    np.testing.assert_allclose(diffusion['residence'],.625)
    np.testing.assert_allclose(diffusion['exposure'],1.25)
    np.testing.assert_allclose(no_event['residence'],1.)
    np.testing.assert_allclose(no_event['exposure'],1.)
    stopped=initialize_batch(initial,1.,4.,1.5)
    advance_step(law,zero_force,stopped,0.,1.,normal,exponential,np.array([False]),count_rng,event_rng,
                 spatial_stop=4.,velocity_stop=1.5)
    assert stopped['stop_cause'][0]=='velocity_domain'
    assert stopped['end_time'][0]==stopped['residence'][0]==.5
    saved={key:value.copy() for key,value in stopped.items()}
    advance_step(law,zero_force,stopped,1.,1.,normal,exponential,np.array([True]),count_rng,event_rng)
    for key in saved:
        np.testing.assert_array_equal(stopped[key],saved[key])


def test_replacement_selected_after_stop_has_zero_contribution():
    initial=np.zeros((64,6))
    law=ConstantLaw(rate=0.,excluded=1e6)
    pair=run_paired_response(law,zero_force,initial,1.,.1,1237,
                             replacement_steps=np.full(64,9,dtype=int))
    assert not pair['active_at_replacement'].any()
    assert_same_arrays(pair['baseline'],pair['replacement'])
    weighted=run_weighted_response(law,zero_force,initial,1.,.1,1237,
                                   replacement_steps=np.full(64,9,dtype=int))
    for name in ('jump_zero','jump_positive','diffusion_minus'):
        assert_same_arrays(weighted['diffusion_plus'],weighted[name])
    assert not weighted['active_at_replacement'].any()


@pytest.mark.parametrize('mu',[1e-8,.2,3.,25.])
def test_positive_poisson_conditional_moments(mu):
    n=50000
    values=sample_positive_poisson(np.full(n,mu),np.random.default_rng(1238))
    assert np.all(values>=1)
    probability=-np.expm1(-mu)
    mean=mu/probability
    second=(mu+mu*mu)/probability
    variance=max(second-mean*mean,0.)
    assert abs(values.mean()-mean)<6*np.sqrt(variance/n)+1e-10
    # Independently check the zero/positive decomposition of the raw second
    # moment, including probability weight at rates approaching zero.
    fourth=(mu+7*mu**2+6*mu**3+mu**4)/probability
    assert abs(np.mean(values.astype(float)**2)-second)<6*np.sqrt(max(fourth-second**2,0.)/n)+1e-9


def test_positive_poisson_resolves_exhausted_roundoff_tail():
    class UpperEndpointUniform:
        def random(self,n):
            return np.full(n,np.nextafter(1.,0.))

    # Dense means include cases whose accumulated CDF is >4 eps below one
    # when every remaining mass is too small to change it. A recurrence that
    # only clamps deficits <=4 eps stalls forever for these valid inputs.
    means=np.linspace(.001,19.999,40000)
    values=sample_positive_poisson(means,UpperEndpointUniform())
    assert np.all((values>=1)&(values<100))
    assert np.all(values>means)


def test_zero_rate_weighted_response_is_exact_zero_and_baseline_stream_is_preserved():
    initial=np.zeros((64,6));initial[:,0]=.4;initial[:,4]=.3
    law=ConstantLaw(excluded=.4)
    weighted=run_weighted_response(law,lambda x:-x,initial,1.,.1,1239)
    old=run_ensemble(law,lambda x:-x,initial,1.,.1,1239,'diffusion')
    assert_same_arrays(weighted['diffusion_plus'],old)
    for name in ('jump_zero','jump_positive','diffusion_minus'):
        assert_same_arrays(weighted['diffusion_plus'],weighted[name])
    np.testing.assert_array_equal(weighted['jump_zero_weight'],1.)
    np.testing.assert_array_equal(weighted['jump_positive_weight'],0.)
    np.testing.assert_array_equal(weighted['positive_branch_count'],0)


def test_weighted_response_frozen_velocity_cubic_matches_exact_derivative():
    # A single frozen transition has Gaussian and compound-Poisson moments
    # equal through order two. Their third raw moment differs by rate*xi^3*H.
    # No stopping or potential complicates this independent analytical target.
    n=50000;rate=.8;increment=.7;horizon=.3
    initial=np.zeros((n,6));initial[:,3]=.2
    law=ConstantLaw(rate=rate,increment=(increment,0,0))
    result=run_weighted_response(law,zero_force,initial,horizon,horizon,1240,
                                 spatial_stop=1e6,velocity_stop=1e6)
    w0=result['jump_zero_weight'];wp=result['jump_positive_weight']
    def observable(name):
        return result[name]['states'][:,3]**3
    score=(w0*observable('jump_zero')+wp*observable('jump_positive')
           -.5*(observable('diffusion_plus')+observable('diffusion_minus')))
    expected=rate*horizon*increment**3
    assert abs(score.mean()-expected)<6*score.std(ddof=1)/np.sqrt(n)
    np.testing.assert_allclose(w0+wp,1.,atol=2e-16)
    np.testing.assert_allclose(result['selected_poisson_mean'],rate*horizon)
    # The positive branch includes multiple events, never a Bernoulli shortcut.
    assert result['positive_branch_count'].max()>1


def test_paired_zero_mixture_has_exact_zero_difference():
    initial=np.zeros((128,6))
    pair=run_paired_mixture(StateLaw(rate=2.,excluded=.2),zero_force,initial,1.,.1,1241,0.)
    assert_same_arrays(pair['baseline'],pair['mixture'])


def test_invalid_explicit_replacements_and_probabilities_are_rejected():
    initial=np.zeros((2,6));law=ConstantLaw()
    with pytest.raises(ValueError,match='replacement_steps'):
        run_paired_response(law,zero_force,initial,1.,.1,1,replacement_steps=np.array([0,10]))
    with pytest.raises(ValueError,match='replacement_steps'):
        run_weighted_response(law,zero_force,initial,1.,.1,1,replacement_steps=np.array([0.,1.]))
    with pytest.raises(ValueError,match='mixture_probability'):
        run_mixture(law,zero_force,initial,1.,.1,1,-.1)
