"""Moment and event-law checks independent of the transport solver."""
import json
import numpy as np
import pytest
from analysis.matched.kernel import (MatchedKernel, Host, ROOT, host_scales,
    vertex_law, VertexSampler, symmetry_moments, symmetry_samples)


def test_positive_mixture_uses_rate_weighted_event_vertices():
    axes=[np.array([0.,2.]),np.array([0.,2.]),np.array([0.,1.])]
    rate=np.arange(1.,9.).reshape(2,2,2)
    increment=np.arange(1.,9.)[:,None]*np.array([[.02,0.,-.01]])
    A=rate.reshape(-1,1)*increment
    D=rate.reshape(-1,1,1)*increment[:,:,None]*increment[:,None,:]
    class Sampler:
        def __init__(self,i):self.i=i
        def sample(self,n,rng):return np.tile(increment[self.i],(n,1))
    kernel=MatchedKernel.from_arrays(axes,rate,A,D,np.zeros((2,2,2)),Sampler)
    state=np.array([[0,0,1.,2**-.5,0,2**-.5]])
    expected=kernel.moments(state)
    n=120000
    jumps=kernel.sample_events(np.repeat(state,n,axis=0),np.random.default_rng(717))
    targets=[expected['A'][0]/expected['rate'][0],expected['D'][0]/expected['rate'][0]]
    quantities=[jumps,jumps[:,:,None]*jumps[:,None,:]]
    for sample,target in zip(quantities,targets):
        assert np.all(abs(sample.mean(axis=0)-target)<6*sample.std(axis=0)/np.sqrt(n)+1e-13)
    assert np.linalg.eigvalsh(expected['D']).min()>-1e-14
    with pytest.raises(ValueError,match='outside'):
        kernel.moments(np.array([[0,0,3.,0,0,0]]))


@pytest.mark.parametrize('coordinates',[(0.,0.,0.),(0.,1.,.5),(.5,0.,0.),(.5,1.,0.),(.5,1.,.5),(.5,1.,1.)])
def test_physical_vertex_sampler_matches_raw_moments(coordinates):
    p=json.loads((ROOT/'analysis/matched/protocol.json').read_text())
    e=json.loads((ROOT/p['encounter_protocol']).read_text())
    h=Host(dict(np.load(ROOT/p['host_path'],allow_pickle=False)))
    scales=host_scales(h)
    moments,law=vertex_law(h,coordinates,scales,p,e,28)
    draws=VertexSampler(*law,scales[1]).sample(40000,np.random.default_rng(716))
    expected=[moments['A']/moments['rate']/scales[1],moments['D']/moments['rate']/scales[1]**2]
    for sample,target in zip([draws,draws[:,:,None]*draws[:,None,:]],expected):
        assert np.all(abs(sample.mean(axis=0)-target)<6*sample.std(axis=0)/np.sqrt(len(draws))+1e-13)


def test_coordinate_rotation_and_inward_reflection():
    axes=[np.array([0.,2.]),np.array([0.,2.]),np.array([0.,1.])]
    shape=(2,2,2)
    A=np.tile([0.,0.,-1.],shape+(1,))
    D=np.tile(np.diag([1.,1.,2.]),shape+(1,1))
    k=MatchedKernel.from_arrays(axes,np.ones(shape),A,D,np.zeros(shape))
    states=np.array([[1.,0,0,1.,0,0],[1.,0,0,-1.,0,0],[0,1.,0,0,1.,0]])
    m=k.moments(states)
    np.testing.assert_allclose(m['A'],[[-1,0,0],[1,0,0],[0,-1,0]])
    np.testing.assert_allclose(m['D'][0],np.diag([2,1,1]))
