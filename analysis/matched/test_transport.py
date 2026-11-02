"""Independent analytical limits for matched collisional transport."""
import numpy as np
import pytest

from analysis.matched.transport import run_ensemble, sphere_residence_fraction


class ConstantLaw:
    def __init__(self,rate=0.,increment=(1.,0.,0.),excluded=0.,D=None):
        self.rate,self.increment,self.excluded = rate,np.asarray(increment),excluded
        self.D = rate*np.outer(self.increment,self.increment) if D is None else np.asarray(D)

    def moments(self,states):
        n=len(states)
        return {'rate':np.full(n,self.rate),'excluded':np.full(n,self.excluded),
                'A':np.tile(self.rate*self.increment,(n,1)),
                'D':np.tile(self.D,(n,1,1))}

    def sample_events(self,states,rng):
        return np.tile(self.increment,(len(states),1))


def zero_force(x):
    return np.zeros_like(x)


def test_chord_occupation_inside_outside_tangent_and_stationary():
    initial=np.array([[-2.,0,0],[0,0,0],[-2,1,0],[1,0,0],[-1,0,0],[2,0,0]])
    final=np.array([[2.,0,0],[.5,0,0],[2,1,0],[1,0,0],[1,0,0],[3,0,0]])
    np.testing.assert_allclose(sphere_residence_fraction(initial,final),[.5,1.,0.,0.,1.,0.],atol=1e-14)


@pytest.mark.parametrize('method',['diffusion','jump'])
def test_compound_poisson_and_gaussian_exact_frozen_moments(method):
    n=80000
    out=run_ensemble(ConstantLaw(rate=4.),zero_force,np.zeros((n,6)),1.,1.,
                     seed=2026103001,method=method,spatial_stop=1e6,velocity_stop=1e6)
    velocity=out['states'][:,3]
    # For Poisson(4), E[N]=4, E[N^2]=20 and E[N^4]=756. The Gaussian
    # with matching mean/variance has raw fourth moment 688.
    assert abs(velocity.mean()-4.) < 6*np.sqrt(4/n)
    fourth=756 if method=='jump' else 688
    assert abs(np.mean(velocity**2)-20.) < 6*np.sqrt((fourth-400)/n)
    np.testing.assert_array_equal(out['states'][:,4:],0.)
    assert out['alive'].all()
    if method=='jump':
        assert out['elastic_events'].max() > 1
        np.testing.assert_array_equal(velocity,out['elastic_events'])


def test_harmonic_null_control_and_chord_time_refinement():
    initial=np.array([[2.,0,0,0,0,0]])
    expected=2*np.pi/3
    errors=[]
    for dt in [.04,.02,.01]:
        out=run_ensemble(ConstantLaw(),lambda x:-x,initial,2*np.pi,dt,7,'diffusion',
                         spatial_stop=3.,velocity_stop=3.)
        errors.append(abs(out['residence'][0]-expected))
        assert out['stop_cause'][0]=='horizon'
    assert errors[-1] < 3e-5
    assert errors[-1] < errors[0]


@pytest.mark.parametrize('method',['diffusion','jump'])
def test_exact_excluded_clock_and_stopped_occupation(method):
    n=40000
    out=run_ensemble(ConstantLaw(excluded=2.),zero_force,np.zeros((n,6)),1.,.1,
                     2026103002,method,spatial_stop=2.,velocity_stop=2.,
                     density=lambda x:np.full(len(x),2.))
    expected=(1-np.exp(-2))/2
    second=(1-3*np.exp(-2))/2
    assert abs(out['residence'].mean()-expected) < 6*np.sqrt((second-expected**2)/n)
    np.testing.assert_allclose(out['residence'],out['end_time'],atol=1e-14)
    np.testing.assert_allclose(out['exposure'],2*out['residence'],atol=1e-14)
    excluded=out['stop_cause']=='excluded_channel'
    assert abs(excluded.mean()-(1-np.exp(-2))) < 6*np.sqrt(np.exp(-2)*(1-np.exp(-2))/n)
    assert np.all(out['alive']==~excluded)
    assert not out['elastic_events'].any()


def test_spatial_stop_is_first_chord_exit_and_initial_domain_checks():
    initial=np.array([[0.,0,0,1.,0,0],[2.,0,0,0,0,0],[0,0,0,4.,0,0]])
    out=run_ensemble(ConstantLaw(),zero_force,initial,5.,.7,8,'jump',
                     spatial_stop=2.,velocity_stop=3.)
    np.testing.assert_allclose(out['end_time'],[2.,0.,0.],atol=1e-13)
    np.testing.assert_allclose(out['residence'],[1.,0.,0.],atol=1e-13)
    assert list(out['stop_cause'])==['spatial_domain','spatial_domain','velocity_domain']
    np.testing.assert_allclose(out['states'][0,:3],[2.,0,0],atol=1e-13)
    assert not out['alive'].any()


def test_velocity_exit_at_midpoint_censors_future_residence():
    # A deterministic Gaussian drift supplies exactly one known increment.
    law=ConstantLaw(rate=1.,increment=(2.,0,0),D=np.zeros((3,3)))
    out=run_ensemble(law,zero_force,np.zeros((1,6)),1.,1.,9,'diffusion',
                     spatial_stop=4.,velocity_stop=1.)
    assert out['stop_cause'][0]=='velocity_domain'
    assert out['end_time'][0]==out['residence'][0]==.5
    np.testing.assert_array_equal(out['states'][0,:3],0.)


def test_reproducibility_and_no_input_mutation():
    initial=np.zeros((100,6));saved=initial.copy()
    law=ConstantLaw(rate=2.,excluded=.3)
    first=run_ensemble(law,zero_force,initial,1.,.1,10,'jump',spatial_stop=20.,velocity_stop=20.)
    second=run_ensemble(law,zero_force,initial,1.,.1,10,'jump',spatial_stop=20.,velocity_stop=20.)
    for key in ('states','residence','end_time','stop_cause','elastic_events'):
        np.testing.assert_array_equal(first[key],second[key])
    np.testing.assert_array_equal(initial,saved)


def test_non_psd_diffusion_is_rejected():
    law=ConstantLaw(D=np.diag([1.,-.001,1.]))
    with pytest.raises(ValueError,match='positive semidefinite'):
        run_ensemble(law,zero_force,np.zeros((1,6)),1.,.1,4,'diffusion')


@pytest.mark.parametrize('method',['diffusion','jump'])
def test_anisotropic_raw_moments_with_nonzero_initial_velocity(method):
    # A genuinely three-dimensional positive jump measure checks rotations and
    # off-diagonal D terms. D is the raw second jump moment, not the central
    # variance of a single jump; Poisson fluctuations supply the extra term.
    jumps=np.array([[1.,.3,-.2],[-.4,.8,.1],[.2,-.7,.9],[-.2,-.1,-.3]])
    weights=np.array([.4,.3,.2,.1])
    rate,horizon,n=3.,.4,50000
    class DiscreteLaw:
        def moments(self,states):
            size=len(states)
            A=rate*np.einsum('i,ij->j',weights,jumps)
            D=rate*np.einsum('i,ij,ik->jk',weights,jumps,jumps)
            return {'rate':np.full(size,rate),'excluded':np.zeros(size),
                    'A':np.tile(A,(size,1)),'D':np.tile(D,(size,1,1))}

        def sample_events(self,states,rng):
            return jumps[rng.choice(len(weights),size=len(states),p=weights)]

    initial=np.zeros((n,6)); initial[:,3:]=[.3,-.2,.1]
    out=run_ensemble(DiscreteLaw(),zero_force,initial,horizon,horizon,
                     2026110101,method,spatial_stop=1e6,velocity_stop=1e6)
    projections=np.concatenate((np.eye(3),[[1,1,0],[1,0,1],[0,1,1]]))
    for direction in projections:
        z=jumps@direction
        mean=initial[0,3:]@direction+horizon*rate*(weights@z)
        variance=horizon*rate*(weights@z**2)
        third=horizon*rate*(weights@z**3) if method=='jump' else 0.
        fourth=horizon*rate*(weights@z**4) if method=='jump' else 0.
        raw_second=mean**2+variance
        raw_fourth=mean**4+6*mean**2*variance+4*mean*third+3*variance**2+fourth
        observed=out['states'][:,3:]@direction
        assert abs(observed.mean()-mean) < 6*np.sqrt(variance/n)
        assert abs(np.mean(observed**2)-raw_second) < 6*np.sqrt((raw_fourth-raw_second**2)/n)
    # Both half-drifts bracket exactly one velocity update in this frozen step.
    expected_position=.5*horizon*(initial[:,3:]+out['states'][:,3:])
    np.testing.assert_allclose(out['states'][:,:3],expected_position,atol=1e-14)


def test_competing_spatial_and_excluded_stops_preserve_weighted_accounting():
    n,horizon,outer,rate=40000,2.,.6,2.
    initial=np.zeros((n,6)); initial[:,3]=1.
    out=run_ensemble(ConstantLaw(excluded=rate),zero_force,initial,horizon,.23,
                     2026110102,'jump',spatial_stop=outer,velocity_stop=3.,
                     density=lambda x:1+x[:,0],region_radius=.4)
    elapsed=out['end_time']
    np.testing.assert_allclose(out['states'][:,0],elapsed,atol=1e-14)
    np.testing.assert_allclose(out['residence'],np.minimum(elapsed,.4),atol=1e-14)
    # Midpoint integration is exact for this linear density, including the
    # shortened boundary drift, and must exclude every post-stop contribution.
    np.testing.assert_allclose(out['exposure'],elapsed+.5*elapsed**2,atol=1e-14)
    np.testing.assert_allclose(out['normalized_residence'],out['residence']/horizon)
    np.testing.assert_allclose(out['normalized_exposure'],out['exposure']/horizon)
    spatial=out['stop_cause']=='spatial_domain'
    excluded=out['stop_cause']=='excluded_channel'
    assert np.all(spatial ^ excluded)
    assert not out['alive'].any()
    assert not out['elastic_events'].any()
    survival=np.exp(-rate*outer)
    assert abs(spatial.mean()-survival) < 6*np.sqrt(survival*(1-survival)/n)
    mean=(1-survival)/rate
    second=2*(1-survival*(1+rate*outer))/rate**2
    assert abs(elapsed.mean()-mean) < 6*np.sqrt((second-mean**2)/n)
