"""Independent time integration, analytic limits and finite-window checks."""
import math
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.special import gamma
from recoil_dynamics import PlummerOrbit, core_exit_threshold


def integrate_reference(eta, horizon, core_radius=1.):
    def rhs(t, y):
        x, v, exposure = y
        return [v, -x/(1+x*x)**1.5, (1+x*x)**(-2.5)]
    def central(t, y):
        return y[0]
    central.direction = -1
    def core_crossing(t, y):
        return abs(y[0]) - core_radius
    result = solve_ivp(rhs, [0,horizon], [0,math.sqrt(2)*eta,0],
                       method='DOP853', rtol=2e-11, atol=2e-12,
                       dense_output=True, max_step=.07, events=[central,core_crossing])
    assert result.success
    edges = np.r_[0,result.t_events[1],horizon]
    core_time = 0.
    for a,b in zip(edges[:-1],edges[1:]):
        if abs(result.sol((a+b)/2)[0]) < core_radius:
            core_time += b-a
    return result, core_time


@pytest.mark.parametrize('eta',[.2,.5,.8,.95])
def test_first_return_against_ode_and_conserved_energy(eta):
    orb=PlummerOrbit(eta)
    sol,_=integrate_reference(eta,orb.period*1.1)
    assert sol.t_events[0][0] == pytest.approx(orb.period,rel=3e-8)
    energy=.5*sol.y[1]**2-1/np.sqrt(1+sol.y[0]**2)
    assert np.max(np.abs(energy-(eta*eta-1))) < 2e-9


@pytest.mark.parametrize('eta,horizon,core',[(.3,1.7,.1),(.3,25.,.3),(.8,7.,.3),
                                           (.8,40.,1.),(.99,100.,1.)])
def test_finite_exposure_and_event_based_core_time_against_ode(eta,horizon,core):
    orb=PlummerOrbit(eta)
    sol,core_time=integrate_reference(eta,horizon,core)
    assert orb.exposure(horizon)==pytest.approx(sol.y[2,-1],rel=2e-8,abs=2e-9)
    assert orb.core_time(horizon,core)==pytest.approx(core_time,rel=2e-8,abs=2e-9)


def test_stationary_and_zero_window():
    orb=PlummerOrbit(0)
    assert orb.period is None and orb.apocenter==0
    for h in [0,1,1000]:
        assert orb.exposure(h)==h
        assert orb.core_time(h,.1)==h
    with pytest.raises(ValueError): orb.cycle_exposure()
    with pytest.raises(ValueError): orb.cycle_core_time(1)
    assert PlummerOrbit(.7).exposure(0)==0
    assert PlummerOrbit(.7).core_time(0,1)==0


@pytest.mark.parametrize('core',[.1,.3,1.])
def test_core_exit_threshold_and_all_inside(core):
    eta=core_exit_threshold(core)
    assert PlummerOrbit(eta).apocenter==pytest.approx(core,rel=1e-13)
    orb=PlummerOrbit(eta*(1-1e-5))
    assert orb.core_time(17.3,core)==17.3
    outside=PlummerOrbit(eta*(1+1e-5))
    assert outside.cycle_core_time(core)<outside.period


def test_analytic_harmonic_and_parabolic_limits():
    small=PlummerOrbit(1e-5)
    assert small.period==pytest.approx(math.pi,rel=2e-10)
    assert small.cycle_exposure()/small.period==pytest.approx(1,rel=3e-10)
    large=PlummerOrbit(.9999)
    asymptote=math.pi/math.sqrt(2)*(1-large.eta**2)**(-1.5)
    assert large.period/asymptote==pytest.approx(1,rel=2e-5)
    exposure_limit=math.sqrt(math.pi/2)*gamma(7/4)/gamma(9/4)
    assert large.cycle_exposure()==pytest.approx(exposure_limit,rel=2e-4)
    density_limit=gamma(7/4)/(math.sqrt(math.pi)*gamma(9/4))
    actual=large.cycle_exposure()/large.period/(1-large.eta**2)**1.5
    assert actual==pytest.approx(density_limit,rel=2e-4)


@pytest.mark.parametrize('eta',[.1,.8,.99])
def test_period_accounting_bounds_and_continuity(eta):
    orb=PlummerOrbit(eta)
    for n in [1,2,10]:
        assert orb.exposure(n*orb.period)==pytest.approx(n*orb.cycle_exposure(),rel=2e-12)
        assert orb.core_time(n*orb.period,.3)==pytest.approx(n*orb.cycle_core_time(.3),rel=2e-12)
    for h in [1,3,10,100,1000]:
        assert 0<=orb.exposure(h)<=h*(1+1e-12)
        assert 0<=orb.core_time(h,.3)<=h*(1+1e-12)
    # Both integrals are 1-Lipschitz in the horizon, including crossing times.
    tc=orb.cycle_core_time(.3)/2
    eps=1e-6
    for h in [tc,orb.period/2,orb.period-tc,orb.period,2*orb.period]:
        for fn in [orb.exposure,lambda t:orb.core_time(t,.3)]:
            delta=fn(h+eps)-fn(max(0,h-eps))
            assert -1e-10<=delta<=2*eps+1e-9


def test_exposure_time_derivative_against_orbit_density():
    orb=PlummerOrbit(.8)
    sol,_=integrate_reference(.8,7)
    for h in [1.,4.,6.]:
        dx=1e-5
        derivative=(orb.exposure(h+dx)-orb.exposure(h-dx))/(2*dx)
        density=(1+sol.sol(h)[0]**2)**(-2.5)
        assert derivative==pytest.approx(density,rel=1e-6,abs=2e-9)


@pytest.mark.parametrize('eta',[.99,.9999])
def test_quadrature_tolerance_convergence(eta):
    base=PlummerOrbit(eta)
    tight=PlummerOrbit(eta,epsabs=1e-13,epsrel=1e-13)
    assert base.period==pytest.approx(tight.period,rel=1e-10)
    assert base.cycle_exposure()==pytest.approx(tight.cycle_exposure(),rel=1e-10)
    for h in [1.,100.]:
        assert base.exposure(h)==pytest.approx(tight.exposure(h),rel=1e-9,abs=1e-10)


def test_first_return_comparator_changes_error_sign():
    orb=PlummerOrbit(.99)
    early=orb.period/2
    late=10*orb.period
    assert max(early-orb.period,0)-orb.exposure(early)<0
    assert max(late-orb.period,0)-orb.exposure(late)>0


@pytest.mark.parametrize('eta',[-.1,1.,2.,float('nan'),float('inf')])
def test_invalid_eta(eta):
    with pytest.raises(ValueError): PlummerOrbit(eta)


def test_invalid_horizons_radii_and_tolerances():
    orb=PlummerOrbit(.5)
    for h in [-1,float('nan'),float('inf')]:
        with pytest.raises(ValueError): orb.exposure(h)
        with pytest.raises(ValueError): orb.core_time(h,1)
    for c in [0,-1,float('nan'),float('inf')]:
        with pytest.raises(ValueError): orb.core_time(1,c)
        with pytest.raises(ValueError): core_exit_threshold(c)
    with pytest.raises(ValueError): PlummerOrbit(.5,epsabs=0)
    with pytest.raises(ValueError): PlummerOrbit(.5,epsrel=-1)
