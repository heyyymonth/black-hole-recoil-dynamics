"""Independent velocity quadrature and time-integrated encounter checks.

The production velocity moments use beta functions. The references below
instead integrate the normalized speed density numerically; selected tests
also integrate polar angle directly. The orbit reference solves the signed
equations of motion with two accumulated encounter integrals as extra states.
It does not use the production orbital energy quadrature or period folding.
"""
import copy
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from recoil_dynamics.encounters import (EncounterOrbit, experiment_rows,
                                       relative_speed_moments)


def speed_density(y):
    return 512 / (7 * math.pi) * y**2 * (1 - y*y)**3.5


def velocity_reference(speed, escape_speed, direct_angle=False):
    """Numerically integrate speed, optionally also the polar angle."""
    ratio = speed / escape_speed
    edges = [0., ratio, 1.] if 0 < ratio < 1 else [0., 1.]

    def at_y(y, inverse):
        target_speed = escape_speed * y
        if direct_angle:
            def angular(theta):
                # This form avoids cancellation for equal speeds/alignment.
                w = math.hypot(speed-target_speed,
                               2 * math.sqrt(speed*target_speed) * math.sin(theta/2))
                weight = .5 * math.sin(theta)
                return weight / w if inverse else weight * w
            moment = quad(angular, 0., math.pi, epsabs=2e-11,
                          epsrel=2e-11, limit=150)[0]
        else:
            hi, lo = max(speed, target_speed), min(speed, target_speed)
            moment = 1 / hi if inverse else hi + lo*lo / (3*hi)
        return speed_density(y) * moment

    return tuple(math.fsum(
        quad(lambda y: at_y(y, inverse), lo, hi,
             epsabs=2e-11, epsrel=2e-11, limit=150)[0]
        for lo, hi in zip(edges[:-1], edges[1:]))
        for inverse in (False, True))


def integrate_reference(eta, horizon):
    def rhs(_time, state):
        x, velocity, _geometric, _focusing = state
        radius_factor = math.sqrt(1 + x*x)
        density = radius_factor**-5
        moments = velocity_reference(abs(velocity), math.sqrt(2/radius_factor))
        return [velocity, -x/radius_factor**3,
                density*moments[0], density*moments[1]]

    solution = solve_ivp(rhs, (0., horizon), [0., math.sqrt(2)*eta, 0., 0.],
                         method="DOP853", rtol=2e-11, atol=2e-12,
                         max_step=.08)
    assert solution.success
    return solution.y[2:, -1]


def test_reference_speed_density_normalization_and_analytic_moments():
    checks = [(0, 1.), (1, 1024/(693*math.pi)),
              (2, .25), (-1, 512/(63*math.pi))]
    for power, expected in checks:
        value = quad(lambda y: speed_density(y)*y**power, 0., 1.,
                     epsabs=1e-12, epsrel=1e-12)[0]
        assert value == pytest.approx(expected, rel=2e-12)


@pytest.mark.parametrize("ratio", [0., .05, .3, .8, .99, 1., 1.7])
def test_beta_moments_against_direct_speed_and_angle_quadrature(ratio):
    escape = 1.37
    expected = velocity_reference(ratio*escape, escape, direct_angle=True)
    actual = relative_speed_moments(ratio*escape, escape)
    np.testing.assert_allclose(actual, expected, rtol=2e-9, atol=2e-11)


@pytest.mark.parametrize("escape", [.2, math.sqrt(2), 3.])
def test_stationary_and_fast_speed_analytic_limits(escape):
    mean, reciprocal = relative_speed_moments(0., escape)
    assert mean == pytest.approx(escape*1024/(693*math.pi), rel=2e-12)
    assert reciprocal == pytest.approx(512/(63*math.pi)/escape, rel=2e-12)
    for speed in [escape, 2*escape, 10*escape]:
        mean, reciprocal = relative_speed_moments(speed, escape)
        assert mean == pytest.approx(speed+escape**2/(12*speed), rel=2e-12)
        assert reciprocal == pytest.approx(1/speed, rel=2e-12)


@pytest.mark.parametrize("eta", [0., .3, .8, .9])
def test_finite_window_against_independent_time_and_velocity_integration(eta):
    expected = integrate_reference(eta, 10.)
    actual = EncounterOrbit(eta).integrals(10.)
    np.testing.assert_allclose(actual, expected, rtol=1e-8, atol=2e-10)


def test_stationary_integrals_and_zero_window():
    expected_rate = np.array([math.sqrt(2)*1024/(693*math.pi),
                              512/(63*math.pi)/math.sqrt(2)])
    for horizon in [0., 1., 10., 100.]:
        np.testing.assert_allclose(EncounterOrbit(0.).integrals(horizon),
                                   horizon*expected_rate, rtol=2e-12, atol=0.)
    for eta in [.3, .8, .9]:
        np.testing.assert_array_equal(EncounterOrbit(eta).integrals(0.), [0., 0.])


@pytest.mark.parametrize("eta", [.3, .8, .9])
def test_period_accumulation_and_partial_cycle_folding(eta):
    encounter = EncounterOrbit(eta)
    period = encounter.orbit.period
    cycle = np.array(encounter.integrals(period))
    for count in [2, 7]:
        np.testing.assert_allclose(encounter.integrals(count*period), count*cycle,
                                   rtol=2e-10, atol=2e-11)
    for fraction in [.19, .73]:
        partial = np.array(encounter.integrals(fraction*period))
        np.testing.assert_allclose(encounter.integrals((3+fraction)*period),
                                   3*cycle+partial, rtol=2e-10, atol=2e-11)
    # Radial density and speed are symmetric about the apocenter.
    np.testing.assert_allclose(encounter.integrals(period/2), cycle/2,
                               rtol=2e-10, atol=2e-11)


@pytest.mark.parametrize("eta", [.3, .8, .9])
def test_quadrature_tolerance_convergence(eta):
    ordinary = EncounterOrbit(eta).integrals(10.)
    tighter = EncounterOrbit(eta, epsabs=1e-12, epsrel=1e-12).integrals(10.)
    np.testing.assert_allclose(ordinary, tighter, rtol=1e-8, atol=2e-10)


@pytest.mark.parametrize("speed,escape", [(-1., 1.), (math.nan, 1.),
                                           (math.inf, 1.), (1., 0.),
                                           (1., -1.), (1., math.nan),
                                           (1., math.inf)])
def test_invalid_velocity_arguments(speed, escape):
    with pytest.raises(ValueError):
        relative_speed_moments(speed, escape)


def test_invalid_horizons():
    encounter = EncounterOrbit(.3)
    for horizon in [-1., math.nan, math.inf]:
        with pytest.raises(ValueError):
            encounter.integrals(horizon)


@pytest.fixture(scope="module")
def declared_experiment():
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root / "experiments/e003.json").read_text())
    return config, experiment_rows(config)


def test_target_count_scales_absolute_depth_not_proxy_comparison(declared_experiment):
    config, (_scales, baseline) = declared_experiment
    changed = copy.deepcopy(config)
    changed["target_binary_count"] *= 3
    _, rows = experiment_rows(changed)
    for actual, original in zip(rows, baseline):
        for key in ["tau_geometric", "tau_focusing", "tau_total",
                    "tau_stationary", "density_proxy_tau"]:
            assert actual[key] == pytest.approx(3*original[key], rel=2e-12)
        for key in ["encounter_ratio", "proxy_relative_error", "focus_fraction"]:
            assert actual[key] == pytest.approx(original[key], rel=2e-12, abs=2e-14)


def test_passage_radius_has_distinct_geometric_and_focusing_scalings(declared_experiment):
    config, (_scales, baseline) = declared_experiment
    changed = copy.deepcopy(config)
    changed["passage_radius_au"] *= 2
    _, rows = experiment_rows(changed)
    for actual, original in zip(rows, baseline):
        assert actual["tau_geometric"] == pytest.approx(4*original["tau_geometric"], rel=2e-12)
        assert actual["tau_focusing"] == pytest.approx(2*original["tau_focusing"], rel=2e-12)


def test_total_encounter_mass_scales_focusing_only(declared_experiment):
    config, (_scales, baseline) = declared_experiment
    changed = copy.deepcopy(config)
    changed["remnant_mass_msun"] *= 2
    ratio = ((changed["remnant_mass_msun"] + changed["target_binary_mass_msun"])
             / (config["remnant_mass_msun"] + config["target_binary_mass_msun"]))
    _, rows = experiment_rows(changed)
    for actual, original in zip(rows, baseline):
        assert actual["tau_geometric"] == pytest.approx(original["tau_geometric"], rel=2e-12)
        assert actual["tau_focusing"] == pytest.approx(ratio*original["tau_focusing"], rel=2e-12)


def test_optical_depth_against_physical_rate_units(declared_experiment):
    config, (scales, rows) = declared_experiment
    constants = config["illustrative_constants"]
    radius_pc = config["passage_radius_au"]*constants["au_in_km"]/constants["pc_in_km"]
    area_pc2 = math.pi*radius_pc**2
    focusing_speed2 = (2*constants["G_pc_kms2_per_msun"]
                       * (config["remnant_mass_msun"]+config["target_binary_mass_msun"])
                       / radius_pc)
    velocity_unit = scales["velocity_unit_kms"]
    seconds_per_time_unit = scales["dynamical_time_myr"]*constants["seconds_per_myr"]
    for row in rows:
        # n (pc^-3) * area (pc^2) * speed (km/s) * time (s) / (km/pc).
        integrated_speed = (velocity_unit*row["geometric_integral"]
                            + focusing_speed2/velocity_unit*row["focusing_integral"])
        expected = (scales["central_binary_density_pc_minus3"] * area_pc2
                    * integrated_speed * seconds_per_time_unit / constants["pc_in_km"])
        assert row["tau_total"] == pytest.approx(expected, rel=2e-12)
    stationary = next(row for row in rows if row["eta"] == 0)
    assert stationary["encounter_ratio"] == pytest.approx(1., rel=2e-12)
    assert stationary["proxy_relative_error"] == pytest.approx(0., abs=2e-12)


def test_near_escape_case_is_rejected_by_declared_apocenter_filter(declared_experiment):
    config, _ = declared_experiment
    changed = copy.deepcopy(config)
    changed["eta"] = [.99]
    with pytest.raises(ValueError, match="apocenter"):
        experiment_rows(changed)
