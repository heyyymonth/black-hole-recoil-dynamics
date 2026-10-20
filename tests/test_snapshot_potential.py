"""Independent checks of the frozen-shell potential and return integrals."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from recoil_dynamics.plummer import PlummerOrbit
from recoil_dynamics.snapshot_potential import ShellPotential


def single_shell_return(radius, mass, kick):
    """Closed form: central free flight plus a radial Kepler arc."""
    binding = mass / radius - kick * kick / 2
    turning = mass / binding
    z = radius / turning
    outside = mass / (math.sqrt(2) * binding ** 1.5) * (
        math.pi / 2 - math.asin(math.sqrt(z)) + math.sqrt(z * (1 - z)))
    return 2 * (radius / kick + outside)


def plummer_shells(count):
    # Midpoint quantiles of the analytic cumulative mass r^3/(1+r^2)^1.5.
    quantile = (np.arange(count) + 0.5) / count
    q23 = quantile ** (2 / 3)
    return ShellPotential(np.sqrt(q23 / (1 - q23)), np.full(count, 1 / count))


def test_potential_direct_sum_and_exterior_kepler():
    r = np.array([0.15, 0.4, 1.2, 3.0])
    m = np.array([0.1, 0.2, 0.4, 0.3])
    model = ShellPotential(r, m)
    queries = np.array([0, 0.01, 0.15, 0.3, 0.4, 1, 1.2, 2, 3, 10])
    direct = np.array([-sum(mi / max(ri, x) for ri, mi in zip(r, m)) for x in queries])
    np.testing.assert_allclose(model.potential(queries), direct, rtol=2e-15, atol=2e-15)
    assert model.central_potential == pytest.approx(direct[0])
    assert model.potential(10) == pytest.approx(-model.total_mass / 10)
    assert model.acceleration(10) == pytest.approx(-model.total_mass / 100)
    assert model.acceleration(0) == 0
    assert model.acceleration(0.1) == 0
    assert model.acceleration(r[0]) == 0


def test_acceleration_from_independent_potential_finite_difference():
    model = ShellPotential([0.2, 0.7, 1.1, 2.8], [0.15, 0.25, 0.3, 0.3])
    for radius in (0.05, 0.4, 0.9, 1.7, 4.0):
        step = 1e-5 * max(radius, 0.1)
        numerical = -(model.potential(radius + step) - model.potential(radius - step)) / (2 * step)
        assert model.acceleration(radius) == pytest.approx(numerical, rel=2e-9, abs=2e-10)


def test_sorting_duplicate_combination_and_shell_phase_independence():
    # No angular phase information enters a spherical monopole. Reordering and
    # splitting mass at a fixed radius must leave every result unchanged.
    original = ShellPotential([0.2, 0.8, 2], [0.2, 0.5, 0.3])
    reordered = ShellPotential([2, 0.8, 0.2, 0.8], [0.3, 0.15, 0.2, 0.35])
    np.testing.assert_allclose(original.potential(np.linspace(0, 4, 101)),
                               reordered.potential(np.linspace(0, 4, 101)), atol=1e-15)
    for eta in (0.1, 0.5, 0.95):
        kick = eta * original.escape_speed
        assert original.apocenter(kick) == pytest.approx(reordered.apocenter(kick), rel=2e-14)
        assert original.first_return(kick) == pytest.approx(reordered.first_return(kick), rel=2e-13)


@pytest.mark.parametrize("radius,mass", [(1.0, 1.0), (0.4, 2.3)])
@pytest.mark.parametrize("eta", [0.01, 0.3, 0.8, 0.99, 0.9999])
def test_single_shell_analytic_return(radius, mass, eta):
    model = ShellPotential([radius], [mass])
    kick = eta * model.escape_speed
    expected = single_shell_return(radius, mass, kick)
    assert model.apocenter(kick) == pytest.approx(radius / (1 - eta * eta), rel=3e-12)
    assert model.first_return(kick) == pytest.approx(expected, rel=3e-10)


def test_multishell_return_independent_radius_quadrature():
    radii = np.array([0.1, 0.3, 0.6, 1.4, 2.5])
    masses = np.array([0.1, 0.15, 0.2, 0.25, 0.3])
    model = ShellPotential(radii, masses)
    for eta in (0.2, 0.55, 0.85, 0.98):
        kick = eta * model.escape_speed
        turning = model.apocenter(kick)
        energy = model.central_potential + kick * kick / 2
        def independent_integrand(radius):
            phi = -sum(mass / max(shell, radius) for shell, mass in zip(radii, masses))
            return 1 / math.sqrt(2 * (energy - phi))
        boundaries = np.concatenate(([0], radii[radii < turning], [turning]))
        reference = 2 * sum(quad(independent_integrand, lo, hi, epsabs=1e-9,
                                 epsrel=1e-10, limit=200)[0]
                            for lo, hi in zip(boundaries[:-1], boundaries[1:]))
        assert model.first_return(kick) == pytest.approx(reference, rel=3e-8)
        assert model.potential(turning) == pytest.approx(energy, rel=5e-14, abs=5e-14)


def test_high_resolution_plummer_recovery_and_discretization_convergence():
    coarse, fine = plummer_shells(2048), plummer_shells(32768)
    assert fine.central_potential == pytest.approx(-1, rel=5e-4)
    for eta in (0.3, 0.6, 0.9):
        # Same physical kick in units G=M=b=1, not the shell escape fraction.
        kick = math.sqrt(2) * eta
        reference = PlummerOrbit(eta).period
        coarse_error = abs(coarse.first_return(kick) / reference - 1)
        fine_error = abs(fine.first_return(kick) / reference - 1)
        assert fine_error < coarse_error
        assert fine_error < 2e-3


def test_time_and_length_scaling():
    original = ShellPotential([0.1, 0.8, 2], [0.2, 0.5, 0.3])
    length, mass = 7.0, 2.5
    scaled = ShellPotential(original.radii * length, original.masses * mass)
    for eta in (0.2, 0.8):
        kick = eta * original.escape_speed
        scaled_kick = kick * math.sqrt(mass / length)
        assert scaled.apocenter(scaled_kick) == pytest.approx(length * original.apocenter(kick))
        assert scaled.first_return(scaled_kick) == pytest.approx(
            math.sqrt(length ** 3 / mass) * original.first_return(kick), rel=2e-12)


def test_innermost_force_free_region_and_stationary_exception():
    model = ShellPotential([0.1, 0.8, 2], [0.2, 0.5, 0.3])
    assert model.apocenter(0) == 0
    with pytest.raises(ValueError, match="stationary"):
        model.first_return(0)
    for eta in (1e-4, 1e-6):
        kick = eta * model.escape_speed
        assert model.first_return(kick) / (2 * model.radii[0] / kick) == pytest.approx(1, rel=2e-7)


@pytest.mark.parametrize("radii,masses", [([], []), ([0], [1]), ([-1], [1]),
    ([1], [0]), ([1], [-1]), ([1, 2], [1]), ([[1]], [[1]]),
    ([math.inf], [1]), ([1], [math.nan])])
def test_invalid_shells(radii, masses):
    with pytest.raises(ValueError):
        ShellPotential(radii, masses)


def test_invalid_evaluation_kicks_and_tolerances():
    model = ShellPotential([1], [1])
    for radius in (-1, math.nan, math.inf):
        with pytest.raises(ValueError):
            model.potential(radius)
        with pytest.raises(ValueError):
            model.acceleration(radius)
    for kick in (-1, math.nan, math.inf, model.escape_speed, 2 * model.escape_speed):
        with pytest.raises(ValueError):
            model.apocenter(kick)
        with pytest.raises(ValueError):
            model.first_return(kick)
    for tolerance in (-1, 0, math.nan, math.inf):
        with pytest.raises(ValueError):
            model.first_return(0.5, epsrel=tolerance)
        with pytest.raises(ValueError):
            model.first_return(0.5, epsabs=tolerance)
