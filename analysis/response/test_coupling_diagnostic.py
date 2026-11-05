"""Independent exact-noise-moment checks of the local coupling diagnostic."""
import json

import numpy as np
import pytest

from analysis.response.coupling_diagnostic import (
    OUT, check_factor, contrast, factors, verify,
)


def test_covariance_equivalence_does_not_determine_common_noise_coupling():
    D = np.diag([1., 2., 3.])
    first = factors(D)
    second = factors(D)
    second['eigenvectors'] = (-np.asarray(second['eigenvectors'])).tolist()
    second['eigenfactor'] = (-np.asarray(second['eigenfactor'])).tolist()
    # A simultaneous eigenvector sign reversal leaves the unique PSD root fixed.
    check_factor(first); check_factor(second)
    result = contrast(first, second)
    assert result['D_frobenius_difference'] == 0.
    assert result['symmetric_root_common_noise_mse_per_unit_dt'] == 0.
    np.testing.assert_allclose(result['eigenfactor_common_noise_mse_per_unit_dt'], 4*np.trace(D))


def test_common_noise_mse_matches_exact_second_moment_cubature():
    first = factors(np.array([[2., .3, .1], [.3, 1., .2], [.1, .2, 3.]]))
    second = factors(np.array([[2.1, .2, .1], [.2, .9, .25], [.1, .25, 3.2]]))
    # These six nodes exactly integrate every quadratic under E[ZZ^T]=I.
    nodes = np.r_[np.sqrt(3)*np.eye(3), -np.sqrt(3)*np.eye(3)]
    result = contrast(first, second)
    dt = .037
    for label in ('eigenfactor', 'symmetric_root'):
        difference = np.asarray(second[label])-np.asarray(first[label])
        increments = np.sqrt(dt)*nodes@difference.T
        measured_mse = np.mean(np.sum(increments**2, axis=1))
        np.testing.assert_allclose(measured_mse,
            dt*result[label+'_common_noise_mse_per_unit_dt'], rtol=1e-13)


def test_factor_verification_detects_wrong_covariance():
    record = factors(np.diag([1., 2., 3.]))
    record['eigenfactor'][0][0] += .1
    with pytest.raises(AssertionError):
        check_factor(record)


def test_retained_deterministic_diagnostic_verifies_without_reconstruction():
    assert verify(json.loads(OUT.read_text()))
