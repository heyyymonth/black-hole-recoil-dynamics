"""Decision and decomposition checks, independent of physical trajectories."""
import json

import numpy as np
import pytest

from analysis.response.analyze import PROTOCOL, classification, contribution_bins


@pytest.mark.parametrize('interval,expected', [
    ((-.0009, .0008), 'within_tolerance'),
    ((-.001, .001), 'within_tolerance'),
    ((.0011, .0012), 'outside_tolerance'),
    ((-.002, -.0011), 'outside_tolerance'),
    ((.0009, .0011), 'inconclusive'),
    ((-.0011, -.0009), 'inconclusive'),
    ((-.002, .002), 'inconclusive'),
    # Touching the acceptance boundary does not demonstrate disagreement.
    ((.001, .002), 'inconclusive'),
    ((-.002, -.001), 'inconclusive'),
])
def test_endpoint_decisions_use_entire_interval(interval, expected):
    assert classification(interval, .001) == expected


def test_signed_bins_partition_retained_score_with_full_sample_denominator():
    coordinate = np.array([0., .25, .5, 1., 2., 4., 0.])
    active = np.array([True, True, True, True, True, True, False])
    score = np.array([7., -14., 21., -28., 35., -42., 49.])
    result = contribution_bins({'active_at_replacement': active}, score,
                               [0., .5, 1., 2., 4.], coordinate)
    assert [r['count'] for r in result['bins']] == [2, 1, 1, 2]
    np.testing.assert_allclose([r['contribution'] for r in result['bins']], [-1., 3., -4., -1.])
    assert result['inactive_count'] == 1
    assert result['inactive_contribution'] == 7.
    np.testing.assert_allclose(sum(r['contribution'] for r in result['bins'])
                               + result['inactive_contribution'], score.mean())
    assert result['denominator'] == len(score)


def test_registered_family_covers_all_reported_inferential_contrasts():
    p = json.loads(PROTOCOL.read_text())
    launches = len(p['launches'])
    settings = p['numerical_settings']
    responses = launches*len(settings)
    validation = launches*sum(len(p['validation']['primary_a_values']
        if setting['role'] == 'primary' else p['validation']['control_a_values'])
        for setting in settings)
    # Grid and timestep contrasts for both response and endpoint mean shift.
    refinements = launches*2*2
    family = responses+2*validation+refinements
    assert (responses, validation, refinements, family) == (9, 21, 12, 63)
    assert family <= p['reporting']['bonferroni_family_size']


def test_coordinates_outside_bins_are_explicitly_retained():
    score = np.array([2., -6., 12., -4.])
    result = contribution_bins({'active_at_replacement': np.array([True, True, True, False])},
                               score, [0., 1., 2.], np.array([.5, -1., 3., 0.]))
    assert result['outside_count'] == 2
    assert result['outside_contribution'] == 1.5
    np.testing.assert_allclose(sum(row['contribution'] for row in result['bins'])
        + result['outside_contribution']+result['inactive_contribution'], score.mean())
