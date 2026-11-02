"""Descriptive uncertainty for independent matched-transport ensembles.

Normal and Wilson intervals are approximations. The union-DKW envelope is a
finite-sample confidence bound for the population CDF sup distance under iid
sampling within each ensemble. None of these summaries decides numerical
convergence or physical diffusion adequacy. Samples include stopped paths with
the full declared horizon in their denominator; numerical stops are not escape.
"""
import math

import numpy as np


DEFAULT_STOP_CAUSES = ('horizon', 'excluded_channel', 'spatial_domain',
                       'velocity_domain')
DEFAULT_NUMERICAL_STOP_CAUSES = ('spatial_domain', 'velocity_domain')


def _confidence(value):
    value = float(value)
    if not math.isfinite(value) or not 0 < value < 1:
        raise ValueError('confidence must be finite and strictly between 0 and 1')
    return value


def _normal_critical(confidence):
    """Central normal quantile by bisection, avoiding local module shadowing.

    This file is named statistics.py, so importing stdlib statistics from a
    script in this directory could import this module instead. erfc retains
    tail accuracy for confidence levels close to one; erf handles small levels.
    """
    confidence = _confidence(confidence)
    low, high = 0., 40.
    for _ in range(100):
        mid = (low + high) / 2
        below = (math.erf(mid / math.sqrt(2)) < confidence if confidence <= .5
                 else math.erfc(mid / math.sqrt(2)) > 1 - confidence)
        if below:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def _samples(values, name):
    samples = np.asarray(values, dtype=float)
    if samples.ndim != 1 or samples.size < 2 or not np.isfinite(samples).all():
        raise ValueError(f'{name} must be a finite one-dimensional sample of at least two paths')
    return samples


def mean_summary(values, confidence=.95):
    """Unbiased sample variance and an approximate iid normal mean interval."""
    confidence = _confidence(confidence)
    values = _samples(values, 'values')
    mean = float(values.mean())
    variance = float(values.var(ddof=1))
    se = math.sqrt(variance / len(values))
    half = _normal_critical(confidence) * se
    return {'n': int(len(values)), 'mean': mean, 'sample_variance': variance,
            'standard_error': se, 'confidence': confidence,
            'interval': [mean - half, mean + half],
            'minimum': float(values.min()), 'maximum': float(values.max()),
            'interval_method': 'normal approximation using independent paths',
            'zero_variance_caveat': 'A constant observed sample does not establish zero population variance'}


def wilson_probability(count, n, confidence=.95):
    """Two-sided Wilson score interval (approximate coverage, no continuity correction)."""
    confidence = _confidence(confidence)
    if (isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer))
            or n < 1 or isinstance(count, (bool, np.bool_))
            or not isinstance(count, (int, np.integer)) or not 0 <= count <= n):
        raise ValueError('count and n must be integers with n > 0 and 0 <= count <= n')
    count, n = int(count), int(n)
    p = count / n
    z = _normal_critical(confidence)
    denominator = 1 + z*z/n
    center = (p + z*z/(2*n)) / denominator
    half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denominator
    return {'count': count, 'n': n, 'estimate': p, 'confidence': confidence,
            'interval': [max(0., center-half), min(1., center+half)],
            'interval_method': 'Wilson score, approximate coverage'}


def empirical_cdf_distance(first, second, confidence=.95):
    """Empirical two-CDF sup and a union-DKW population-distance envelope.

    Each empirical CDF receives failure probability alpha/2. On their joint
    coverage event, |D_empirical-D_population| <= epsilon_first+epsilon_second.
    This is not a p-value, and the confidence envelope need not be tight.
    """
    confidence = _confidence(confidence)
    first, second = _samples(first, 'first'), _samples(second, 'second')
    support = np.unique(np.concatenate((first, second)))
    cdf_first = np.searchsorted(np.sort(first), support, side='right') / len(first)
    cdf_second = np.searchsorted(np.sort(second), support, side='right') / len(second)
    distance = float(np.max(np.abs(cdf_first-cdf_second)))
    log_bound = math.log(4 / (1-confidence))
    eps_first = math.sqrt(log_bound / (2*len(first)))
    eps_second = math.sqrt(log_bound / (2*len(second)))
    radius = eps_first+eps_second
    return {'estimate': distance, 'n_first': int(len(first)),
            'n_second': int(len(second)), 'confidence': confidence,
            'epsilon_first': eps_first, 'epsilon_second': eps_second,
            'radius': radius,
            'interval': [max(0., distance-radius), min(1., distance+radius)],
            'interval_method': 'union-DKW bound for population CDF sup distance',
            'assumptions': 'iid paths within each ensemble; all stopped paths retained'}


def compare_samples(first, second, confidence=.95):
    """Independent first-minus-second mean contrast and CDF discrepancy.

    Signed and absolute mean intervals use normal approximation; the absolute
    interval is the image of the signed interval under absolute value. Neither
    interval treats paths with equal index as paired observations.
    """
    confidence = _confidence(confidence)
    first, second = _samples(first, 'first'), _samples(second, 'second')
    first_summary = mean_summary(first, confidence)
    second_summary = mean_summary(second, confidence)
    difference = first_summary['mean']-second_summary['mean']
    se = math.sqrt(first_summary['standard_error']**2 +
                   second_summary['standard_error']**2)
    half = _normal_critical(confidence)*se
    low, high = difference-half, difference+half
    absolute = [0. if low <= 0 <= high else min(abs(low), abs(high)),
                max(abs(low), abs(high))]
    return {'first': first_summary, 'second': second_summary,
            'mean_difference': {'estimate': difference, 'standard_error': se,
                'interval': [low, high], 'absolute_estimate': abs(difference),
                'absolute_interval': absolute, 'confidence': confidence,
                'direction': 'first minus second',
                'interval_method': 'normal approximation for independent ensembles'},
            'empirical_cdf_distance': empirical_cdf_distance(first, second, confidence)}


def summarize_samples(residence, exposure, stop_time, stop_cause, events,
                      horizon, confidence=.95, *,
                      numerical_stop_causes=DEFAULT_NUMERICAL_STOP_CAUSES,
                      known_stop_causes=DEFAULT_STOP_CAUSES):
    """Summarize full-horizon normalized rewards and stopping diagnostics.

    ``residence`` and ``exposure`` are already normalized by the full horizon;
    stop_time and horizon use the same time unit. Exposure need not be bounded
    by one for a generic normalized density field. The supplied events are
    resolved finite-event counts, not the scattering count of a diffusion arm.

    For numerical-domain stops, b=(H-t_stop)/H bounds each omitted future
    residence fraction. Its empirical mean is reported separately from the
    probability-based upper bound E[b] <= P(numerical stop). The Wilson upper
    endpoint is an approximate-confidence bound on that probability, not an
    exact finite-sample confidence guarantee or a physical escape fraction.
    """
    confidence = _confidence(confidence)
    horizon = float(horizon)
    if not math.isfinite(horizon) or horizon <= 0:
        raise ValueError('horizon must be positive and finite')
    residence = _samples(residence, 'residence')
    exposure = _samples(exposure, 'exposure')
    stop_time = _samples(stop_time, 'stop_time')
    events = _samples(events, 'events')
    causes = np.asarray(stop_cause)
    n = len(residence)
    if (len(exposure) != n or len(stop_time) != n or len(events) != n
            or causes.ndim != 1 or len(causes) != n):
        raise ValueError('all path arrays must have the same one-dimensional shape')
    if (any(not isinstance(c, str) or not c for c in causes.tolist())
            or any(not isinstance(c, str) or not c for c in known_stop_causes)
            or any(not isinstance(c, str) or not c for c in numerical_stop_causes)):
        raise ValueError('stop causes must be nonempty strings')
    known = tuple(known_stop_causes)
    numerical = tuple(numerical_stop_causes)
    if len(set(known)) != len(known) or len(set(numerical)) != len(numerical):
        raise ValueError('registered stop causes must be unique')
    if not set(causes.tolist()).issubset(known) or not set(numerical).issubset(known):
        raise ValueError('all observed and numerical stop causes must be registered')
    tolerance = 1e-10
    normalized_time = stop_time / horizon
    if (np.any(normalized_time < -tolerance) or np.any(normalized_time > 1+tolerance)
            or np.any(abs(normalized_time[causes == 'horizon']-1) > tolerance)):
        raise ValueError('stop times must lie within the horizon and horizon stops must finish it')
    if (np.any(residence < -tolerance) or np.any(residence > 1+tolerance)
            or np.any(residence > normalized_time+tolerance)):
        raise ValueError('residence fractions must lie between zero and elapsed time divided by horizon')
    if np.any(exposure < 0):
        raise ValueError('exposure must be nonnegative')
    if np.any(events < 0) or np.any(events != np.floor(events)):
        raise ValueError('events must be nonnegative integer counts')
    probability = {cause: wilson_probability(int(np.sum(causes == cause)), n, confidence)
                   for cause in known}
    numerical_mask = np.isin(causes, numerical)
    # Times within roundoff of either endpoint were accepted above; clipping
    # preserves the exact pathwise bound b <= indicator despite that tolerance.
    missing = np.clip(1-normalized_time, 0., 1.)*numerical_mask
    missing_summary = mean_summary(missing, confidence)
    numerical_probability = wilson_probability(int(numerical_mask.sum()), n, confidence)
    return {'n': n, 'horizon': horizon, 'confidence': confidence,
            'residence': mean_summary(residence, confidence),
            'exposure': mean_summary(exposure, confidence),
            'stop_probabilities': probability,
            'numerical_domain_missing_residence': {
                'remaining_time_bound': missing_summary,
                'stop_probability': numerical_probability,
                'probability_based_mean_upper_bound': numerical_probability['interval'][1],
                'numerical_stop_causes': list(numerical),
                'interpretation': 'Omitted residence <= remaining-time bound <= numerical-stop indicator pathwise',
                'upper_bound_method': 'Wilson probability upper endpoint; approximate coverage'},
            'event_counts': {**mean_summary(events, confidence),
                'total': sum(int(value) for value in events),
                'interpretation': 'Supplied resolved finite-event counts; zero in diffusion does not mean no scattering'},
            'inference_scope': 'Descriptive marginal intervals, not simultaneous across observables or comparisons'}
