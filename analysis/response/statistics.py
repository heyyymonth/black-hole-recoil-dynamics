"""Sampling uncertainty for the discrete residence-response experiment.

Each response sample must be an unbiased complete-transition insertion score,
either N times a normalized reward difference or the probability-weighted
four-branch contrast. Independent rows may share random numbers within their
own continuations. Validation is a separate batch, independent of that response
batch. Reusing one validation baseline or one response estimate across mixture
parameters creates dependence between reported errors. Bonferroni coverage
does not require independence of those errors, but its normal component
intervals remain approximations. No interval here bounds discretization error,
numerical-domain censoring, physical uncertainty, or extrapolation bias.
"""
import math

import numpy as np

from analysis.matched.statistics import _confidence, _normal_critical, _samples


def _family(family_size, confidence):
    confidence = _confidence(confidence)
    if (isinstance(family_size, (bool, np.bool_))
            or not isinstance(family_size, (int, np.integer)) or family_size < 1):
        raise ValueError('family_size must be a positive integer fixed before evaluation')
    family_size = int(family_size)
    component = 1-(1-confidence)/family_size
    if component == 1:
        raise ValueError('family_size is too large for floating-point confidence resolution')
    return confidence, family_size, component


def _moments(values, name):
    values = _samples(values, name)
    mean, variance = float(values.mean()), float(values.var(ddof=1))
    if not math.isfinite(mean) or not math.isfinite(variance):
        raise ValueError(f'{name} moments overflow; use correctly normalized finite samples')
    return values, mean, variance, variance/len(values)


def _interval(estimate, variance_of_mean, confidence, family_size):
    confidence, family_size, component = _family(family_size, confidence)
    se = math.sqrt(variance_of_mean)
    half = _normal_critical(confidence)*se
    adjusted_half = _normal_critical(component)*se
    low, high = estimate-adjusted_half, estimate+adjusted_half
    absolute_low = 0. if low <= 0 <= high else min(abs(low), abs(high))
    return {'estimate': float(estimate), 'standard_error': se,
            'variance_of_estimate': float(variance_of_mean),
            'marginal_confidence': confidence,
            'marginal_interval': [estimate-half, estimate+half],
            'family_confidence': confidence, 'family_size': family_size,
            'component_confidence': component,
            'interval': [low, high], 'absolute_estimate': abs(estimate),
            'absolute_interval': [absolute_low, max(abs(low), abs(high))],
            'interval_method': 'Normal approximation with Bonferroni component level',
            'coverage_scope': 'Approximate family coverage only if every selected inferential comparison is counted in the predeclared family'}


def _normalized_rewards(values, name):
    values = _samples(values, name)
    if np.any(values < -1e-10) or np.any(values > 1+1e-10):
        raise ValueError(f'{name} must be full-horizon normalized residence in [0,1]')
    return values


def replacement_response_samples(diffusion_values, replacement_values, n_steps):
    """N*(Q_replacement-Q_diffusion) for paired Q=C/T, without survivor conditioning.

    Uniformly randomized insertion among all N complete transitions, including
    their integrated rewards and absorbing outcomes, is required by the
    estimator identity. Other insertion designs need their own weights.
    """
    if (isinstance(n_steps, (bool, np.bool_))
            or not isinstance(n_steps, (int, np.integer)) or n_steps < 1):
        raise ValueError('n_steps must be a positive integer')
    diffusion_values = _normalized_rewards(diffusion_values, 'diffusion_values')
    replacement_values = _normalized_rewards(replacement_values, 'replacement_values')
    if diffusion_values.shape != replacement_values.shape:
        raise ValueError('replacement and diffusion arrays must contain aligned pairs')
    return int(n_steps)*(replacement_values-diffusion_values)


def summarize_response(response_samples, confidence=.95, family_size=1):
    """Mean slope and sampling intervals from independent insertion scores.

    A large zero-contribution fraction is retained explicitly. Estimated normal
    uncertainty can be unreliable when a few rare insertions dominate; zero
    observed variance is not evidence of an identically zero response.
    """
    values, mean, variance, mean_variance = _moments(response_samples, 'response_samples')
    return {**_interval(mean, mean_variance, confidence, family_size),
            'n_pairs': int(len(values)), 'sample_variance': variance,
            'minimum': float(values.min()), 'maximum': float(values.max()),
            'zero_contribution_count': int(np.sum(values == 0)),
            'zero_contribution_fraction': float(np.mean(values == 0)),
            'normalization': 'Complete-transition response score for full-horizon normalized reward C/T; estimator weights supplied by caller',
            'interpretation': 'Derivative at a=0 of the discrete whole-transition mixture; not an endpoint identity',
            'rare_contribution_caveat': 'Normal intervals use observed variance; rare unseen contributions and zero observed variance require separate diagnostics'}


def prediction_error(values_a, baseline_values, response_samples, a,
                     confidence=.95, family_size=1, paired_validation=False):
    """Estimate [mu(a)-mu(0)]-a*S with response/validation batch independence.

    Default: mixture and validation-baseline arms are independent, and may have
    different sample counts. With paired_validation=True their rows must be
    aligned independent pairs with correct marginal laws; uncertainty uses the
    variance of paired differences. The same baseline may be reused across a,
    and the same S may enter every prediction. The Bonferroni family size must
    count all inferential residual comparisons across a, launches and numerical
    settings selected for the simultaneous statement. It is not inferred from
    how many times this function is called.

    The response batch must be independent of BOTH validation arrays. Sharing a
    response baseline with validation requires covariance terms not supplied
    by this API. Zero is excluded as a because it is the reference, not an
    independent validation contrast with itself.
    """
    a = float(a)
    if not math.isfinite(a) or not 0 < a <= 1:
        raise ValueError('a must be finite and in (0,1]')
    if not isinstance(paired_validation, (bool, np.bool_)):
        raise ValueError('paired_validation must be an explicit boolean')
    values_a = _normalized_rewards(values_a, 'values_a')
    baseline_values = _normalized_rewards(baseline_values, 'baseline_values')
    _, mean_a, variance_a, mean_variance_a = _moments(values_a, 'values_a')
    _, mean_0, variance_0, mean_variance_0 = _moments(baseline_values, 'baseline_values')
    slope_values, slope, slope_variance, slope_mean_variance = _moments(response_samples, 'response_samples')
    shift = mean_a-mean_0
    if paired_validation:
        if values_a.shape != baseline_values.shape:
            raise ValueError('paired validation arrays must have equal shape')
        _, shift, shift_variance, shift_mean_variance = _moments(
            values_a-baseline_values, 'paired_validation_difference')
        covariance = float(np.cov(values_a, baseline_values, ddof=1)[0, 1])
    else:
        shift_mean_variance = mean_variance_a+mean_variance_0
        shift_variance, covariance = None, None
    predicted_shift = a*slope
    predicted_variance = a*a*slope_mean_variance
    error = shift-predicted_shift
    result = _interval(error, shift_mean_variance+predicted_variance,
                       confidence, family_size)
    return {**result, 'a': a,
            'direction': '[mixture mean minus diffusion mean] minus a times response slope',
            'observed_shift': float(shift),
            'predicted_shift': float(predicted_shift),
            'observed_shift_standard_error': math.sqrt(shift_mean_variance),
            'predicted_shift_standard_error': math.sqrt(predicted_variance),
            'variance_components': {'validation_shift': float(shift_mean_variance),
                                    'response_prediction': float(predicted_variance)},
            'validation': {'n_mixture': int(len(values_a)),
                'n_baseline': int(len(baseline_values)),
                'mixture_mean': mean_a, 'baseline_mean': mean_0,
                'mixture_sample_variance': variance_a,
                'baseline_sample_variance': variance_0,
                'paired': bool(paired_validation),
                'paired_difference_sample_variance': shift_variance,
                'paired_sample_covariance': covariance},
            'response': {'n_pairs': int(len(slope_values)), 'slope': slope,
                         'sample_variance': slope_variance,
                         'variance_of_slope_estimate': slope_mean_variance},
            'assumptions': 'Response batch independent of validation batch; independent rows within each batch; shared validation baseline across a allowed',
            'numerical_scope': 'Sampling error for the fixed numerical process; no time/grid/censoring allowance or physical-accuracy verdict'}
