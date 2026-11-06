"""Sampling intervals and bounded-continuation control for an independent retest.

Normal intervals are approximations using observed variance; they do not
protect against unseen rare scores. The primary continuation control uses an
exact binomial occurrence bound multiplied by a prespecified worst-case score.
The bound requires a fixed number of iid bounded rows. Each row may itself
be a correlated response quartet or a paired validation observation. It cannot certify the continuation model,
physical calibration, or time/grid convergence.
"""
import math

import numpy as np
from analysis.matched.statistics import _normal_critical


def _samples(values, name):
    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError(name+' must contain at least two finite scalar observations')
    return values


def _moments(values, name):
    values = _samples(values, name)
    with np.errstate(over='ignore', invalid='ignore'):
        mean = float(values.mean())
        variance = float(values.var(ddof=1))
    if not math.isfinite(mean) or not math.isfinite(variance):
        raise ValueError(name+' moments overflow; check normalization')
    return values, mean, variance


def _probability(value, name):
    if isinstance(value, (bool, np.bool_)) or not np.isscalar(value):
        raise ValueError(name+' must be a finite scalar in (0,1)')
    value = float(value)
    if not math.isfinite(value) or not 0 < value < 1:
        raise ValueError(name+' must be a finite scalar in (0,1)')
    return value


def _family(confidence, family_size):
    confidence = _probability(confidence, 'confidence')
    if (isinstance(family_size, (bool, np.bool_))
            or not isinstance(family_size, (int, np.integer)) or family_size < 1):
        raise ValueError('family_size must be a positive integer')
    family_size = int(family_size)
    component = 1-(1-confidence)/family_size
    if component == 1.:
        raise ValueError('family confidence is too close to one for floating-point resolution')
    return confidence, family_size, component


def normal_interval(estimate, variance_of_mean, confidence=.975, family_size=64):
    """Two-sided normal interval with a prespecified Bonferroni component level."""
    confidence, family_size, component = _family(confidence, family_size)
    estimate, variance_of_mean = float(estimate), float(variance_of_mean)
    if (not math.isfinite(estimate) or not math.isfinite(variance_of_mean)
            or variance_of_mean < 0):
        raise ValueError('estimate and nonnegative sampling variance must be finite')
    critical = _normal_critical(component)
    se = math.sqrt(variance_of_mean)
    interval = [estimate-critical*se, estimate+critical*se]
    if not np.isfinite(interval).all():
        raise ValueError('normal interval overflows')
    return {'estimate':estimate, 'standard_error':se,
        'variance_of_estimate':variance_of_mean, 'interval':interval,
        'half_width':critical*se, 'normal_critical':critical,
        'family_confidence':confidence, 'family_size':family_size,
        'component_confidence':component,
        'method':'Normal approximation with Bonferroni component level',
        'scope':'Sampling uncertainty for fixed numerical laws; approximate coverage, no protection against unseen rare contributions'}


def mean_interval(samples, confidence=.975, family_size=64):
    """Mean and observed sampling uncertainty for independent identically distributed rows."""
    values, mean, variance = _moments(samples, 'samples')
    return {**normal_interval(mean, variance/len(values), confidence, family_size),
        'n':len(values), 'mean':mean, 'sample_variance':variance,
        'minimum':float(values.min()), 'maximum':float(values.max()),
        'zero_count':int(np.count_nonzero(values == 0))}


def contrast(first, second, confidence=.975, family_size=64):
    """Difference of independent sample means; no paired covariance is assumed."""
    x, mx, vx = _moments(first, 'first')
    y, my, vy = _moments(second, 'second')
    parts = {'first':vx/len(x), 'second':vy/len(y)}
    return {**normal_interval(mx-my, sum(parts.values()), confidence, family_size),
        'n_first':len(x), 'n_second':len(y), 'first_mean':mx, 'second_mean':my,
        'variance_components':parts,
        'assumptions':'Independent batches and independent identically distributed rows within each batch'}


def residual(validation_difference, response_scores, a, confidence=.975, family_size=64):
    """Independent validation shift minus a times the first-response estimate.

    Validation rows are already paired mixture-minus-diffusion differences.
    Their covariance is retained in that difference. The validation batch
    must be independent of the response-score batch. Reusing a prediction
    across several mixture values is allowed only with the full inferential
    family declared; the resulting residuals are not mutually independent.
    """
    if isinstance(a, (bool, np.bool_)) or not np.isscalar(a):
        raise ValueError('a must be a finite scalar in (0,1]')
    a = float(a)
    if not math.isfinite(a) or not 0 < a <= 1:
        raise ValueError('a must be a finite scalar in (0,1]')
    v, mv, vv = _moments(validation_difference, 'validation_difference')
    s, ms, vs = _moments(response_scores, 'response_scores')
    parts = {'validation_shift':vv/len(v), 'response_prediction':a*a*vs/len(s)}
    return {**normal_interval(mv-a*ms, sum(parts.values()), confidence, family_size),
        'a':a, 'observed_shift':mv, 'predicted_shift':a*ms, 'response_mean':ms,
        'n_validation':len(v), 'n_prediction':len(s), 'variance_components':parts,
        'assumptions':'Validation and prediction batches independent; iid rows within each batch; paired validation covariance already included'}


def classify(interval, tolerance):
    """Closed tolerance band: touching from inside passes; crossing is inconclusive."""
    interval = np.asarray(interval, dtype=float)
    if (interval.shape != (2,) or not np.isfinite(interval).all()
            or interval[0] > interval[1]):
        raise ValueError('interval must be two finite ordered endpoints')
    if isinstance(tolerance, (bool, np.bool_)) or not np.isscalar(tolerance):
        raise ValueError('tolerance must be a finite nonnegative scalar')
    tolerance = float(tolerance)
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError('tolerance must be a finite nonnegative scalar')
    low, high = interval
    if low >= -tolerance and high <= tolerance:
        return 'within_tolerance'
    if low > tolerance or high < -tolerance:
        return 'outside_tolerance'
    return 'inconclusive'


def binomial_upper(count, n, *, delta):
    """One-sided Clopper--Pearson upper limit for an iid Bernoulli probability.

    Solve Pr_p[Binomial(n,p) <= count] = delta for count < n. Zero
    occurrences use the stable closed form; all occurrences give one.
    Other counts use a log-sum-exp binomial CDF and return the upper root
    bracket with a small outward floating-point guard. The statistical
    result is exact; numerical evaluation is ordinary float64 arithmetic,
    not a formal interval-arithmetic certificate.
    """
    delta = _probability(delta, 'delta')
    for value, name in ((n, 'n'), (count, 'count')):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise ValueError(name+' must be an integer')
    n, count = int(n), int(count)
    if n < 1 or not 0 <= count <= n:
        raise ValueError('require n >= 1 and 0 <= count <= n')
    if count == n:
        return {'upper_probability':1., 'root_bracket':[1.,1.],
                'numerical_guard':0., 'count':count, 'n':n, 'delta':delta}
    if count == 0:
        value = -math.expm1(math.log(delta)/n)
        upper = float(np.nextafter(value, 1.))
        return {'upper_probability':upper, 'root_bracket':[value,upper],
                'numerical_guard':upper-value, 'count':count, 'n':n, 'delta':delta}
    j = np.arange(count+1, dtype=float)
    log_choose = np.r_[0., np.cumsum(np.log(n-j[1:]+1)-np.log(j[1:]))]
    target = math.log(delta)

    def log_cdf(p):
        terms = log_choose+j*math.log(p)+(n-j)*math.log1p(-p)
        largest = float(terms.max())
        return largest+math.log(float(np.exp(terms-largest).sum()))

    low, high = 0., 1.
    for _ in range(80):
        mid = (low+high)/2
        if mid == low or mid == high:
            break
        if log_cdf(mid) > target:
            low = mid
        else:
            high = mid
    # Outward padding is far below scientific budgets, and protects the final
    # bracket against small CDF/quantile roundoff in the tested n<=65536 range.
    guard = 1e-12
    upper = float(np.nextafter(min(1., high+guard), 1.))
    return {'upper_probability':upper, 'root_bracket':[low,high],
            'numerical_guard':upper-high, 'count':count, 'n':n, 'delta':delta}


def continuation_upper(values, maximum, *, delta):
    """Bound E[B] by M Pr(B>0), then use a one-sided exact binomial limit.

    Each row must be iid, 0<=B<=M, with M and sample size fixed before the
    outcomes. The binomial count is the number of nonzero row-level bounds,
    not a count of correlated branches or encounter events. Supplying a
    separately allocated delta gives a simultaneous population bound by a
    union bound. Zero observed contributions do not imply zero expectation.
    """
    if isinstance(maximum, (bool, np.bool_)) or not np.isscalar(maximum):
        raise ValueError('maximum must be a finite positive deterministic bound')
    maximum = float(maximum)
    if not math.isfinite(maximum) or maximum <= 0:
        raise ValueError('maximum must be a finite positive deterministic bound')
    values = _samples(values, 'values')
    if np.any(values < 0) or np.any(values > maximum):
        raise ValueError('values fall outside the declared deterministic range [0,maximum]')
    count = int(np.count_nonzero(values)); n = len(values)
    limit = binomial_upper(count, n, delta=delta)
    return {'n':n, 'sample_mean_bound':float(values.mean()),
        'nonzero_rows':count, 'maximum':maximum, 'delta':float(delta),
        'confidence':1-float(delta), 'upper_bound':maximum*limit['upper_probability'],
        'occurrence_probability_limit':limit,
        'method':'Clopper-Pearson one-sided occurrence bound times deterministic maximum',
        'reference':'https://doi.org/10.1093/biomet/26.4.404',
        'scope':'Fixed-n iid bounded-score population-mean upper bound; branch dependence retained within rows; no continuum or physical-calibration claim'}


def tail_diagnostics(values):
    """Descriptive influence in observed scores; no population-tail inference."""
    values, mean, _ = _moments(values, 'values')
    absolute = np.abs(values)
    centered = (values-mean)**2
    count = max(1, int(np.ceil(.01*len(values))))
    total_abs = float(absolute.sum()); total_squared = float(centered.sum())
    if not np.isfinite([total_abs, total_squared]).all():
        raise ValueError('tail diagnostics overflow')
    return {'n':len(values), 'largest_one_percent_count':count,
        'largest_one_percent_absolute_score_share':float(np.sort(absolute)[-count:].sum()/total_abs) if total_abs else 0.,
        'largest_single_variance_share':float(centered.max()/total_squared) if total_squared else None,
        'largest_one_percent_variance_share':float(np.sort(centered)[-count:].sum()/total_squared) if total_squared else None,
        'scope':'Observed-sample concentration only; no trimming and no inference about unseen tails'}
