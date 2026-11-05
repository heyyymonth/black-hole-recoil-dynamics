"""Sampling and efficiency estimates for the fixed precision pilot.

Rows are independent quartets. Correlation between branches within a quartet
is already retained in its response score. Time-stratified designs have four
fixed equal-probability strata; their sampling variance must not be estimated
by treating the pooled scores as identically distributed.

Bootstrap intervals are exploratory percentile intervals, with a nominal
Bonferroni component level. They do not provide exact simultaneous coverage,
account for unseen rare contributions, or include timing uncertainty.
"""
from analysis.matched.statistics import _normal_critical

import numpy as np


STRATUM_WEIGHTS = np.full(4, .25)


def _groups(scores, strata):
    values = np.asarray(scores, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.isfinite(values).all():
        raise ValueError('scores must be a finite one-dimensional array with at least two rows')
    if strata is None:
        return values, [(None, 1., values)]
    labels = np.asarray(strata)
    if (labels.shape != values.shape or labels.dtype.kind not in 'iu'
            or set(np.unique(labels).tolist()) != {0, 1, 2, 3}):
        raise ValueError('strata must contain integer labels 0, 1, 2, 3, each represented')
    groups = [(int(label), float(weight), values[labels == label])
              for label, weight in enumerate(STRATUM_WEIGHTS)]
    if any(len(group) < 2 for _, _, group in groups):
        raise ValueError('each fixed stratum requires at least two rows')
    return values, groups


def estimate(scores, strata=None):
    """Estimate the response mean and its design-correct sampling variance.

    ``strata=None`` denotes independent uniform-time insertion. Otherwise
    integer labels 0..3 identify equal-probability time intervals, irrespective
    of the number of rows assigned to each. All rows, including stopped-prefix
    zero scores, remain in their assigned stratum.
    """
    values, groups = _groups(scores, strata)
    components = []
    for label, weight, group in groups:
        components.append({'stratum': label, 'weight': weight, 'n': len(group),
            'mean': float(group.mean()), 'sample_variance': float(group.var(ddof=1))})
    mean = sum(item['weight']*item['mean'] for item in components)
    variance = sum(item['weight']**2*item['sample_variance']/item['n']
                   for item in components)
    absolute = np.abs(values)
    total_absolute = float(absolute.sum())
    largest_count = max(1, int(np.ceil(.01*len(values))))
    top_share = (float(np.partition(absolute, len(values)-largest_count)[-largest_count:].sum())
                 / total_absolute if total_absolute else 0.)
    return {'n': len(values), 'mean': float(mean), 'variance_of_mean': float(variance),
        'standard_error': float(np.sqrt(variance)),
        'equivalent_score_variance': float(len(values)*variance),
        'minimum': float(values.min()), 'maximum': float(values.max()),
        'zero_count': int(np.count_nonzero(values == 0)),
        'absolute_score_quantiles': {str(q): float(np.quantile(absolute, q))
                                     for q in (.5, .9, .99, .999)},
        'largest_one_percent_absolute_score_share': top_share,
        'design': 'uniform' if strata is None else 'four_equal_probability_time_strata',
        'strata': components,
        'scope': 'Sampling variance for independent rows under the declared fixed design; no numerical or physical error allowance'}


def _confidence(confidence, family_size):
    if not np.isfinite(confidence) or not 0 < confidence < 1:
        raise ValueError('confidence must lie strictly between zero and one')
    if type(family_size) is not int or family_size < 1:
        raise ValueError('family_size must be a positive integer')
    alpha = (1-confidence)/family_size
    return alpha, _normal_critical(1-alpha)


def compare(first, second, *, confidence=.95, family_size=1):
    """Difference of two independent arm estimates, not an equivalence test."""
    alpha, critical = _confidence(confidence, family_size)
    means = np.array([first['mean'], second['mean']], dtype=float)
    variances = np.array([first['variance_of_mean'], second['variance_of_mean']], dtype=float)
    if not np.isfinite(means).all() or not np.isfinite(variances).all() or np.any(variances < 0):
        raise ValueError('means must be finite and variances finite and nonnegative')
    difference = float(means[0]-means[1])
    variance = float(variances.sum())
    se = float(np.sqrt(variance))
    return {'difference': difference, 'variance_of_difference': variance,
        'standard_error': se, 'interval': [difference-critical*se, difference+critical*se],
        'confidence': confidence, 'family_size': family_size,
        'component_confidence': 1-alpha,
        'interval_method': 'Normal approximation with Bonferroni component level',
        'scope': 'Independent arm difference; an interval containing zero does not establish equivalence'}


def _bootstrap_variances(groups, reps, rng, batch_size=64):
    """Resample each design stratum independently, with bounded working memory."""
    output = np.zeros(reps)
    for _, weight, values in groups:
        count = len(values)
        for start in range(0, reps, batch_size):
            size = min(batch_size, reps-start)
            sample = values[rng.integers(0, count, size=(size, count))]
            output[start:start+size] += weight**2*sample.var(axis=1, ddof=1)/count
    return output


def efficiency_bootstrap(scores_a, strata_a, scores_b, strata_b, *, reps=8192,
                         seed=2026110501, family_size=6, confidence=.95,
                         cost_a=None, cost_b=None):
    """Exploratory variance-efficiency ratio A/B for equal sample counts.

    A ratio above one favors B (lower variance). With both costs supplied,
    ``cost_adjusted`` is Var(mean_A)*cost_A / (Var(mean_B)*cost_B).
    Costs are total CPU time for all retained rows in that arm and are held
    fixed within bootstrap resamples: no timing confidence interval is claimed.
    A zero denominator in any resample makes the interval unavailable, with
    the count retained. No invalid replicates are silently discarded.
    """
    a, groups_a = _groups(scores_a, strata_a)
    b, groups_b = _groups(scores_b, strata_b)
    if len(a) != len(b):
        raise ValueError('efficiency comparisons require equal total sample counts')
    if type(reps) is not int or reps < 2:
        raise ValueError('reps must be an integer of at least two')
    alpha, _ = _confidence(confidence, family_size)
    if (cost_a is None) != (cost_b is None):
        raise ValueError('supply both total CPU costs or neither')
    if cost_a is not None and (not np.isfinite([cost_a, cost_b]).all()
                              or min(cost_a, cost_b) <= 0):
        raise ValueError('CPU costs must be finite and positive')
    ea, eb = estimate(a, strata_a), estimate(b, strata_b)
    result = {'n_per_arm': len(a), 'bootstrap_replicates': reps,
        'confidence': confidence, 'family_size': family_size,
        'component_confidence': 1-alpha,
        'interval_method': 'Independent within-design percentile bootstrap with nominal Bonferroni component level',
        'ratio_direction': 'A variance divided by B variance; values above one favor B',
        'scope': 'Exploratory approximate intervals; unseen rare scores and timing uncertainty are not covered',
        'variance_ratio': {'estimate': None, 'interval': None},
        'cost_adjusted': None,
        'invalid_bootstrap_denominators': None,
        'unavailable_reason': None}
    if eb['variance_of_mean'] == 0:
        result['unavailable_reason'] = 'Observed denominator variance is zero'
        return result
    ratio = ea['variance_of_mean']/eb['variance_of_mean']
    result['variance_ratio']['estimate'] = float(ratio)
    child_a, child_b = np.random.SeedSequence(seed).spawn(2)
    va = _bootstrap_variances(groups_a, reps, np.random.default_rng(child_a))
    vb = _bootstrap_variances(groups_b, reps, np.random.default_rng(child_b))
    invalid = int(np.count_nonzero(vb == 0))
    result['invalid_bootstrap_denominators'] = invalid
    if invalid:
        result['unavailable_reason'] = 'At least one bootstrap denominator variance is zero; no replicates discarded'
    else:
        interval = np.quantile(va/vb, [alpha/2, 1-alpha/2])
        result['variance_ratio']['interval'] = interval.tolist()
    if cost_a is not None:
        factor = float(cost_a/cost_b)
        interval = result['variance_ratio']['interval']
        result['cost_adjusted'] = {'estimate': float(ratio*factor),
            'interval': None if interval is None else [float(x*factor) for x in interval],
            'total_cpu_seconds_a': float(cost_a), 'total_cpu_seconds_b': float(cost_b),
            'timing_scope': 'Conditional on recorded total CPU costs; timing variability excluded'}
    return result
