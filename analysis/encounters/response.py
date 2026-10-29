"""First response by a single transition replacement.

For N steps with transition kernel P_a=(1-a)P_D+a P_J, differentiating
the N-fold product at a=0 gives a sum of N terms. Each term has one P_J-P_D
insertion between a diffusion prefix and suffix. Uniformly sampling the
insertion step K therefore gives the unbiased estimator

    N * (observable_with_one_jump_step - diffusion_observable).

This is exact for the declared DISCRETE transition-mixture family. It is
not an exact finite-time continuous-generator response: step refinement is
required. A prediction for the all-jump endpoint is a first-order
extrapolation, not an identity. Shared noise is variance reduction only.
"""

import numpy as np


def _positive_integer(value, name):
    if isinstance(value, (bool, np.bool_)) or int(value) != value or value < 1:
        raise ValueError(f'{name} must be a positive integer')
    return int(value)


def finite_state_response(diffusion, jump, initial, reward_diffusion,
                          reward_jump, n_steps, terminal_reward=None):
    """Exact baseline and derivative for finite-state transition matrices.

    ``reward_*[i]`` is the expected additive reward of one corresponding
    transition starting at state i. Transition matrices are row stochastic.
    The terminal reward, if supplied, is common to both descriptions.
    """
    n_steps = _positive_integer(n_steps, 'n_steps')
    pd = np.asarray(diffusion, dtype=float)
    pj = np.asarray(jump, dtype=float)
    if pd.ndim != 2 or pd.shape[0] != pd.shape[1] or pj.shape != pd.shape:
        raise ValueError('transition matrices must have the same square shape')
    n_states = pd.shape[0]
    for transition in [pd, pj]:
        if (not np.isfinite(transition).all() or np.any(transition < 0)
                or not np.allclose(transition.sum(axis=1), 1., rtol=0, atol=1e-12)):
            raise ValueError('transition matrices must be row stochastic')
    p0 = np.asarray(initial, dtype=float)
    rd = np.asarray(reward_diffusion, dtype=float)
    rj = np.asarray(reward_jump, dtype=float)
    terminal = (np.zeros(n_states) if terminal_reward is None
                else np.asarray(terminal_reward, dtype=float))
    for vector in [p0, rd, rj, terminal]:
        if vector.shape != (n_states,) or not np.isfinite(vector).all():
            raise ValueError('initial and reward vectors must be finite and state sized')
    if np.any(p0 < 0) or not np.isclose(p0.sum(), 1., rtol=0, atol=1e-12):
        raise ValueError('initial law must be a probability distribution')
    value = terminal.copy()
    response = np.zeros(n_states)
    for _ in range(n_steps):
        response = rj - rd + (pj - pd) @ value + pd @ response
        value = rd + pd @ value
    return {
        'baseline': float(p0 @ value),
        'derivative': float(p0 @ response),
        'baseline_by_initial_state': value,
        'derivative_by_initial_state': response,
        'interpretation': 'exact derivative of the discrete transition mixture at zero',
    }


def randomized_replacement_response(diffusion_step, jump_step, initial_states,
                                    n_steps, dt, seed, replacement_steps=None,
                                    initial_time=0., terminal_reward=None):
    """Estimate the first response using paired diffusion continuations.

    A callback has signature ``step(states, time, dt, rng)`` and returns
    ``(new_states, integrated_reward)``. The input ensemble is N-by-d, but
    each callback receives ONE path as a 1-by-d array and returns a length-1
    reward already integrated over this step.
    Absorption must be represented in state and callbacks (zero subsequent
    reward). A callback must not change its input. Each row represents an
    independent path; callbacks must not couple different rows.

    One uniformly random replacement step is chosen independently for each
    path by default. ``replacement_steps`` permits a caller-supplied schedule
    of length N; its weighting and uncertainty interpretation then belong to
    the caller. Uniform randomized stratification preserves the derivative,
    while an arbitrary selected step estimates a different quantity.

    Independent per-path, per-step RNG streams are shared only within each
    baseline/replacement pair, including its jump insertion. State-dependent
    draw counts therefore cannot couple different paths. Identical callbacks
    cancel path by path. Callbacks must use only their supplied RNG and must
    not retain cross-call state. The replacement schedule has a separate RNG.

    ``normalized_response`` divides the raw derivative by the full horizon
    N*dt, including paths killed early; it does not renormalize survivors.
    Reported standard errors assume independent sampled paths. They exclude
    time discretization and physical model uncertainty.
    """
    n_steps = _positive_integer(n_steps, 'n_steps')
    if not np.isfinite([dt, initial_time]).all() or dt <= 0:
        raise ValueError('dt must be positive and initial_time finite')
    states = np.asarray(initial_states, dtype=float)
    if states.ndim != 2 or len(states) < 2 or not np.isfinite(states).all():
        raise ValueError('initial_states must be a finite N-by-d array with N>=2')
    n_paths = len(states)
    schedule_seed, path_seed = np.random.SeedSequence(seed).spawn(2)
    schedule_rng = np.random.default_rng(schedule_seed)
    path_rng = np.random.default_rng(path_seed)
    if replacement_steps is None:
        insertion = schedule_rng.integers(0, n_steps, size=n_paths)
        schedule_description = 'independent uniform step for each path'
    else:
        insertion = np.asarray(replacement_steps)
        if (insertion.shape != (n_paths,) or not np.isfinite(insertion).all()
                or np.any(insertion != insertion.astype(int))
                or np.any(insertion < 0) or np.any(insertion >= n_steps)):
            raise ValueError('replacement_steps must contain one valid integer per path')
        insertion = insertion.astype(int)
        schedule_description = 'caller supplied; caller must establish weighting and uncertainty'

    baseline = states.copy()
    replaced = states.copy()
    observable_d = np.zeros(n_paths)
    observable_replaced = np.zeros(n_paths)

    def evaluate(callback, incoming, time, callback_seeds):
        outputs = np.empty_like(incoming)
        rewards = np.empty(len(incoming))
        for row, callback_seed in enumerate(callback_seeds):
            output, reward = callback(incoming[row:row+1].copy(), time, dt,
                                     np.random.default_rng(callback_seed))
            output, reward = np.asarray(output, dtype=float), np.asarray(reward, dtype=float)
            if (output.shape != (1, incoming.shape[1]) or reward.shape != (1,)
                    or not np.isfinite(output).all() or not np.isfinite(reward).all()):
                raise ValueError('step must return finite same-shaped states and one reward per path')
            outputs[row], rewards[row] = output[0], reward[0]
        return outputs, rewards

    for step_index in range(n_steps):
        time = initial_time + step_index * dt
        common_seeds = path_rng.integers(0, 2**63, size=n_paths)
        next_d, reward_d = evaluate(diffusion_step, baseline, time, common_seeds)
        next_replaced, reward_replaced = evaluate(
            diffusion_step, replaced, time, common_seeds)
        selected = insertion == step_index
        if np.any(selected):
            jump_states, jump_reward = evaluate(
                jump_step, replaced[selected], time, common_seeds[selected])
            next_replaced[selected] = jump_states
            reward_replaced[selected] = jump_reward
        # The same state and same path-specific RNG share the prefix exactly.
        future = insertion > step_index
        next_replaced[future] = next_d[future]
        reward_replaced[future] = reward_d[future]
        baseline, replaced = next_d, next_replaced
        observable_d += reward_d
        observable_replaced += reward_replaced

    if terminal_reward is not None:
        for state, observable in [(baseline, observable_d),
                                  (replaced, observable_replaced)]:
            for row in range(n_paths):
                reward = np.asarray(terminal_reward(
                    state[row:row+1].copy(), initial_time + n_steps * dt), dtype=float)
                if reward.shape != (1,) or not np.isfinite(reward).all():
                    raise ValueError('terminal reward must return one finite value per path')
                observable[row] += reward[0]
    raw = n_steps * (observable_replaced - observable_d)
    normalized = raw / (n_steps * dt)
    return {
        'raw_response': raw,
        'normalized_response': normalized,
        'raw_mean': float(np.mean(raw)),
        'raw_standard_error': float(np.std(raw, ddof=1) / np.sqrt(n_paths)),
        'normalized_mean': float(np.mean(normalized)),
        'normalized_standard_error': float(np.std(normalized, ddof=1) / np.sqrt(n_paths)),
        'baseline_observable': observable_d,
        'replaced_observable': observable_replaced,
        'replacement_steps': insertion,
        'n_paths': n_paths,
        'n_steps': n_steps,
        'dt': float(dt),
        'horizon': float(n_steps * dt),
        'schedule': schedule_description,
        'randomness': 'independent per-path per-step streams shared only within a pair',
        'interpretation': 'first response of a discrete transition mixture; endpoint extrapolation untested',
    }
