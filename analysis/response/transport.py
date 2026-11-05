"""Complete-step mixtures and paired single-replacement transport.

The transition below is the same frozen-law numerical transition as the
published matched pilot. It includes both gravity half-flows, accumulated
residence/exposure, and every stopping rule. A mixture chooses between complete
diffusion and jump transitions, rather than mixing their moments or replacing
only a final velocity. This module does not modify the hash-bound pilot solver.

Within one response pair, fixed-shape normal and unit-exponential arrays are
shared row by row, even after one branch stops. Finite-event counts and samples
have separate RNG streams. Pairing is asserted only within that constructed
pair, never between independent runs or timesteps.
"""
import numpy as np

from analysis.matched.transport import (
    _checked_moments, _exit_fraction, _gaussian_increment,
    sphere_residence_fraction,
)


def _validate_controls(horizon,dt,spatial_stop,velocity_stop,region_radius,event_batch_size):
    values=[horizon,dt,spatial_stop,velocity_stop,region_radius]
    if not np.isfinite(values).all() or min(values)<=0:
        raise ValueError('horizon, dt and radii/speed limits must be positive and finite')
    if not isinstance(event_batch_size,(int,np.integer)) or event_batch_size<1:
        raise ValueError('event_batch_size must be a positive integer')


def initialize_batch(initial_states,horizon,spatial_stop=4.,velocity_stop=4.):
    """Create a stopped Markov state including all additive observables."""
    states=np.asarray(initial_states,dtype=float).copy()
    if (states.ndim!=2 or states.shape[1]!=6 or not len(states)
            or not np.isfinite(states).all()):
        raise ValueError('initial_states must be finite nonempty N by 6')
    if not np.isfinite([horizon,spatial_stop,velocity_stop]).all() or min(horizon,spatial_stop,velocity_stop)<=0:
        raise ValueError('horizon and domain limits must be positive and finite')
    n=len(states)
    batch={'states':states,'alive':np.ones(n,dtype=bool),
           'stop_cause':np.full(n,'horizon',dtype='U24'),
           'end_time':np.full(n,float(horizon)),
           'residence':np.zeros(n),'exposure':np.zeros(n),
           'elastic_events':np.zeros(n,dtype=np.int64),
           'evaluated_steps':np.zeros(n,dtype=np.int64)}
    outside=np.linalg.norm(states[:,:3],axis=1)>=spatial_stop
    _stop(batch,np.flatnonzero(outside),0.,'spatial_domain')
    outside_v=batch['alive']&(np.linalg.norm(states[:,3:],axis=1)>=velocity_stop)
    _stop(batch,np.flatnonzero(outside_v),0.,'velocity_domain')
    return batch


def _stop(batch,indices,time,label):
    batch['alive'][indices]=False
    batch['stop_cause'][indices]=label
    batch['end_time'][indices]=time


def advance_step(kernel,acceleration,batch,time,interval,normal,exponential,
                 jump_mask,count_rng,event_rng,spatial_stop=4.,velocity_stop=4.,
                 density=None,region_radius=1.,event_batch_size=65536,forced_counts=None):
    """Mutate one complete stopped transition using caller-supplied randomness.

    ``normal`` is N by 3, ``exponential`` is N unit-exponential draws, and
    ``jump_mask`` is N booleans. Dead rows are absorbing. The caller draws all
    rows before either branch advances so state-dependent deaths cannot shift
    the paired Gaussian or excluded-channel streams. Each active row freezes
    its own coefficients at its own initial state.

    Spatial exits use exact straight-drift crossing and occupation. The
    excluded clock shortens the full step. Resolved velocity/spatial validity
    stops take precedence over excluded stopping; no post-stop reward accrues.
    Density exposure uses midpoint quadrature separately on each clipped drift.
    ``forced_counts``, when supplied, is an N-row nonnegative integer array
    used only on jump rows that reach the collision update. It supports explicit
    probability-weighted conditional-count branches, not an unconditional law.
    """
    states=batch['states'];alive=batch['alive'];n=len(states)
    normal=np.asarray(normal,dtype=float)
    exponential=np.asarray(exponential,dtype=float)
    jump_mask=np.asarray(jump_mask)
    if (normal.shape!=(n,3) or exponential.shape!=(n,) or jump_mask.shape!=(n,)
            or jump_mask.dtype.kind!='b' or not np.isfinite(normal).all()
            or not np.isfinite(exponential).all() or np.any(exponential<0)):
        raise ValueError('finite aligned N-row normal/exponential arrays and boolean jump_mask required')
    if not np.isfinite([time,interval]).all() or time<0 or interval<=0:
        raise ValueError('step time must be nonnegative and interval positive')
    if forced_counts is not None:
        forced_counts=np.asarray(forced_counts)
        if forced_counts.shape!=(n,) or forced_counts.dtype.kind not in 'iu' or np.any(forced_counts<0):
            raise ValueError('forced_counts must contain N nonnegative integers')

    def force(position):
        value=np.asarray(acceleration(position.copy()),dtype=float)
        if value.shape!=position.shape or not np.isfinite(value).all():
            raise ValueError('acceleration must return finite N by 3 vectors')
        return value

    def flow(indices,elapsed,starting_time):
        if not len(indices):
            return
        initial_x=states[indices,:3].copy()
        half_v=states[indices,3:]+.5*elapsed[:,None]*force(initial_x)
        states[indices,3:]=half_v
        velocity_exit=np.linalg.norm(half_v,axis=1)>=velocity_stop
        if np.any(velocity_exit):
            _stop(batch,indices[velocity_exit],starting_time[velocity_exit],'velocity_domain')
        keep=~velocity_exit
        indices,elapsed,starting_time=indices[keep],elapsed[keep],starting_time[keep]
        initial_x,half_v=initial_x[keep],half_v[keep]
        if not len(indices):
            return
        proposed_x=initial_x+elapsed[:,None]*half_v
        fraction,spatial_exit=_exit_fraction(initial_x,proposed_x,spatial_stop)
        final_x=initial_x+fraction[:,None]*(proposed_x-initial_x)
        duration=elapsed*fraction
        batch['residence'][indices]+=duration*sphere_residence_fraction(initial_x,final_x,region_radius)
        if density is not None:
            weight=np.asarray(density((initial_x+final_x)/2),dtype=float)
            if weight.shape!=(len(indices),) or not np.isfinite(weight).all() or np.any(weight<0):
                raise ValueError('density must return a finite nonnegative weight per position')
            batch['exposure'][indices]+=duration*weight
        states[indices,:3]=final_x
        if np.any(spatial_exit):
            _stop(batch,indices[spatial_exit],starting_time[spatial_exit]+duration[spatial_exit],
                  'spatial_domain')
        keep=~spatial_exit
        indices,elapsed,starting_time=indices[keep],elapsed[keep],starting_time[keep]
        if len(indices):
            states[indices,3:]+=.5*elapsed[:,None]*force(states[indices,:3])
            velocity_exit=np.linalg.norm(states[indices,3:],axis=1)>=velocity_stop
            if np.any(velocity_exit):
                _stop(batch,indices[velocity_exit],starting_time[velocity_exit]+elapsed[velocity_exit],
                      'velocity_domain')

    ids=np.flatnonzero(alive)
    if not len(ids):
        return
    initial=states[ids].copy()
    moments=_checked_moments(kernel,initial)
    batch['evaluated_steps'][ids]+=1
    killing_time=np.divide(exponential[ids],moments['excluded'],
                           out=np.full(len(ids),np.inf),where=moments['excluded']>0)
    elapsed=np.minimum(interval,killing_time)
    excluded=killing_time<=interval
    flow(ids,elapsed/2,np.full(len(ids),time))
    keep=alive[ids]
    collision_ids=ids[keep]
    if len(collision_ids):
        increments=np.zeros((len(collision_ids),3))
        jump=jump_mask[collision_ids]
        diffusion=~jump
        if np.any(diffusion):
            increments[diffusion]=_gaussian_increment(moments['A'][keep][diffusion],
                moments['D'][keep][diffusion],elapsed[keep][diffusion],normal[collision_ids[diffusion]])
        if np.any(jump):
            counts=(count_rng.poisson(moments['rate'][keep][jump]*elapsed[keep][jump])
                    if forced_counts is None else forced_counts[collision_ids[jump]])
            repeated=np.repeat(np.flatnonzero(jump),counts)
            frozen_states=initial[keep]
            for start in range(0,len(repeated),event_batch_size):
                rows=repeated[start:start+event_batch_size]
                change=np.asarray(kernel.sample_events(frozen_states[rows],event_rng),dtype=float)
                if change.shape!=(len(rows),3) or not np.isfinite(change).all():
                    raise ValueError('sample_events must return a finite N by 3 increment')
                np.add.at(increments,rows,change)
            batch['elastic_events'][collision_ids[jump]]+=counts
        states[collision_ids,3:]+=increments
        velocity_exit=np.linalg.norm(states[collision_ids,3:],axis=1)>=velocity_stop
        if np.any(velocity_exit):
            _stop(batch,collision_ids[velocity_exit],time+elapsed[keep][velocity_exit]/2,'velocity_domain')
    keep=alive[ids]
    flow(ids[keep],elapsed[keep]/2,time+elapsed[keep]/2)
    killed=alive[ids]&excluded
    if np.any(killed):
        _stop(batch,ids[killed],time+elapsed[killed],'excluded_channel')


def _finalize(batch,horizon,dt,method,density):
    residence=batch['residence'];exposure=batch['exposure']
    if (np.any(residence < -1e-12) or np.any(residence>batch['end_time']+1e-10)
            or not np.isfinite(batch['states']).all() or not np.isfinite(exposure).all()):
        raise ArithmeticError('invalid stopped trajectory accounting')
    return {**batch,'normalized_residence':residence/horizon,
            'normalized_exposure':exposure/horizon,'exposure_evaluated':density is not None,
            'n_paths':len(residence),'n_steps':int(np.ceil(horizon/dt)),
            'horizon':float(horizon),'dt':float(dt),'method':method,
            'scheme':'initial-state frozen collision law, two VV half-flows, midpoint velocity update'}


def _streams(seed):
    return [np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(5)]


def run_mixture(kernel,acceleration,initial_states,horizon,dt,seed,mixture_probability,
                spatial_stop=4.,velocity_stop=4.,density=None,region_radius=1.,
                event_batch_size=65536):
    """Run the per-step Markov mixture (1-a) P_diffusion + a P_jump.

    Bernoulli choices are independent across rows/steps and independent of every
    transition draw. The endpoints reproduce the original solver's random
    streams and transition law. Intermediate a does not scale A, D, or Gamma.
    """
    _validate_controls(horizon,dt,spatial_stop,velocity_stop,region_radius,event_batch_size)
    if not np.isscalar(mixture_probability) or not np.isfinite(mixture_probability) or not 0<=mixture_probability<=1:
        raise ValueError('mixture_probability must be a finite scalar in [0,1]')
    batch=initialize_batch(initial_states,horizon,spatial_stop,velocity_stop)
    n=len(batch['states']);steps=int(np.ceil(horizon/dt))
    normal_rng,kill_rng,count_rng,event_rng,choice_rng=_streams(seed)
    for step in range(steps):
        if not np.any(batch['alive']):
            break
        time=step*dt;interval=min(dt,horizon-time)
        if interval<=0:
            break
        normal=normal_rng.standard_normal((n,3));exponential=kill_rng.exponential(size=n)
        jump=choice_rng.random(n)<mixture_probability
        advance_step(kernel,acceleration,batch,time,interval,normal,exponential,jump,count_rng,event_rng,
                     spatial_stop,velocity_stop,density,region_radius,event_batch_size)
    result=_finalize(batch,horizon,dt,'mixture',density)
    result['mixture_probability']=float(mixture_probability)
    return result


def run_paired_response(kernel,acceleration,initial_states,horizon,dt,seed,
                        replacement_steps=None,spatial_stop=4.,velocity_stop=4.,
                        density=None,region_radius=1.,event_batch_size=65536):
    """Pair diffusion with diffusion having exactly one complete jump step.

    If unspecified, each row chooses K uniformly from *all scheduled steps*,
    including steps after an eventual stop. A dead row is absorbing and its
    replacement contribution is zero; such rows must not be discarded. For
    M scheduled steps, M times the paired observable difference is an unbiased
    estimator of the derivative at a=0 of this finite-step mixture's expectation.
    This is a finite-step identity, not a claim about its continuum limit.

    Explicit replacement_steps support deterministic controls or independently
    registered designs. They do not by themselves justify a uniform-K estimator.
    """
    _validate_controls(horizon,dt,spatial_stop,velocity_stop,region_radius,event_batch_size)
    baseline=initialize_batch(initial_states,horizon,spatial_stop,velocity_stop)
    replacement=initialize_batch(initial_states,horizon,spatial_stop,velocity_stop)
    n=len(baseline['states']);steps=int(np.ceil(horizon/dt))
    normal_rng,kill_rng,count_rng,event_rng,choice_rng=_streams(seed)
    randomized=replacement_steps is None
    if randomized:
        selected=choice_rng.integers(0,steps,size=n)
    else:
        selected=np.asarray(replacement_steps)
        if (selected.shape!=(n,) or selected.dtype.kind not in 'iu'
                or np.any(selected<0) or np.any(selected>=steps)):
            raise ValueError('replacement_steps must contain one valid integer step per row')
        selected=selected.copy()
    active_at_replacement=np.zeros(n,dtype=bool)
    replacement_states=np.zeros((n,6))
    no_jumps=np.zeros(n,dtype=bool)
    for step in range(steps):
        if not np.any(baseline['alive']|replacement['alive']):
            break
        time=step*dt;interval=min(dt,horizon-time)
        if interval<=0:
            break
        normal=normal_rng.standard_normal((n,3));exponential=kill_rng.exponential(size=n)
        jump=selected==step
        active_at_replacement[jump]=replacement['alive'][jump]
        active=jump&replacement['alive']
        replacement_states[active]=replacement['states'][active]
        advance_step(kernel,acceleration,baseline,time,interval,normal,exponential,no_jumps,
                     count_rng,event_rng,spatial_stop,velocity_stop,density,region_radius,event_batch_size)
        advance_step(kernel,acceleration,replacement,time,interval,normal,exponential,jump,
                     count_rng,event_rng,spatial_stop,velocity_stop,density,region_radius,event_batch_size)
    return {'baseline':_finalize(baseline,horizon,dt,'diffusion',density),
            'replacement':_finalize(replacement,horizon,dt,'single_jump_replacement',density),
            'replacement_steps':selected,'active_at_replacement':active_at_replacement,
            'replacement_states':replacement_states,
            'replacement_selection':'uniform_scheduled_step' if randomized else 'explicit',
            'n_steps':steps,
            'coupling':'Shared fixed-row Gaussian and excluded-clock draws; independent count/event streams'}


def run_paired_mixture(kernel,acceleration,initial_states,horizon,dt,seed,mixture_probability,
                       spatial_stop=4.,velocity_stop=4.,density=None,region_radius=1.,
                       event_batch_size=65536):
    """Pair the complete-step mixture with diffusion using aligned row noise.

    The Bernoulli choices are independent of transition randomness. This is a
    validation ensemble and must use a seed independent of response estimation
    if a subsequent comparison assumes independent response/validation errors.
    """
    _validate_controls(horizon,dt,spatial_stop,velocity_stop,region_radius,event_batch_size)
    if not np.isscalar(mixture_probability) or not np.isfinite(mixture_probability) or not 0<=mixture_probability<=1:
        raise ValueError('mixture_probability must be a finite scalar in [0,1]')
    baseline=initialize_batch(initial_states,horizon,spatial_stop,velocity_stop)
    mixture=initialize_batch(initial_states,horizon,spatial_stop,velocity_stop)
    n=len(baseline['states']);steps=int(np.ceil(horizon/dt))
    normal_rng,kill_rng,count_rng,event_rng,choice_rng=_streams(seed)
    no_jumps=np.zeros(n,dtype=bool)
    for step in range(steps):
        if not np.any(baseline['alive']|mixture['alive']):
            break
        time=step*dt;interval=min(dt,horizon-time)
        if interval<=0:
            break
        normal=normal_rng.standard_normal((n,3));exponential=kill_rng.exponential(size=n)
        jump=choice_rng.random(n)<mixture_probability
        advance_step(kernel,acceleration,baseline,time,interval,normal,exponential,no_jumps,
                     count_rng,event_rng,spatial_stop,velocity_stop,density,region_radius,event_batch_size)
        advance_step(kernel,acceleration,mixture,time,interval,normal,exponential,jump,
                     count_rng,event_rng,spatial_stop,velocity_stop,density,region_radius,event_batch_size)
    return {'baseline':_finalize(baseline,horizon,dt,'diffusion',density),
            'mixture':_finalize(mixture,horizon,dt,'mixture',density),
            'mixture_probability':float(mixture_probability),
            'coupling':'Shared fixed-row Gaussian and excluded-clock draws; independent Bernoulli/count/event streams'}


def sample_positive_poisson(mu,rng):
    """Sample Poisson(mu) conditioned positive; mu=0 has zero branch weight.

    Small rates use inverse-CDF recurrence starting at mu/expm1(mu), avoiding
    cancellation and rare rejection. Large rates use unconditional Poisson
    draws with zero draws rejected, avoiding underflow of the first mass. Zero
    input returns a count of zero solely as a placeholder for a zero-weight
    branch; it does not define a positive conditional distribution at mu=0.
    """
    mu=np.asarray(mu,dtype=float)
    if mu.ndim!=1 or not np.isfinite(mu).all() or np.any(mu<0):
        raise ValueError('mu must be a finite nonnegative vector')
    result=np.zeros(len(mu),dtype=np.int64)
    small=(mu>0)&(mu<20.)
    if np.any(small):
        means=mu[small]
        target=rng.random(len(means))
        probability=means/np.expm1(means)
        cumulative=probability.copy()
        counts=np.ones(len(means),dtype=np.int64)
        remaining=target>cumulative
        while np.any(remaining):
            counts[remaining]+=1
            probability[remaining]*=means[remaining]/counts[remaining]
            previous=cumulative.copy()
            cumulative[remaining]+=probability[remaining]
            # Remove only machine-roundoff deficits once the tail is below
            # floating-point resolution; this is not a finite-count cap. A
            # cumulative sum may stall more than four eps below one because
            # of accumulated rounding. Once its decreasing tail cannot alter
            # the sum, continuing the recurrence cannot resolve that deficit.
            exhausted=remaining&(counts>means)&(cumulative==previous)
            resolved=((1-cumulative)<=4*np.finfo(float).eps)|exhausted
            cumulative[resolved]=1.
            remaining=target>cumulative
        result[small]=counts
    large=mu>=20.
    if np.any(large):
        means=mu[large]
        counts=rng.poisson(means)
        missing=counts==0
        while np.any(missing):
            counts[missing]=rng.poisson(means[missing])
            missing=counts==0
        result[large]=counts
    return result


def run_weighted_response(kernel,acceleration,initial_states,horizon,dt,seed,
                          replacement_steps=None,spatial_stop=4.,velocity_stop=4.,
                          density=None,region_radius=1.,event_batch_size=65536):
    """Four-branch single-replacement estimator with conditional jump counts.

    At the selected complete transition, decompose its compound-Poisson law
    into zero and positive counts with weights exp(-mu) and -expm1(-mu), where
    mu=Gamma*min(interval,T_excluded) at the common initial state and shared
    excluded clock. Subtract the average of Gaussian +Z and -Z transitions.
    Afterward every branch follows diffusion with shared row-aligned noise.

    Uniform selection over all M scheduled steps makes M times this weighted
    contrast unbiased for the derivative of the finite-step mixture at zero.
    The negative Gaussian branch changes only the selected step's normal draw.
    Stopping, prefix rewards, and within-step rewards remain part of each branch.
    This is intended to reduce variance; it neither changes the derivative nor removes the
    need to test timestep, kernel-grid, and finite-a remainder sensitivity.
    """
    _validate_controls(horizon,dt,spatial_stop,velocity_stop,region_radius,event_batch_size)
    names=('jump_zero','jump_positive','diffusion_plus','diffusion_minus')
    batches={name:initialize_batch(initial_states,horizon,spatial_stop,velocity_stop) for name in names}
    n=len(batches['diffusion_plus']['states']);steps=int(np.ceil(horizon/dt))
    normal_rng,kill_rng,count_rng,event_rng,choice_rng=_streams(seed)
    randomized=replacement_steps is None
    if randomized:
        selected=choice_rng.integers(0,steps,size=n)
    else:
        selected=np.asarray(replacement_steps)
        if (selected.shape!=(n,) or selected.dtype.kind not in 'iu'
                or np.any(selected<0) or np.any(selected>=steps)):
            raise ValueError('replacement_steps must contain one valid integer step per row')
        selected=selected.copy()
    weight_positive=np.zeros(n);weight_zero=np.ones(n);mu_selected=np.zeros(n)
    active_at_replacement=np.zeros(n,dtype=bool)
    replacement_states=np.zeros((n,6))
    positive_counts=np.zeros(n,dtype=np.int64)
    zero_counts=np.zeros(n,dtype=np.int64);no_jumps=np.zeros(n,dtype=bool)
    for step in range(steps):
        if not np.any(np.logical_or.reduce([batch['alive'] for batch in batches.values()])):
            break
        time=step*dt;interval=min(dt,horizon-time)
        if interval<=0:
            break
        normal=normal_rng.standard_normal((n,3));exponential=kill_rng.exponential(size=n)
        jump=selected==step
        active=jump&batches['diffusion_plus']['alive']
        active_at_replacement[active]=True
        replacement_states[active]=batches['diffusion_plus']['states'][active]
        if np.any(active):
            # Before a row's replacement all four states coincide exactly.
            state=batches['diffusion_plus']['states'][active]
            moment=_checked_moments(kernel,state)
            killing_time=np.divide(exponential[active],moment['excluded'],
                out=np.full(len(state),np.inf),where=moment['excluded']>0)
            mu=moment['rate']*np.minimum(interval,killing_time)
            mu_selected[active]=mu
            weight_zero[active]=np.exp(-mu)
            weight_positive[active]=-np.expm1(-mu)
            positive_counts[active]=sample_positive_poisson(mu,count_rng)
        negative_normal=normal.copy();negative_normal[jump]*=-1
        for name in names:
            update_normal=negative_normal if name=='diffusion_minus' else normal
            update_jump=jump if name.startswith('jump_') else no_jumps
            counts=(zero_counts if name=='jump_zero' else
                    positive_counts if name=='jump_positive' else None)
            advance_step(kernel,acceleration,batches[name],time,interval,update_normal,exponential,
                         update_jump,count_rng,event_rng,spatial_stop,velocity_stop,density,
                         region_radius,event_batch_size,forced_counts=counts)
    return {**{name:_finalize(batch,horizon,dt,name,density) for name,batch in batches.items()},
            'jump_zero_weight':weight_zero,'jump_positive_weight':weight_positive,
            'selected_poisson_mean':mu_selected,'positive_branch_count':positive_counts,
            'replacement_steps':selected,'active_at_replacement':active_at_replacement,
            'replacement_states':replacement_states,
            'replacement_selection':'uniform_scheduled_step' if randomized else 'explicit',
            'n_steps':steps,
            'coupling':'Four branches share prefix and continuation Gaussian/excluded draws; Gaussian sign flips only at K; zero/positive compound counts are probability weighted'}
