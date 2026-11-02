"""Matched frozen-law jump/diffusion transport in dimensionless coordinates.

Each step freezes the collision law at its initial state. An independent
exponential excluded-channel clock shortens the step to delta. Two conservative
velocity-Verlet half-flows surround either a Gaussian increment (A delta,
D delta) or the sum of Poisson(Gamma delta) finite increments. Thus the frozen
velocity updates have exactly matching mean and covariance; the evolving
process still requires time-step and kernel-grid refinement.

Gravity drifts are straight segments of the declared splitting. Their spatial
exits and sphere occupation are calculated geometrically. Velocity-domain exits
are detected at resolved kick/flow endpoints, not by a continuous stochastic
bridge. All stops retain the full requested horizon as the normalization.
"""
import numpy as np


def sphere_residence_fraction(start, end, radius=1.):
    """Fraction of a straight segment strictly inside a sphere, vectorized.

    Tangencies have zero measure. A constant segment on the boundary is outside;
    a chord with both endpoints on the boundary has its interior counted.
    """
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    if (start.ndim != 2 or start.shape[1] != 3 or end.shape != start.shape
            or not np.isfinite(start).all() or not np.isfinite(end).all()
            or not np.isfinite(radius) or radius <= 0):
        raise ValueError('finite N by 3 endpoints and a positive radius required')
    difference = end-start
    a = np.einsum('ij,ij->i',difference,difference)
    b = np.einsum('ij,ij->i',start,difference)
    c = np.einsum('ij,ij->i',start,start)-radius*radius
    result = np.where(a == 0, (c < 0).astype(float), 0.)
    discriminant = b*b-a*c
    selected = (a > 0) & (discriminant > 0)
    if np.any(selected):
        aa, bb, cc = a[selected],b[selected],c[selected]
        root = np.sqrt(discriminant[selected])
        q = -bb-np.where(bb >= 0,root,-root)
        first = q/aa
        second = np.divide(cc,q,out=np.zeros_like(q),where=q != 0)
        lower,upper = np.minimum(first,second),np.maximum(first,second)
        result[selected] = np.maximum(np.minimum(upper,1.)-np.maximum(lower,0.),0.)
    return result


def _exit_fraction(start, end, radius):
    """First outer-sphere exit for a segment whose initial point is inside."""
    delta = end-start
    a = np.einsum('ij,ij->i',delta,delta)
    b = np.einsum('ij,ij->i',start,delta)
    c = np.einsum('ij,ij->i',start,start)-radius*radius
    crossed = np.einsum('ij,ij->i',end,end) >= radius*radius
    fraction = np.ones(len(start))
    if np.any(crossed):
        aa,bb,cc = a[crossed],b[crossed],c[crossed]
        discriminant = np.maximum(bb*bb-aa*cc,0.)
        root = np.sqrt(discriminant)
        # Stable positive root of a f^2+2b f+c=0 with c<0.
        value = np.empty_like(root)
        outward = bb >= 0
        value[outward] = -cc[outward]/(bb[outward]+root[outward])
        value[~outward] = (-bb[~outward]+root[~outward])/aa[~outward]
        fraction[crossed] = np.clip(value,0.,1.)
    return fraction,crossed


def _checked_moments(kernel, states):
    supplied = kernel.moments(states.copy())
    n = len(states)
    expected = {'rate':(n,), 'excluded':(n,), 'A':(n,3), 'D':(n,3,3)}
    moments = {}
    for key,shape in expected.items():
        value = np.asarray(supplied[key],dtype=float)
        if value.shape != shape or not np.isfinite(value).all():
            raise ValueError(f'kernel {key} must be finite and shaped {shape}')
        moments[key] = value
    if np.any(moments['rate'] < 0) or np.any(moments['excluded'] < 0):
        raise ValueError('encounter and excluded rates must be nonnegative')
    D = moments['D']
    scale = np.maximum(np.max(abs(D),axis=(1,2)),np.finfo(float).tiny)
    if np.any(np.max(abs(D-D.transpose(0,2,1)),axis=(1,2)) > 1e-12*scale):
        raise ValueError('diffusion tensor is not symmetric')
    moments['D'] = (D+D.transpose(0,2,1))/2
    return moments


def _gaussian_increment(A, D, elapsed, normal):
    eigenvalues,eigenvectors = np.linalg.eigh(D)
    scale = np.maximum(np.max(abs(eigenvalues),axis=1),np.finfo(float).tiny)
    if np.any(eigenvalues[:,0] < -1e-12*scale):
        raise ValueError('diffusion tensor is not positive semidefinite')
    roots = np.sqrt(np.maximum(eigenvalues,0.)*elapsed[:,None])
    return A*elapsed[:,None]+np.einsum('nij,nj->ni',eigenvectors,roots*normal)


def run_ensemble(kernel, acceleration, initial_states, horizon, dt, seed, method,
                 spatial_stop=4., velocity_stop=4., density=None, region_radius=1.,
                 event_batch_size=65536):
    """Return stopped paths and occupation for independent trajectory samples.

    States are (x,y,z,vx,vy,vz) in consistent dimensionless units. ``kernel``
    supplies ``moments(Nx6)`` with rate, excluded, A and raw D per unit time,
    and ``sample_events(Nx6,rng)`` returning one finite increment per row.
    ``acceleration`` maps Nx3 positions to Nx3 vectors. Optional ``density``
    maps Nx3 positions to a nonnegative exposure weight, sampled at the midpoint
    of each modeled gravity drift. Its absence is explicitly recorded.

    Stops: spatial_domain at a drift's first exit; velocity_domain at a resolved
    gravity kick or collision update; excluded_channel at its sampled elapsed
    time. At equal times spatial/velocity validity stops precede excluded stops.
    Initial points at either domain boundary stop immediately. A horizon result
    retains alive=True. No stopping label asserts irreversible physical escape.

    Gaussian and excluded-clock arrays have fixed shape on every step, including
    dead paths. Jump counts and finite-event sampling use separate streams.
    Different arms/seeds are independent ensembles: no paired error reduction
    or per-path coupling across runs is asserted by this function.
    """
    states = np.asarray(initial_states,dtype=float).copy()
    if (states.ndim != 2 or states.shape[1] != 6 or len(states) == 0
            or not np.isfinite(states).all()):
        raise ValueError('initial_states must be finite nonempty N by 6')
    controls = [horizon,dt,spatial_stop,velocity_stop,region_radius]
    if not np.isfinite(controls).all() or min(controls) <= 0:
        raise ValueError('horizon, dt and radii/speed limits must be positive and finite')
    if method not in ('diffusion','jump'):
        raise ValueError('method must be diffusion or jump')
    if (not isinstance(event_batch_size,(int,np.integer)) or event_batch_size < 1):
        raise ValueError('event_batch_size must be a positive integer')
    n = len(states)
    streams = np.random.SeedSequence(seed).spawn(4)
    normal_rng,kill_rng,count_rng,event_rng = [np.random.default_rng(s) for s in streams]
    alive = np.ones(n,dtype=bool)
    cause = np.full(n,'horizon',dtype='U24')
    end_time = np.full(n,float(horizon))
    residence = np.zeros(n)
    exposure = np.zeros(n)
    event_count = np.zeros(n,dtype=np.int64)
    evaluated_steps = np.zeros(n,dtype=np.int64)
    n_steps = int(np.ceil(horizon/dt))

    def stop(indices, time, label):
        alive[indices] = False
        cause[indices] = label
        end_time[indices] = time

    outside = np.linalg.norm(states[:,:3],axis=1) >= spatial_stop
    stop(np.flatnonzero(outside),0.,'spatial_domain')
    outside_v = alive & (np.linalg.norm(states[:,3:],axis=1) >= velocity_stop)
    stop(np.flatnonzero(outside_v),0.,'velocity_domain')

    def force(position):
        value = np.asarray(acceleration(position.copy()),dtype=float)
        if value.shape != position.shape or not np.isfinite(value).all():
            raise ValueError('acceleration must return finite N by 3 vectors')
        return value

    def flow(indices, elapsed, starting_time):
        """One VV gravity flow, retaining its stopped piecewise-linear drift."""
        if not len(indices):
            return
        initial_x = states[indices,:3].copy()
        half_v = states[indices,3:]+.5*elapsed[:,None]*force(initial_x)
        states[indices,3:] = half_v
        velocity_exit = np.linalg.norm(half_v,axis=1) >= velocity_stop
        if np.any(velocity_exit):
            stop(indices[velocity_exit],starting_time[velocity_exit],'velocity_domain')
        keep = ~velocity_exit
        indices,elapsed,starting_time = indices[keep],elapsed[keep],starting_time[keep]
        initial_x,half_v = initial_x[keep],half_v[keep]
        if not len(indices):
            return
        proposed_x = initial_x+elapsed[:,None]*half_v
        fraction,spatial_exit = _exit_fraction(initial_x,proposed_x,spatial_stop)
        final_x = initial_x+fraction[:,None]*(proposed_x-initial_x)
        duration = elapsed*fraction
        residence[indices] += duration*sphere_residence_fraction(initial_x,final_x,region_radius)
        if density is not None:
            weight = np.asarray(density((initial_x+final_x)/2),dtype=float)
            if weight.shape != (len(indices),) or not np.isfinite(weight).all() or np.any(weight < 0):
                raise ValueError('density must return a finite nonnegative weight per position')
            exposure[indices] += duration*weight
        states[indices,:3] = final_x
        if np.any(spatial_exit):
            stop(indices[spatial_exit],starting_time[spatial_exit]+duration[spatial_exit],
                 'spatial_domain')
        keep = ~spatial_exit
        indices,elapsed,starting_time = indices[keep],elapsed[keep],starting_time[keep]
        if len(indices):
            states[indices,3:] += .5*elapsed[:,None]*force(states[indices,:3])
            velocity_exit = np.linalg.norm(states[indices,3:],axis=1) >= velocity_stop
            if np.any(velocity_exit):
                stop(indices[velocity_exit],starting_time[velocity_exit]+elapsed[velocity_exit],
                     'velocity_domain')

    for step in range(n_steps):
        if not np.any(alive):
            break
        time = step*dt
        interval = min(dt,horizon-time)
        if interval <= 0:
            break
        normal = normal_rng.standard_normal((n,3))
        exponential = kill_rng.exponential(size=n)
        ids = np.flatnonzero(alive)
        initial = states[ids].copy()
        moments = _checked_moments(kernel,initial)
        evaluated_steps[ids] += 1
        killing_time = np.divide(exponential[ids],moments['excluded'],
                                 out=np.full(len(ids),np.inf),where=moments['excluded'] > 0)
        elapsed = np.minimum(interval,killing_time)
        excluded = killing_time <= interval
        flow(ids,elapsed/2,np.full(len(ids),time))
        keep = alive[ids]
        collision_ids = ids[keep]
        if len(collision_ids):
            if method == 'diffusion':
                increments = _gaussian_increment(moments['A'][keep],moments['D'][keep],
                                                  elapsed[keep],normal[collision_ids])
            else:
                counts = count_rng.poisson(moments['rate'][keep]*elapsed[keep])
                increments = np.zeros((len(collision_ids),3))
                repeated = np.repeat(np.arange(len(collision_ids)),counts)
                frozen_states = initial[keep]
                for start in range(0,len(repeated),event_batch_size):
                    rows = repeated[start:start+event_batch_size]
                    change = np.asarray(kernel.sample_events(frozen_states[rows],event_rng),dtype=float)
                    if change.shape != (len(rows),3) or not np.isfinite(change).all():
                        raise ValueError('sample_events must return a finite N by 3 increment')
                    np.add.at(increments,rows,change)
                event_count[collision_ids] += counts
            states[collision_ids,3:] += increments
            velocity_exit = np.linalg.norm(states[collision_ids,3:],axis=1) >= velocity_stop
            if np.any(velocity_exit):
                stop(collision_ids[velocity_exit],time+elapsed[keep][velocity_exit]/2,
                     'velocity_domain')
        keep = alive[ids]
        flow(ids[keep],elapsed[keep]/2,time+elapsed[keep]/2)
        killed = alive[ids] & excluded
        if np.any(killed):
            stop(ids[killed],time+elapsed[killed],'excluded_channel')

    if (np.any(residence < -1e-12) or np.any(residence > end_time+1e-10)
            or not np.isfinite(states).all() or not np.isfinite(exposure).all()):
        raise ArithmeticError('invalid stopped trajectory accounting')
    return {'states':states,'residence':residence,'normalized_residence':residence/horizon,
            'exposure':exposure,'normalized_exposure':exposure/horizon,
            'exposure_evaluated':density is not None,'end_time':end_time,
            'stop_cause':cause,'alive':alive,'elastic_events':event_count,
            'evaluated_steps':evaluated_steps,'n_paths':n,'n_steps':n_steps,
            'horizon':float(horizon),'dt':float(dt),'method':method,
            'scheme':'initial-state frozen collision law, two VV half-flows, midpoint velocity update',
            'interpretation':'independent conditional ensembles; time/grid and stopping refinement required'}
