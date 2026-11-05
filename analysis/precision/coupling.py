"""Isolated covariance-factor choices for the precision comparison.

The published response transport is loaded into a private module namespace.
Only that instance's Gaussian-increment function is replaced for the symmetric
arm. The complete transition, all random streams, rewards, and stopping rules
therefore remain those of the frozen source. No imported published module is
mutated. Changing the factor changes path coupling, not the marginal Gaussian
law; individual trajectories need not be equal across the two factors.
"""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from types import ModuleType

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FROZEN_TRANSPORT = ROOT / 'analysis/response/transport.py'


def symmetric_root(D):
    """Return the unique symmetric PSD root of each finite symmetric 3x3 tensor.

    As in the frozen Gaussian update, eigenvalues below minus 1e-12 times
    the spectral scale are rejected; negative values within that roundoff
    tolerance are clipped to zero. A singular or zero PSD tensor is permitted.
    The root is continuous on PSD matrices; uniform Lipschitz continuity near
    a zero eigenvalue is not asserted.
    """
    D = np.asarray(D, dtype=float)
    if D.ndim != 3 or D.shape[1:] != (3, 3) or not np.isfinite(D).all():
        raise ValueError('diffusion tensors must be finite N by 3 by 3 arrays')
    scale = np.maximum(np.max(abs(D), axis=(1, 2)), np.finfo(float).tiny)
    transpose = D.transpose(0, 2, 1)
    if np.any(np.max(abs(D-transpose), axis=(1, 2)) > 1e-12*scale):
        raise ValueError('diffusion tensor is not symmetric')
    values, vectors = np.linalg.eigh((D+transpose)/2)
    scale = np.maximum(np.max(abs(values), axis=1), np.finfo(float).tiny)
    if np.any(values[:, 0] < -1e-12*scale):
        raise ValueError('diffusion tensor is not positive semidefinite')
    factor = vectors*np.sqrt(np.maximum(values, 0.))[:, None, :]
    return factor @ vectors.transpose(0, 2, 1)


def symmetric_gaussian_increment(A, D, elapsed, normal):
    """Produce A*delta + sqrt(delta)*D^(1/2)*Z without changing the drift."""
    D = np.asarray(D, dtype=float)
    root = symmetric_root(D)
    A, elapsed, normal = (np.asarray(value, dtype=float)
                          for value in (A, elapsed, normal))
    n = len(root)
    if (A.shape != (n, 3) or normal.shape != (n, 3)
            or elapsed.shape != (n,) or not np.isfinite(A).all()
            or not np.isfinite(normal).all() or not np.isfinite(elapsed).all()
            or np.any(elapsed < 0)):
        raise ValueError('finite N-row drift/noise and nonnegative durations required')
    return A*elapsed[:, None] + np.sqrt(elapsed)[:, None]*np.einsum(
        'nij,nj->ni', root, normal)


def make_transport(factor='eigen') -> ModuleType:
    """Load an independent complete solver using ``eigen`` or ``symmetric``.

    The source is never rewritten and this module is not installed under the
    canonical response-module name. Each call returns an isolated namespace,
    so simultaneous callers cannot change one another's transition function.
    The eigen arm keeps the exact original Gaussian function and RNG sequence.
    """
    if factor not in ('eigen', 'symmetric'):
        raise ValueError('factor must be eigen or symmetric')
    spec = spec_from_file_location('analysis.precision._frozen_transport',
                                   FROZEN_TRANSPORT)
    if spec is None or spec.loader is None:
        raise ImportError('Cannot load frozen response transport')
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    if factor == 'symmetric':
        module._gaussian_increment = symmetric_gaussian_increment
    module.precision_factor = factor
    return module
