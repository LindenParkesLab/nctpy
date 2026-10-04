"""Controllability metrics for each node of a network.

Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia, Vettel, Miller, Grafton & Bassett,
Nature Communications 6:8414, 2015.
"""

from typing import cast

import numpy as np
import numpy.typing as npt
import scipy.linalg as la
from scipy.linalg import schur

from nctpy._validation import _check_system
from nctpy.energies import FloatArray, _as_float, gramian


def ave_control(A_norm: npt.ArrayLike, system: str | None = None) -> FloatArray:
    """Compute the average controllability of each node.

    The average controllability of node i is the trace of the controllability Gramian when input enters at
    node i alone (Gu et al. 2015): how much energy input at node i spreads into the network. With
    ``A[i, j]`` the connection from node j to node i, that is the sum over k >= 0 of ``||A^k e_i||^2`` in
    discrete time (an infinite horizon), and the integral over [0, 1] of ``||e^{At} e_i||^2`` in continuous
    time, where e_i is the i-th unit vector. Equivalently, the diagonal of the Gramian of ``A^T`` with B = I;
    for a symmetric (undirected) A, the diagonal of the Gramian of A.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system. Required.

    Returns
    -------
    ac : (N,) ndarray
        Average controllability of each node.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    """
    _check_system(system)
    A_norm = _as_float(A_norm)
    A_T = np.ascontiguousarray(A_norm.T)  # input at node i spreads through column i of A_norm

    if system == "continuous":
        return cast(FloatArray, gramian(A_T, T=1, system=system)).diagonal()
    # the infinite sum of (A^T)^k A^k solves X = A^T X A + I
    return np.diag(la.solve_discrete_lyapunov(A_T, np.eye(A_norm.shape[0])))


def modal_control(A_norm: npt.ArrayLike, *, system: str | None = "discrete") -> FloatArray:
    """Compute the modal controllability of each node.

    Modal controllability measures a node's ability to control the system's modes, weighted towards the
    fast-decaying ones, which are hard to reach. With ``A_norm = U T U^T`` its real Schur decomposition and
    t = diag(T):

    - discrete time (Gu et al. 2015, following Pasqualetti et al. 2014): ``sum over j of U_ij^2 (1 - t_j^2)``;
    - continuous time (Karrer et al. 2020, Eq. 8, adapted from the discrete form): ``sum over j of
      U_ij^2 (1 - exp(t_j))``.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix, normalised for the time system given by `system`.
    system : {'discrete', 'continuous'}, default 'discrete', keyword-only
        Time system that A_norm was normalised for. The default, discrete time, was the only option in earlier
        versions, and its results are unchanged.

    Returns
    -------
    phi : (N,) ndarray
        Modal controllability of each node. For a stable symmetric A_norm every value lies between 0 and 1.

    Raises
    ------
    Exception
        If `system` is None or not one of the two options.

    Notes
    -----
    Modal controllability is defined for undirected (symmetric) connectomes, for which the Schur decomposition is
    the eigendecomposition: U holds the eigenvectors and t the eigenvalues. For a directed A it is not, so the
    values are an approximation and depend on the order of the nodes.

    References
    ----------
    Gu, S. et al. Controllability of structural brain networks. Nat. Commun. 6, 8414 (2015).
    Karrer, T. M. et al. A practical guide to methodological considerations in the controllability of structural
    brain networks. J. Neural Eng. 17, 026031 (2020).
    """
    _check_system(system)
    A_norm = _as_float(A_norm)
    T, U = schur(A_norm, "real")
    t = np.diag(T)
    if system == "continuous":
        return np.sum(U**2 * (1 - np.exp(t)), axis=1)
    return np.sum(U**2 * (1 - t**2), axis=1)
