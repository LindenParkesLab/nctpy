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


def modal_control(A_norm: npt.ArrayLike) -> FloatArray:
    """Compute the modal controllability of each node, for a discrete-time system.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix, for a discrete-time system.

    Returns
    -------
    phi : (N,) ndarray
        Modal controllability of each node: ``sum over j of U_ij^2 (1 - t_j^2)``, where ``A_norm = U T U^T`` is the
        real Schur decomposition and t = diag(T).

    Notes
    -----
    Modal controllability is defined for undirected (symmetric) connectomes, for which the Schur decomposition is
    the eigendecomposition: U holds the eigenvectors and t the eigenvalues. For a directed A it is not, so the
    values are an approximation and depend on the order of the nodes.
    """
    A_norm = _as_float(A_norm)
    T, U = schur(A_norm, "real")
    t = np.diag(T)
    return np.sum(U**2 * (1 - t**2), axis=1)
