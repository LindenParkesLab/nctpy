"""Controllability metrics for each node of a network.

Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia, Vettel, Miller, Grafton & Bassett,
Nature Communications 6:8414, 2015.
"""

from typing import cast

import numpy as np
import numpy.typing as npt
from scipy.linalg import schur

from nctpy._validation import _check_system
from nctpy.energies import FloatArray, _as_float, gramian


def ave_control(A_norm: npt.ArrayLike, system: str | None = None) -> FloatArray:
    """Compute the average controllability of each node.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system. Required.

    Returns
    -------
    ac : (N,) ndarray
        Average controllability of each node. In continuous time, the diagonal of the controllability Gramian
        over T = 1 (see :func:`nctpy.energies.gramian`). In discrete time, ``sum over j of U_ij^2 / (1 - t_j^2)``,
        where ``A_norm = U T U^T`` is the real Schur decomposition and t = diag(T): the diagonal of the
        infinite-horizon Gramian when A_norm is symmetric.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    """
    _check_system(system)
    A_norm = _as_float(A_norm)

    if system == "continuous":
        return cast(FloatArray, gramian(A_norm, T=1, system=system)).diagonal()

    T, U = schur(A_norm, "real")
    t = np.diag(T)
    return np.sum(U**2 / (1 - t**2), axis=1)


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
    """
    A_norm = _as_float(A_norm)
    T, U = schur(A_norm, "real")
    t = np.diag(T)
    return np.sum(U**2 * (1 - t**2), axis=1)
