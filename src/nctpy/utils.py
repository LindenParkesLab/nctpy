"""Normalisation, brain states, null-model p-values and other helpers."""

from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests

from nctpy._validation import _check_system


def matrix_normalization(
    A: npt.ArrayLike,
    system: str | None = None,
    c: float = 1,
    l: float | None = None,
    *,
    zero_diagonal: bool = False,
    decay: npt.ArrayLike | None = None,
) -> npt.NDArray[np.float64]:
    """Normalise a structural connectome A for modelling linear dynamics.

    ``A_norm = A / (c + l) - diag(decay)`` for continuous-time systems (Kim et al., Nat Commun 2025, Eq. 4) and
    ``A_norm = A / (c + l)`` for discrete-time systems. By default l is the spectral radius of A (its largest
    absolute eigenvalue) and decay is 1 at every node, i.e. ``A / (c + l) - I``.

    Parameters
    ----------
    A : (N, N) array_like
        Adjacency matrix representing a structural connectome.
    system : {'continuous', 'discrete'}
        Time system to normalise A for. Required.
    c : float, default 1
        Normalisation constant.
    l : float, optional
        Fixed spectral radius to normalise by, in place of A's own. Use one l across several connectomes (e.g.
        subjects) so that they share a normalisation, typically the maximum spectral radius over all of them.

        The normalised system is guaranteed to be stable when c + l exceeds the spectral radius of A, which the
        default l always satisfies for c > 0. An l below A's own spectral radius forfeits that guarantee. No check
        is made: an unstable system still returns values.
    zero_diagonal : bool, default False, keyword-only
        Set the diagonal of (a copy of) A to zero before normalising. The model assumes A has no self-connections:
        each node's own dynamics come from the normalisation (in continuous time, the ``- I`` term), not from
        diag(A). A connectome with a non-zero diagonal is otherwise used as given; nctpy never removes
        self-connections unless asked. The protocol paper's code uses the PNC connectome with its diagonal intact,
        so the default stays False. The array passed in is never modified.
    decay : float or (N,) array_like, optional, keyword-only
        Continuous time only: the decay rate of each node, subtracted from the diagonal (Eq. 4), so a larger value
        means stronger self-inhibition (dx_i/dt = -decay_i x_i + ...). A scalar applies to every node. The default,
        None, is 1 at every node, as in the protocol paper. With the default l and c > 0, a uniform decay of at
        least 1 keeps the system stable, as does decay >= 1 at every node for a symmetric A; smaller (or negative)
        rates forfeit that guarantee.
        Nothing is checked: an unstable system still returns values.

        Without ``decay``, self-connections already make the effective decay rate non-uniform: node i's is
        ``1 - A_ii / (c + l)``. Decay rates fitted by optimising them (Kim et al., 2025) assume a connectome without
        self-connections, so that every node starts from the same baseline. A value v reported in that paper's
        Fig. 2B corresponds to ``decay = 1 - v`` here: the fitted diagonal of A_norm is v - 1.

    Returns
    -------
    A_norm : (N, N) ndarray
        Normalised adjacency matrix.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    ValueError
        If `decay` is given for a discrete-time system, or is neither a scalar nor one value per node.
    """
    _check_system(system, see="function help")
    if decay is not None and system == "discrete":
        raise ValueError("decay applies to continuous-time systems only (Eq. 4); discrete-time normalisation has none")
    A = np.asarray(A)
    A = A.astype(np.result_type(A.dtype, np.float64), copy=False)
    if zero_diagonal:
        A = A.copy()
        np.fill_diagonal(A, 0)
    if l is None:
        l = np.abs(np.linalg.eig(A)[0]).max()
    A_norm = A / (c + l)
    if system == "continuous":
        A_norm = A_norm - (np.eye(A.shape[0]) if decay is None else np.diag(_decay_rates(decay, A.shape[0])))
    return A_norm


def _decay_rates(decay: npt.ArrayLike, n_nodes: int) -> npt.NDArray[np.float64]:
    """decay as one float per node: a scalar is broadcast; anything else must have one value per node."""
    rates = np.asarray(decay, dtype=float)
    if rates.ndim == 0:
        return np.full(n_nodes, float(rates))
    if rates.shape != (n_nodes,):
        raise ValueError(f"decay must be a scalar or have one value per node ({n_nodes},); got shape {rates.shape}")
    return rates


def get_p_val_string(p_val: float) -> str:
    """Format a p-value for a matplotlib label: '-log10(p)>25' for 0, scientific below 0.05, else 3 decimals."""
    if p_val == 0.0:
        return r"-log10($\mathit{p}$)>25"
    if p_val < 0.05:
        return rf"$\mathit{{p}}$ = {p_val:0.0e}"
    return rf"$\mathit{{p}}$ = {p_val:.3f}"


def expand_states(states: npt.ArrayLike) -> tuple[npt.NDArray[np.bool_], npt.NDArray[np.bool_]]:
    """Encode every pairwise transition between a set of binary brain states.

    Parameters
    ----------
    states : (N,) array_like of int
        The state each region belongs to, numbered 0 to n_states - 1. Regions cannot belong to more than one
        state. For example, with N = 12, ``states = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2])`` puts the first
        4 regions in state 0, the next 4 in state 1, and the final 4 in state 2.

    Returns
    -------
    x0_mat : (N, n_states**2) ndarray of bool
        Initial states, one transition per column: True marks the regions in that transition's initial state.
    xf_mat : (N, n_states**2) ndarray of bool
        Target states, in the same column order. Column ``i * n_states + j`` is the transition from state i to
        state j, self-transitions included.
    """
    states = np.asarray(states)
    n_states = len(np.unique(states))
    labels = np.arange(n_states)
    x0_mat = states[:, np.newaxis] == np.repeat(labels, n_states)
    xf_mat = states[:, np.newaxis] == np.tile(labels, n_states)
    return x0_mat, xf_mat


def normalize_state(x: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Scale a brain state to unit Euclidean norm.

    Parameters
    ----------
    x : (N,) array_like
        Brain state. Boolean states are accepted and become floats.

    Returns
    -------
    x_norm : (N,) ndarray
        The state divided by its Euclidean norm.
    """
    x = np.asarray(x)
    return x / np.linalg.norm(x, ord=2)


def normalize_weights(x: npt.ArrayLike, rank: bool = True, add_constant: bool = True) -> npt.NDArray[np.float64]:
    """Normalise weights for the diagonal of B.

    By default: (i) rank the data, (ii) rescale to the unit interval, and (iii) add 1, giving weights on [1, 2].
    With rank=False and add_constant=False, only the rescaling is done.

    Parameters
    ----------
    x : (N,) array_like
        Weights for the diagonal of an N x N B matrix.
    rank : bool, default True
        Rank the data first.
    add_constant : bool, default True
        Add 1 after rescaling.

    Returns
    -------
    x : (N,) ndarray
        Normalised weights.
    """
    w = rankdata(x) if rank else np.asarray(x)
    w = (w - min(w)) / (max(w) - min(w))
    if add_constant:
        w = w + 1
    return w


def get_null_p(x: Any, null: npt.ArrayLike, version: str = "standard", abs: bool = False) -> float:
    """Compute a p-value from an empirical null distribution.

    Parameters
    ----------
    x : float
        Observed test statistic.
    null : (n,) array_like
        Null distribution.
    version : {'standard', 'reverse', 'smallest'}, default 'standard'
        'standard' is the fraction of the null at or above x; 'reverse' the fraction at or below x; 'smallest' the
        smaller of the two.
    abs : bool, default False
        Take absolute values of both x and the null first.

    Returns
    -------
    p_val : float

    Raises
    ------
    ValueError
        If `version` is not one of the three options.
    """
    null_arr = np.abs(null) if abs else np.asarray(null)
    if abs:
        x = np.abs(x)

    upper = np.sum(null_arr >= x) / len(null_arr)
    lower = np.sum(x >= null_arr) / len(null_arr)
    if version == "standard":
        return upper
    if version == "reverse":
        return lower
    if version == "smallest":
        return np.min([upper, lower])
    raise ValueError(f"version must be 'standard', 'reverse' or 'smallest', got {version!r}")


def get_fdr_p(p_vals: npt.ArrayLike, alpha: float = 0.05) -> npt.NDArray[np.float64]:
    """Correct p-values for multiple comparisons with the Benjamini-Hochberg false discovery rate.

    Parameters
    ----------
    p_vals : array_like
        p-values, of any shape (e.g. a vector, or a matrix of transitions). All are corrected together.
    alpha : float, default 0.05
        False discovery rate.

    Returns
    -------
    p_fdr : ndarray
        Corrected p-values, in the same shape as `p_vals`.
    """
    p_vals = np.asarray(p_vals)
    return multipletests(p_vals.ravel(), alpha=alpha, method="fdr_bh")[1].reshape(p_vals.shape)


def convert_states_str2int(states_str: Any) -> tuple[npt.NDArray[np.int_], list[Any]]:
    """Encode a list of state names as integers.

    Parameters
    ----------
    states_str : (N,) list of str
        The state each region belongs to, e.g. ``['Vis', 'Vis', 'Vis', 'SomMot', 'SomMot', 'SomMot']``.

    Returns
    -------
    states : (N,) ndarray of int
        The integer code of each region's state.
    state_labels : list
        The state names in alphabetical order; the integer i stands for ``state_labels[i]``. A binary state
        can be extracted like so: ``x0 = states == state_labels.index('SomMot')``.
    """
    labels, states = np.unique(states_str, return_inverse=True)
    return states.reshape(-1).astype(int), list(labels)


def expm(A: npt.ArrayLike) -> npt.NDArray[np.float64]:
    """Compute the matrix exponential by eigendecomposition (spectral mapping theorem), keeping the real part.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix to exponentiate; it must be diagonalisable.

    Returns
    -------
    eA : (N, N) ndarray
        The real part of ``V diag(exp(w)) V^-1``, where A = V diag(w) V^-1.
    """
    w, V = np.linalg.eig(np.asarray(A))
    return np.real(V @ np.diag(np.exp(w)) @ np.linalg.inv(V))
