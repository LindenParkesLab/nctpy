import os
import numpy as np
import scipy as sp
from scipy import stats
from numpy.linalg import eig
from statsmodels.stats import multitest


def matrix_normalization(A, system=None, c=1):
    '''

    Args:
        A: np.array (n_parcels, n_parcels)
            adjacency matrix from structural connectome
        system: str
            options: 'continuous' or 'discrete'. default=None
            string variable that determines whether A is normalized for a continuous-time system or a discrete-time
            system. If normalizing for a continuous-time system, the identity matrix is subtracted.
        c: int
            normalization constant, default=1
    Returns:
        A_norm: np.array (n_parcels, n_parcels)
            normalized adjacency matrix

    '''

    if system == None:
        raise Exception("Time system not specified. "
                        "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                        "(see function help).")

    # eigenvalue decomposition
    w, _ = eig(A)
    l = np.abs(w).max()

    # Matrix normalization for discrete-time systems
    A_norm = A / (c + l)

    if system == 'continuous':
        # for continuous-time systems
        A_norm = A_norm - np.eye(A.shape[0])

    return A_norm


def get_p_val_string(p_val):
    if p_val == 0.0:
        p_str = "-log10($\mathit{:}$)>25".format('{p}')
    elif p_val < 0.05:
        p_str = '$\mathit{:}$ = {:0.0e}'.format('{p}', p_val)
    else:
        p_str = "$\mathit{:}$ = {:.3f}".format('{p}', p_val)

    return p_str


def expand_states(states):
    """
    This function takes a list of integer values that designate a distinct set of binary brain states and returns
    a pair of matrices (x0_mat, xf_mat) that encode all possible pairwise transitions between those states
    Args:
        states: numpy array (N x 1)
            a vector of integers that designate which regions belong to which states. Note, regions cannot belong to
            more than one brain state. For example, assuming N = 12, if:
                states = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2])
            then the first 4 regions belong to state 0, the next 4 to state 1, and the final 4 to state 2

    Returns:
        x0_mat: boolean array (N, n_transitions)
            boolean array of initial states. In each column, True designates regions belonging to a given initial state
        xf_mat: boolean array (N, n_transitions)
            boolean array of target states. In each column, True designates regions belonging to a given target state
    """

    unique, counts = np.unique(states, return_counts=True)
    n_parcels = len(states)
    n_states = len(unique)

    x0_mat = np.zeros((n_parcels, 1)).astype(bool)
    xf_mat = np.zeros((n_parcels, 1)).astype(bool)

    for i in np.arange(n_states):
        for j in np.arange(n_states):
            x0 = states == i
            xf = states == j

            x0_mat = np.append(x0_mat, x0.reshape(-1, 1), axis=1)
            xf_mat = np.append(xf_mat, xf.reshape(-1, 1), axis=1)

    x0_mat = x0_mat[:, 1:]
    xf_mat = xf_mat[:, 1:]

    return x0_mat, xf_mat


def normalize_state(x):
    x_norm = x / np.linalg.norm(x, ord=2)

    return x_norm


def normalize_weights(x, rank=True, add_constant=True):
    if rank:
        # rank data
        x = sp.stats.rankdata(x)

    # rescale to unit interval
    x = (x - min(x)) / (max(x) - min(x))

    if add_constant:
        x = x + 1

    return x


def get_null_p(x, null, version='standard', abs=False):
    if abs:
        x = np.abs(x)
        null = np.abs(null)

    if version == 'standard':
        p_val = np.sum(null >= x) / len(null)
    elif version == 'reverse':
        p_val = np.sum(x >= null) / len(null)
    elif version == 'smallest':
        p_val = np.min([np.sum(null >= x) / len(null),
                        np.sum(x >= null) / len(null)])

    return p_val


def get_fdr_p(p_vals, alpha=0.05):
    if p_vals.ndim == 2:
        do_reshape = True
        dims = p_vals.shape
        p_vals = p_vals.flatten()
    else:
        do_reshape = False

    out = multitest.multipletests(p_vals, alpha=alpha, method='fdr_bh')
    p_fdr = out[1]

    if do_reshape:
        p_fdr = p_fdr.reshape(dims)

    return p_fdr


def convert_states_str2float(states_str):
    n = len(states_str)
    state_labels = list(np.unique(states_str))

    states = np.zeros(n)
    for i, state in enumerate(state_labels):
        for j in np.arange(n):
            if state == states_str[j]:
                states[j] = i

    return states.astype(int), state_labels
