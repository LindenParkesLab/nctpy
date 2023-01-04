import os
import numpy as np
import scipy as sp
from scipy import stats
from numpy.linalg import eig
from statsmodels.stats import multitest


def matrix_normalization(A, system=None, c=1):
    '''This function will normalize A in preparation for modeling linear dynamics.

    Args:
        A (NxN, numpy array): adjacency matrix representing a structural connectome.
        system (str): string variable that determines whether A is normalized for a continuous-time system or a
            discrete-time system. options: 'continuous' or 'discrete'. default=None.
        c (int): normalization constant, default=1.

    Returns:
        A_norm (NxN, numpy array): normalized adjacency matrix.

    '''

    if system is None:
        raise Exception("Time system not specified. "
                        "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                        "(see function help).")
    elif system != 'continuous' and system != 'discrete':
        raise Exception("Incorrect system specification. "
                        "Please specify either 'system=discrete' or 'system=continuous'.")
    else:
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
    """This function takes an array of integer values that designate a distinct set of binary brain states and returns
    a pair of matrices (x0_mat, xf_mat) that encode all possible pairwise transitions between those states.

    Args:
        states (N, numpy array): a vector of integers that designate which regions belong to which states.
            Note, regions cannot belong to more than one brain state.
            For example, assuming N = 12, if states = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2])
            then the first 4 regions belong to state 0, the next 4 to state 1, and the final 4 to state 2.

    Returns:
        x0_mat (Nxn_transitions, numpy boolean array): boolean array of initial states.
            In each column, True designates regions belonging to a given initial state.
        xf_mat (Nxn_transitions, numpy boolean array): boolean array of target states.
            In each column, True designates regions belonging to a given target state.

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
    """This function will normalize a brain state's magnitude using its euclidean norm.

    Args:
        x (N, numpy array): brain state to be normalized.

    Returns:
        x_norm (N, numpy array): normalized brain state.

    """

    x_norm = x / np.linalg.norm(x, ord=2)

    return x_norm


def normalize_weights(x, rank=True, add_constant=True):
    """This function normalizes the weights on B. By default, this involves (i) ranking the data, (ii) rescaling to
    the unit interval, and (iii) adding a constant value. If rank=False and add_constant=False, then only unit rescaling
    is performed.

    Args:
        x (N, numpy array): weights assigned to diagonal of NxN B matrix.
        rank: determines whether normalization includes ranking the data. default=True.
        add_constant: determines whether normalization includes adds a constant to the data. default=True.

    Returns:
        x (N, numpy array): normalized weights.

    """

    if rank:
        # rank data
        x = sp.stats.rankdata(x)

    # rescale to unit interval
    x = (x - min(x)) / (max(x) - min(x))

    if add_constant:
        x = x + 1

    return x


def get_null_p(x, null, version='standard', abs=False):
    """This function will compute p-values using an empirical null distribution.

    Args:
        x: observed test statistic.
        null: null distribution.
        version: determins how p-value will be computed, see below. default='standard'.
        abs: determines whether absolute values are taken for both x and null.

    Returns:
        p_val: p-value from null.

    """

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
    """This function will correct p-values for multiple comparisons using FDR.

    Args:
        p_vals (numpy array): array of p-values. Can be (n,) vector or (nxn) matrix.
        alpha (float): false discovery rate. default=0.05

    Returns:
        p_fdr (numpy array): corrected p-values.

    """

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


def convert_states_str2int(states_str):
    """This function takes a list of strings that designate a distinct set of binary brain states and returns
    a numpy array of integers encoding those states alongside a list of keys for those integers.

    Args:
        states_str (N, list): a list of strings that designate which regions belong to which states.
            For example, states = ['Vis', 'Vis', 'Vis', 'SomMot', 'SomMot', 'SomMot']

    Returns:
        states (N, numpy array): array of integers denoting which node belongs to which state.
        state_labels (n_states, list): list of keys corresponding to integers.
            For example, if state_labels[1] = 'SomMot' then the integer 1 in `states` corresponds to 'SomMot'.
            Together, a binary state can be extracted like so: x0 = states == state_labels.index('SomMot')

    """

    n_states = len(states_str)
    state_labels = list(np.unique(states_str))

    states = np.zeros(n_states)
    for i, state in enumerate(state_labels):
        for j in np.arange(n_states):
            if state == states_str[j]:
                states[j] = i

    return states.astype(int), state_labels


def expm(A):
    """This function computes the matrix exponential using eigen decomposition as per the Spectral Mapping Theorem

    Args:
        A: matrix to exponentiate

    Returns:
        eA: the matrix exponential
    """

    A_eig = np.linalg.eig(A)  # get eigenvalues and eigenvectors of A
    eA = np.diag(np.exp(A_eig[0]))  # put exponentiated eigenvalues on diagonal of a matrix
    eA = np.matmul(A_eig[1], eA)  # multiply by eigenvectors
    eA = np.matmul(eA, np.linalg.inv(A_eig[1]))  # multiple by inverse of eigenvectors
    eA = np.real(eA)  # retain real components

    return eA
