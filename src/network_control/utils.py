import numpy as np
import scipy as sp
from scipy.linalg import svd

def rank_to_normal(data, c, n):
    # Standard quantile function
    data = (data - c) / (n - 2 * c + 1)
    return sp.stats.norm.ppf(data)


def rank_int(data, c=3.0 / 8):
    if data.ndim > 1:
        do_reshape = True
        dims = data.shape
        data = data.flatten()
    else:
        do_reshape = False

    # Set seed
    np.random.seed(0)

    # Get rank, ties are averaged
    data = sp.stats.rankdata(data, method="average")

    # Convert rank to normal distribution
    transformed = rank_to_normal(data=data, c=c, n=len(data))

    if do_reshape:
        transformed = transformed.reshape(dims)

    return transformed


def matrix_normalization(A, c=1):
    '''

    Args:
        A: np.array (n_parcels, n_parcels)
            adjacency matrix from structural connectome
        c: int
            normalization constant, default=1

    Returns:
        A_norm: np.array (n_parcels, n_parcels)
            normalized adjacency matrix

    '''
    u, s, vt = svd(A)  # singluar value decomposition
    A_norm = A / (c + s[0])  # Matrix normalization

    return A_norm


def get_p_val_string(p_val):
    if p_val == 0.0:
        p_str = "-log10($\mathit{:}$)>25".format('{p}')
    elif p_val < 0.001:
        p_str = '$\mathit{:}$ < 0.001'.format('{p}')
    elif p_val >= 0.001 and p_val < 0.05:
        p_str = '$\mathit{:}$ < 0.05'.format('{p}')
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
