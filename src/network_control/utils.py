import numpy as np
import scipy as sp
import scipy.linalg as la
from scipy.linalg import svd
from scipy.linalg import eig
from numpy import matmul as mm
from scipy.linalg import expm as expm
from numpy import transpose as tp


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


def matrix_normalization(A, version=None, c=1):
    '''

    Args:
        A: np.array (n_parcels, n_parcels)
            adjacency matrix from structural connectome
        version: str
            options: 'continuous' or 'discrete'. default=None
            string variable that determines whether A is normalized for a continuous-time system or a discrete-time
            system. If normalizing for a continuous-time system, the identity matrix is subtracted.
        c: int
            normalization constant, default=1
    Returns:
        A_norm: np.array (n_parcels, n_parcels)
            normalized adjacency matrix

    '''

    if version == 'continuous':
        print("Normalizing A for a continuous-time system")
    elif version == 'discrete':
        print("Normalizing A for a discrete-time system")
    elif version == None:
        raise Exception("Time system not specified. "
                        "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                        "(see function help).")

    # singluar value decomposition
    u, s, vt = svd(A)

    # Matrix normalization for discrete-time systems
    A_norm = A / (c + s[0])

    if version == 'continuous':
        # for continuous-time systems
        A_norm = A_norm - np.eye(A.shape[0])

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

