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