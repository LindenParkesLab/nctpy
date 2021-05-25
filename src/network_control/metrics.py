import numpy as np
from scipy.linalg import svd, schur

def ave_control(A, c=1):
    """ Returns values of AVERAGE CONTROLLABILITY for each node in a network, given the adjacency matrix for that
    network. Average controllability measures the ease by which input at that node can steer the system into many
    easily-reachable states.

    Args:
        A: np.array (n_parcels, n_parcels)
            Adjacency matrix from structural connectome
        c: int
            normalization constant, default=1

    Returns:
        ac: np.array (n_parcels,)
            vector of average controllability values for each node

    @author lindenmp
    Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
               Vettel, Miller, Grafton & Bassett, Nature Communications
               6:8414, 2015.
    """

    u, s, vt = svd(A)  # singluar value decomposition
    A = A / (c + s[0])  # Matrix normalization
    T, U = schur(A, 'real')  # Schur stability
    midMat = np.multiply(U, U).transpose()
    v = np.matrix(np.diag(T)).transpose()
    N = A.shape[0]
    P = np.diag(1 - np.matmul(v, v.transpose()))
    P = np.tile(P.reshape([N, 1]), (1, N))
    values = sum(np.divide(midMat, P))

    return values


def modal_control(A, c=1):
    """ Returns values of MODAL CONTROLLABILITY for each node in a network, given the adjacency matrix for that network.
    Modal controllability indicates the ability of that node to steer the system into difficult-to-reach states,
    given input at that node.

    Args:
        A: np.array (n_parcels, n_parcels)
            Adjacency matrix from structural connectome
        c: int
            normalization constant, default=1

    Returns:
        ac: np.array (n_parcels,)
            vector of modal controllability values for each node

    @author lindenmp
    Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
               Vettel, Miller, Grafton & Bassett, Nature Communications
               6:8414, 2015.
    """

    u, s, vt = svd(A)  # singluar value decomposition
    A = A / (c + s[0])  # Matrix normalization
    T, U = schur(A, 'real')  # Schur stability
    eigVals = np.diag(T)
    N = A.shape[0]
    phi = np.zeros(N, dtype=float)
    for i in range(N):
        Al = U[i,] * U[i,]
        Ar = (1.0 - np.power(eigVals, 2)).transpose()
        phi[i] = np.matmul(Al, Ar)

    return phi
