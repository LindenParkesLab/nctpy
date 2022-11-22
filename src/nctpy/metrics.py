import numpy as np
from scipy.linalg import schur
from nctpy.energies import gramian


def ave_control(A_norm, system=None):
    """ Returns values of AVERAGE CONTROLLABILITY for each node in a network, given the adjacency matrix for that
    network. Average controllability measures the ease by which input at that node can steer the system into many
    easily-reachable states.

    Args:
        A_norm: np.array (n_parcels, n_parcels)
            Normalized adjacency matrix from structural connectome (see matrix_normalization in utils for example)
        system: str
            options: 'continuous' or 'discrete'. default=None

    Returns:
        ac: np.array (n_parcels,)
            vector of average controllability values for each node

    @author lindenmp
    Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
               Vettel, Miller, Grafton & Bassett, Nature Communications
               6:8414, 2015.
    """

    if system is None:
        raise Exception("Time system not specified. "
                        "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                        "(see matrix_normalization help).")
    elif system != 'continuous' and system != 'discrete':
        raise Exception("Incorrect system specification. "
                        "Please specify either 'system=discrete' or 'system=continuous'.")
    elif system == 'discrete':
        T, U = schur(A_norm, 'real')  # Schur stability
        midMat = np.multiply(U, U).transpose()
        v = np.diag(T)[np.newaxis, :].transpose()
        N = A_norm.shape[0]
        P = np.diag(1 - np.matmul(v, v.transpose()))
        P = np.tile(P.reshape([N, 1]), (1, N))
        ac = sum(np.divide(midMat, P))

        return ac
    elif system == 'continuous':
        G = gramian(A_norm, T=1, system=system)
        ac = G.diagonal()

        return ac


def modal_control(A_norm):
    """ Returns values of MODAL CONTROLLABILITY for each node in a network, given the adjacency matrix for that network.
    Modal controllability indicates the ability of that node to steer the system into difficult-to-reach states,
    given input at that node. Expects input to be a DISCRETE system

    Args:
        A_norm: np.array (n_parcels, n_parcels)
            Normalized adjacency matrix from structural connectome (see matrix_normalization in utils for example)

    Returns:
        phi: np.array (n_parcels,)
            vector of modal controllability values for each node

    @author lindenmp
    Reference: Gu, Pasqualetti, Cieslak, Telesford, Yu, Kahn, Medaglia,
               Vettel, Miller, Grafton & Bassett, Nature Communications
               6:8414, 2015.
    """

    T, U = schur(A_norm, 'real')  # Schur stability
    eigVals = np.diag(T)
    N = A_norm.shape[0]
    phi = np.zeros(N, dtype=float)
    for i in range(N):
        Al = U[i, ] * U[i, ]
        Ar = (1.0 - np.power(eigVals, 2)).transpose()
        phi[i] = np.matmul(Al, Ar)

    return phi