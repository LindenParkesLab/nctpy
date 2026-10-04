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


def modal_control(A_norm: npt.ArrayLike, *, system: str | None = "discrete") -> FloatArray:
    """Compute the modal controllability of each node.

    Modal controllability measures a node's ability to control the system's modes, weighted towards the
    fast-decaying ones, which are hard to reach. With ``A_norm = U T U^T`` its real Schur decomposition and
    t = diag(T):

    - discrete time (Gu et al. 2015, following Pasqualetti et al. 2014): ``sum over j of U_ij^2 (1 - t_j^2)``;
    - continuous time (Karrer et al. 2020, Eq. 8, adapted from the discrete form): ``sum over j of
      U_ij^2 (1 - exp(t_j))``.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix, normalised for the time system given by `system`.
    system : {'discrete', 'continuous'}, default 'discrete', keyword-only
        Time system that A_norm was normalised for. The default, discrete time, was the only option in earlier
        versions, and its results are unchanged.

    Returns
    -------
    phi : (N,) ndarray
        Modal controllability of each node. For a stable symmetric A_norm every value lies between 0 and 1.

    Raises
    ------
    Exception
        If `system` is None or not one of the two options.

    Notes
    -----
    Modal controllability is defined for undirected (symmetric) connectomes, for which the Schur decomposition is
    the eigendecomposition: U holds the eigenvectors and t the eigenvalues. For a directed A it is not, so the
    values are an approximation and depend on the order of the nodes.

    References
    ----------
    Gu, S. et al. Controllability of structural brain networks. Nat. Commun. 6, 8414 (2015).
    Karrer, T. M. et al. A practical guide to methodological considerations in the controllability of structural
    brain networks. J. Neural Eng. 17, 026031 (2020).
    """
    _check_system(system)
    A_norm = _as_float(A_norm)
    T, U = schur(A_norm, "real")
    t = np.diag(T)
    if system == "continuous":
        return np.sum(U**2 * (1 - np.exp(t)), axis=1)
    return np.sum(U**2 * (1 - t**2), axis=1)


def boundary_control(
    A: npt.ArrayLike, partition: npt.ArrayLike, system: str | None = None, rho: float = 0.2
) -> FloatArray:
    """Compute the boundary controllability of each node (Gu et al. 2015).

    Boundary controllability ranks nodes by how early they are found on the boundaries between network
    communities, as the communities are split recursively. It follows Gu et al.'s released code
    (github.com/gushiapi/Network-Controllability-Diagnostics, ``getTieredControlStatesSH.m``):

    1. Scale A so that its largest weight is 1. A node is a first-level boundary node if its total weight to
       nodes outside its community is at least `rho`.
    2. Among the communities that still contain nodes without a value, select the least controllable one: the
       community whose nodes, used as control nodes (B = their columns of I), give the infinite-horizon
       controllability Gramian with the smallest eigenvalue.
    3. Split it in two by the sign of the Fiedler vector of its Laplacian. Its boundary nodes are those whose total
       weight to the other half is at least `rho` times the largest weight within the community. The first half
       becomes a new community.
    4. Repeat 2-3 until every node has a value. Nodes found earlier rank higher: the first-level boundary nodes
       get 1, the nodes of the next tier (N - N1)/N, where N1 nodes were found before them, and so on.

    Parameters
    ----------
    A : (N, N) array_like
        Structural connectome, **not** normalised: the algorithm scales and normalises A itself, as Gu et al.'s
        code does, so the result does not depend on the scale of A.
    partition : (N,) array_like
        Community of each node, for example functional systems or a modularity partition. Labels are taken in
        sorted order, which sets the order in which ties are broken.
    system : {'discrete', 'continuous'}
        Time system of the Gramian in step 2. Required. With W = A / max(A) and lambda its spectral radius, the
        system matrix is W / (1 + lambda) in discrete time (Gu et al.) and W / (1 + lambda) - I in continuous time
        (the continuous form of Gu et al.'s Supplementary Methods).
    rho : float, default 0.2
        Threshold ratio for boundary nodes (Gu et al. use 0.2).

    Returns
    -------
    bc : (N,) ndarray
        Boundary controllability of each node, between 0 and 1; 1 for the first-level boundary nodes.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    ValueError
        If `partition` does not have one entry per node.

    Notes
    -----
    This is an implementation of Gu et al.'s released code, and reproduces its behaviour, including:

    - **Selection at rounding level.** The Gramians compared in step 2 are badly conditioned (as Gu et al. note), and
      in large networks their smallest eigenvalues are often below what double precision resolves. Which community is
      split next, and so the order of the tiers, can then depend on rounding: an input changed by rounding error, or
      the same computation with another linear algebra library, can give different values. The arithmetic here
      follows Gu et al.'s code, so it reproduces the code's choices given the same rounding.
    - **Re-marking.** A node found again in a later split is moved to the newest tier, as in the code, so a
      first-level boundary node can end with a value below 1. The paper's text assigns values only to the additional
      boundary nodes of each step.
    - **Parameters.** Gu et al. chose rho = 0.2 for their parcellation (234 regions); at other resolutions, or with
      other edge weights, it may mark many nodes as first-level boundary nodes.

    The Gramian is computed in the Schur basis of the system matrix, as in Gu et al.'s ``Gramian.m``, which is exact
    for a symmetric (undirected) connectome and an approximation for a directed one. Gu et al. stabilise A by
    1 + lambda after scaling it to a largest weight of 1, as in their code; their paper's text describes dividing by
    the mean edge weight instead.

    References
    ----------
    Gu, S. et al. Controllability of structural brain networks. Nat. Commun. 6, 8414 (2015).
    Pasqualetti, F., Zampieri, S. & Bullo, F. Controllability metrics, limitations and algorithms for complex
    networks. IEEE Trans. Control Netw. Syst. 1, 40-52 (2014).
    """
    _check_system(system)
    W = _as_float(A)
    n_nodes = W.shape[0]
    labels = np.asarray(partition)
    if labels.shape != (n_nodes,):
        raise ValueError(f"partition must have one entry per node ({n_nodes}); got shape {labels.shape}")
    community = np.unique(labels, return_inverse=True)[1].reshape(-1)  # 0, 1, ... in sorted label order

    W = W / np.max(W)
    # system matrix in its Schur basis, A_sys = U T U^T (Gu et al.'s Gramian.m)
    A_sys = W / (1 + np.max(np.abs(np.linalg.eigvals(W))))
    if system == "continuous":
        A_sys = A_sys - np.eye(n_nodes)
    T, U = schur(A_sys, "real")
    v = np.diag(T)
    # infinite-horizon Gramian in that basis: discrete sum_k A^k B B^T A^kT, continuous integral of e^{At} B B^T e^{A^Tt}
    denominator = 1 - np.outer(v, v) if system == "discrete" else -np.add.outer(v, v)

    def smallest_gramian_eigenvalue(nodes: npt.NDArray[np.int_]) -> float:
        # U^T B B^T U, with B the columns of I for these nodes, formed as Gu et al.'s code forms it: the smallest
        # eigenvalue can be at the level of rounding error, so the order of operations decides which community is
        # selected (multiplying by the 0/1 mask is exact, so this rounds exactly as U' * (B * B') * U)
        mask = np.zeros(n_nodes)
        mask[nodes] = 1.0
        return float(np.linalg.svd(((U.T * mask) @ U) / denominator, compute_uv=False)[-1])

    # 1. first-level boundary nodes: weight to other communities >= rho (W has largest weight 1)
    tier = np.zeros(n_nodes, dtype=int)  # 1 = most recently found; 0 = not yet found
    for c in np.unique(community):
        inside = community == c
        tier[np.flatnonzero(inside)[W[np.ix_(inside, ~inside)].sum(axis=1) >= rho]] = 1

    eigenvalue: dict[int, float] = {}  # cached per community; a community's nodes change only when it is split
    while np.count_nonzero(tier) < n_nodes:
        candidates = [c for c in np.unique(community) if not np.all(tier[community == c] > 0)]
        for c in candidates:
            if c not in eigenvalue:
                eigenvalue[c] = smallest_gramian_eigenvalue(np.flatnonzero(community == c))
        selected = min(candidates, key=lambda c: eigenvalue[c])  # first of equal minima, like MATLAB's min
        nodes = np.flatnonzero(community == selected)

        boundary, first_half = _fiedler_boundary(W, nodes, rho)
        if boundary.size:
            tier[tier > 0] += 1
            tier[boundary] = 1
        if first_half.size:
            new = int(np.max(community)) + 1
            community[first_half] = new
            del eigenvalue[selected]  # both halves have changed

    # value of tier k: the fraction of nodes found in tiers <= k (Gu et al.'s counts, divided by N)
    counts = np.array([np.count_nonzero(tier <= k) for k in range(np.max(tier) + 1)])
    return counts[tier] / n_nodes


def _fiedler_boundary(
    W: FloatArray, nodes: npt.NDArray[np.int_], rho: float
) -> tuple[npt.NDArray[np.int_], npt.NDArray[np.int_]]:
    """Split a community by its Fiedler vector; return its boundary nodes and the nodes of the first half.

    As Gu et al.'s ``boundaryDetection.m``: a single node is its own boundary; otherwise the halves are the signs of
    the Laplacian's second-smallest singular vector (split by index if one half is empty), and a node is on the
    boundary if its total weight to the other half is at least rho times the community's largest weight.
    """
    if nodes.size == 1:
        return nodes, nodes[:0]
    sub = W[np.ix_(nodes, nodes)]
    threshold = np.max(sub) * rho
    laplacian = np.diag(sub.sum(axis=0)) - sub
    fiedler = np.linalg.svd(laplacian)[0][:, -2]
    first, second = np.flatnonzero(fiedler >= 0), np.flatnonzero(fiedler < 0)
    if first.size == 0 or second.size == 0:
        middle = (1 + nodes.size) // 2
        first, second = np.arange(middle), np.arange(middle, nodes.size)
    on_boundary = np.zeros(nodes.size, dtype=bool)
    on_boundary[first] = sub[np.ix_(first, second)].sum(axis=1) >= threshold
    on_boundary[second] = sub[np.ix_(second, first)].sum(axis=1) >= threshold
    return nodes[on_boundary], nodes[first]
