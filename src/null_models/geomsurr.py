import numpy as np


def rank_reorder(x, scaffold):
    """
    Original authors: M Breakspear, J Roberts
    Translated to Python by Linden Parkes

    If you use this code, please cite the original publication: Roberts et al. (2016) NeuroImage 124:379-393.
    """

    # reorder vector x to have same rank order as vector scaffold
    y = np.zeros((len(x), 2))
    y[:, 0] = scaffold
    y[:, 1] = np.arange(len(x))
    y = y[np.argsort(y[:, 0]), :]
    y[:, 0] = np.sort(x)
    y = y[np.argsort(y[:, 1]), :]

    return y[:, 0]


def strength_correct(W, ss, nreps=9):
    """
    Original authors: M Breakspear, J Roberts
    Translated to Python by Linden Parkes

    If you use this code, please cite the original publication: Roberts et al. (2016) NeuroImage 124:379-393.
    """

    N = len(ss)
    discnodes = np.sum(W, axis=0) == 0
    sW = np.multiply(W, np.repeat(np.divide(ss, np.sum(W, axis=0))[:, np.newaxis], N, axis=1))
    sW[discnodes, :] = 0  # ensure disconnected nodes don't introduce NaNs
    sW = sW + sW.T
    sW = sW / 2
    counter = 0
    while counter < nreps:
        sW = np.multiply(sW, np.repeat(np.divide(ss, np.sum(sW, axis=0))[:, np.newaxis], N, axis=1))
        sW[discnodes, :] = 0
        sW = sW + sW.T
        sW = sW / 2

        counter += 1

    return sW


def geomsurr(W, D, nmean=3, nstd=2, seed=123):
    """
    Random graphs that preserve distance effect
    Note: Wsp and Wssp are generated assuming that W is undirected.

    Original authors: M Breakspear, J Roberts
    Translated to Python by Linden Parkes

    If you use this code, please cite the original publication: Roberts et al. (2016) NeuroImage 124:379-393.

    :param W:
    :param D:
    :param nmean:
    :param nstd:
    :return:
    """

    # set state
    np.random.seed(seed)

    # Check if directed
    drct = 1
    if np.max(W-W.T) == 0:
        drct = 0

    # Preliminaries
    N = W.shape[0]
    Wwp = np.zeros((N, N))

    # ensure no self-connections
    W[np.eye(N) == 1] = 0

    # Only do one triangle if undirected, then replicate
    if drct == 0:
        W = np.tril(W)

    nz = np.where(W != 0)  # Nonzero entries
    w = W[nz]  # Vector of non-zero connections
    d = D[nz]  # Corresponding distances
    logw = np.log(w)  # log-weights

    # 1. remove mean to nmean order
    pfit1 = np.polyfit(d, logw, nmean)
    mnlogw = logw - np.polyval(pfit1, d)

    # 2. adjust variance to nstd order
    pfit2 = np.polyfit(d, np.abs(mnlogw), nstd)
    stdlogw = mnlogw / np.polyval(pfit2, d)

    # 3. Now create surrogate data, adjusted for mean and std
    # Shuffle the old ones
    surr = np.random.permutation(stdlogw)

    # 4. Now put the geometry back in
    # 4.1 Invert
    stdsurr = np.multiply(surr, np.polyval(pfit2, d))  # std
    mnsurr = stdsurr + np.polyval(pfit1, d)  # and mean

    # 4.2 Use surrogate weights as a scaffold to reorder the original weights,
    #    thus preserving the original set of weights but in the new distance-
    #    preserving random order (cf. the "amplitude adjustment" used in
    #    Fourier surrogates for time series analysis)
    surrlogw = rank_reorder(logw, mnsurr)

    # 4.3 Undo logarithm and put into weight-preserving surrogate matrix
    Wwp[nz] = np.exp(surrlogw)

    # Make undirected if W is undirected
    if drct == 0:
        Wwp = Wwp + Wwp.T
        W = W + W.T

    # 5. Adjust node strengths - NOTE: assumes W is undirected
    strengthsW = np.sum(W, axis=0)  # original node strengths
    strengthsWwp = np.sum(Wwp, axis=0)  # new node strengths

    # 5.1 Re-order the old strengths to match new random order in Wwp
    strengthsWsp = rank_reorder(strengthsW, strengthsWwp)

    # 5.2 Adjust strengths to give both original and random strength sequences
    Wsp = strength_correct(Wwp, strengthsWsp)  # orig strengths in new sequence
    Wssp = strength_correct(Wwp, strengthsW)  # orig strengths in orig sequence

    return Wwp, Wsp, Wssp
