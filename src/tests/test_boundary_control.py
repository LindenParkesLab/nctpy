"""nctpy.metrics.boundary_control: boundary controllability (Gu et al., Nat Commun 2015).

The reference is Gu et al.'s released MATLAB code (github.com/gushiapi/Network-Controllability-Diagnostics, commit
cfa9687): getTieredControlStatesSH.m (per subject, given the initial partition), boundaryDetection.m and Gramian.m.
`gu_reference` below transcribes those three files line by line, with MATLAB's semantics (unique() sorts, min()
returns the first minimum, svd() orders singular values descending, find() returns indices in order), and without
the implementation's caching. nctpy must agree with it exactly. In continuous time the reference swaps in the
continuous-time Gramian, which is the only difference between the two time systems.

All networks are small, seeded and synthetic.
"""

import unittest

import numpy as np
from scipy.linalg import schur

from nctpy.metrics import boundary_control


def gu_reference(A, partition_init, threshold, system="discrete"):
    """Line-by-line transcription of Gu et al.'s getTieredControlStatesSH.m for one network."""

    def gramian(W, B):  # Gramian.m: returns the smallest eigenvalue (smeig)
        A_ = W / (1 + np.max(np.abs(np.linalg.eigvals(W))))  # A = W / (1 + eigs(W, 1))
        if system == "continuous":
            A_ = A_ - np.eye(A_.shape[0])
        T, U = schur(A_, "real")  # [U, T] = schur(A, 'real')
        midMat = U.T @ (B @ B.T) @ U  # midMat = U' * (B * B') * U
        v = np.diag(T)  # v = diag(T)
        if system == "continuous":
            P = midMat / (-(v[:, None] + v[None, :]))
        else:
            P = midMat / (1 - np.outer(v, v))  # P = midMat ./ (1 - v * v')
        eigValues = np.linalg.svd(P, compute_uv=False)  # eigValues = svd(P)
        return eigValues[-1]  # smeig = eigValues(end)

    def boundaryDetection(network, subIndex, thresholdRatio):  # boundaryDetection.m
        subNetwork = network[np.ix_(subIndex, subIndex)]
        thr = np.max(subNetwork) * thresholdRatio
        if subNetwork.shape[0] == 1:
            return subIndex, np.array([], dtype=int)
        L = np.diag(subNetwork.sum(axis=0)) - subNetwork
        V = np.linalg.svd(L)[0]  # [V, ~] = svd(L): MATLAB's first output, the left singular vectors
        fV = V[:, -2]  # fV = V(:, end-1)
        ind1 = np.flatnonzero(fV >= 0)
        ind2 = np.flatnonzero(fV < 0)
        if ind1.size == 0 or ind2.size == 0:
            midNum = (1 + fV.size) // 2  # floor((1 + length(fV)) / 2)
            ind1 = np.arange(midNum)
            ind2 = np.arange(midNum, fV.size)
        newPart = np.zeros(fV.size)
        for j in range(ind1.size):
            if subNetwork[ind1[j], ind2].sum() >= thr:
                newPart[ind1[j]] = 1
        for j in range(ind2.size):
            if subNetwork[ind2[j], ind1].sum() >= thr:
                newPart[ind2[j]] = 1
        return subIndex[newPart == 1], subIndex[ind1]

    A = np.asarray(A, dtype=float)
    N = A.shape[0]
    Kset = np.zeros(N, dtype=int)
    I = np.eye(N)  # noqa: E741 (MATLAB's names, kept for a literal transcription)
    cardK = N
    partition = np.unique(partition_init, return_inverse=True)[1].reshape(-1) + 1  # MATLAB labels 1..k
    communities = np.unique(partition)
    A = A / np.max(A)
    for i in communities:
        ind1 = np.flatnonzero(partition == i)
        ind2 = partition != i
        for j in range(ind1.size):
            if A[ind1[j], ind2].sum() >= threshold:
                Kset[ind1[j]] = 1
    while np.sum(Kset > 0) < cardK:
        uniPVal = np.unique(partition)
        eigVals = np.zeros(uniPVal.size)
        for i in range(uniPVal.size):
            if np.all(Kset[partition == uniPVal[i]]):
                eigVals[i] = np.inf
            else:
                eigVals[i] = gramian(A, I[:, partition == uniPVal[i]])
        if np.all(eigVals == np.inf):
            break
        l = int(np.argmin(eigVals))  # noqa: E741  [~, l] = min(eigVals): the first minimum
        indL = np.flatnonzero(partition == uniPVal[l])
        boundInd, newPart = boundaryDetection(A, indL, threshold)
        if boundInd.size > 0:
            Kset[Kset > 0] = Kset[Kset > 0] + 1
            Kset[boundInd] = 1
        if newPart.size > 0:
            partition[newPart] = np.max(uniPVal) + 1
    ind = Kset.copy()
    for i in range(1, np.max(Kset) + 1):
        ind[Kset == i] = np.sum(Kset <= i)
    return ind


def modular_network(n_per_module, seed):
    """Weighted symmetric network with dense modules and sparse links between them."""
    rng = np.random.default_rng(seed)
    sizes = list(n_per_module)
    labels = np.repeat(np.arange(len(sizes)), sizes)
    n = labels.size
    same = np.equal.outer(labels, labels)
    p = np.where(same, 0.8, 0.08)
    A = np.triu((rng.random((n, n)) < p) * rng.lognormal(0, 1, (n, n)), 1)
    return A + A.T, labels


class TestAgainstGuReference(unittest.TestCase):
    """Exact agreement with the transcription of Gu et al.'s code, on many networks and partitions."""

    def test_matches_reference(self):
        cases = 0
        for seed in range(8):
            for sizes in ((6, 7, 5), (10, 10), (5, 8, 6, 9)):
                A, labels = modular_network(sizes, seed)
                for rho in (0.1, 0.2, 0.3):
                    for system in ("discrete", "continuous"):
                        with self.subTest(seed=seed, sizes=sizes, rho=rho, system=system):
                            expected = gu_reference(A, labels, rho, system) / A.shape[0]
                            np.testing.assert_array_equal(boundary_control(A, labels, system, rho=rho), expected)
                            cases += 1
        self.assertEqual(cases, 144)

    def test_matches_reference_with_a_random_partition(self):
        # a partition unrelated to the network's modules, so that splits and ties are exercised differently
        A, _ = modular_network((8, 8, 8), seed=3)
        labels = np.random.default_rng(5).integers(0, 4, A.shape[0])
        for system in ("discrete", "continuous"):
            with self.subTest(system=system):
                expected = gu_reference(A, labels, 0.2, system) / A.shape[0]
                np.testing.assert_array_equal(boundary_control(A, labels, system), expected)


class TestProperties(unittest.TestCase):
    def test_values(self):
        A, labels = modular_network((7, 9, 6), seed=2)
        for system in ("discrete", "continuous"):
            with self.subTest(system=system):
                bc = boundary_control(A, labels, system)
                self.assertTrue(np.all((bc > 0) & (bc <= 1)))
                self.assertEqual(np.max(bc), 1.0)
                # each tier's value is the fraction of nodes found up to and including that tier
                for value in np.unique(bc):
                    self.assertAlmostEqual(value, np.mean(bc <= value))

    def test_scale_invariance(self):
        # A is scaled to a largest weight of 1, so the result does not depend on its scale. The scale here is a power
        # of two, which rescales exactly: an inexact one moves weights by rounding, and with the Gramian's smallest
        # eigenvalues at rounding level that can change which community is selected (see the docstring).
        A, labels = modular_network((6, 6, 6), seed=4)
        for system in ("discrete", "continuous"):
            with self.subTest(system=system):
                np.testing.assert_array_equal(
                    boundary_control(32.0 * A, labels, system), boundary_control(A, labels, system)
                )

    def test_partition_labels(self):
        # any labels in the same sorted order give the same result
        A, labels = modular_network((6, 7, 5), seed=6)
        names = np.array(["Cont", "Default", "Vis"])[labels]
        np.testing.assert_array_equal(boundary_control(A, names, "discrete"), boundary_control(A, labels, "discrete"))
        np.testing.assert_array_equal(
            boundary_control(A, labels + 10, "discrete"), boundary_control(A, labels, "discrete")
        )


class TestInputs(unittest.TestCase):
    def test_system_is_required(self):
        A, labels = modular_network((5, 5), seed=0)
        for system in (None, "cont"):
            with self.subTest(system=system), self.assertRaises(Exception):
                boundary_control(A, labels, system)

    def test_partition_must_have_one_entry_per_node(self):
        A, labels = modular_network((5, 5), seed=0)
        with self.assertRaises(ValueError):
            boundary_control(A, labels[:-1], "discrete")


if __name__ == "__main__":
    unittest.main()
