"""nctpy.utils."""

import unittest

import numpy as np

from nctpy.utils import (
    convert_states_str2int,
    expand_states,
    expm,
    get_fdr_p,
    get_null_p,
    get_p_val_string,
    mask_control_set,
    matrix_normalization,
    normalize_weights,
    random_control_set,
)


def connectome(n=10):
    rng = np.random.default_rng(0)
    A = rng.random((n, n))
    return (A + A.T) / 2


def spectral_radius(A):
    return np.abs(np.linalg.eigvals(A)).max()


class TestMatrixNormalization(unittest.TestCase):
    def test_default_l_is_own_spectral_radius(self):
        A = connectome()
        for system in ("continuous", "discrete"):
            with self.subTest(system=system):
                np.testing.assert_array_equal(
                    matrix_normalization(A, system=system), matrix_normalization(A, system=system, l=spectral_radius(A))
                )

    def test_discrete_spectral_radius(self):
        # A / (c + l) scales the spectral radius to rho / (c + l): below 1 when c + l > rho
        A = connectome()
        rho = spectral_radius(A)
        for c, fixed in ((1, None), (0.5, None), (1, rho / 2), (0, 2 * rho)):
            with self.subTest(c=c, l=fixed):
                expected = rho / (c + (rho if fixed is None else fixed))
                got = spectral_radius(matrix_normalization(A, system="discrete", c=c, l=fixed))
                self.assertAlmostEqual(got, expected)

    def test_shared_l_scales_identically(self):
        # one l across connectomes: each is divided by the same constant
        A1, A2 = connectome(), 2 * connectome()
        shared = max(spectral_radius(A1), spectral_radius(A2))
        N1 = matrix_normalization(A1, system="discrete", l=shared)
        N2 = matrix_normalization(A2, system="discrete", l=shared)
        np.testing.assert_allclose(N2, 2 * N1, rtol=1e-14)

    def test_small_l_can_be_unstable_and_still_returns(self):
        A = connectome()
        A_norm = matrix_normalization(A, system="continuous", c=0, l=spectral_radius(A) / 2)
        self.assertGreater(np.linalg.eigvals(A_norm).real.max(), 0)

    def test_lower_precision_input_computed_in_float64(self):
        A32 = connectome().astype(np.float32)
        for system in ("continuous", "discrete"):
            with self.subTest(system=system):
                got = matrix_normalization(A32, system=system)
                self.assertEqual(got.dtype, np.float64)
                np.testing.assert_array_equal(got, matrix_normalization(A32.astype(np.float64), system=system))


class TestZeroDiagonal(unittest.TestCase):
    """matrix_normalization(..., zero_diagonal=True) applies the no-self-connections assumption."""

    def setUp(self):
        self.A = connectome()
        np.fill_diagonal(self.A, np.arange(1, 11))  # self-connections, as in the protocol paper's PNC connectome

    def test_default_keeps_the_diagonal(self):
        for system in ("continuous", "discrete"):
            with self.subTest(system=system):
                np.testing.assert_array_equal(
                    matrix_normalization(self.A, system=system, zero_diagonal=False),
                    matrix_normalization(self.A, system=system),
                )
                self.assertFalse(
                    np.array_equal(
                        matrix_normalization(self.A, system=system),
                        matrix_normalization(self.A, system=system, zero_diagonal=True),
                    )
                )

    def test_equals_normalising_a_zero_diagonal_copy(self):
        A0 = self.A.copy()
        np.fill_diagonal(A0, 0)
        for system in ("continuous", "discrete"):
            for c, fixed in ((1, None), (0.5, None), (1, 7.0)):
                with self.subTest(system=system, c=c, l=fixed):
                    np.testing.assert_array_equal(
                        matrix_normalization(self.A, system=system, c=c, l=fixed, zero_diagonal=True),
                        matrix_normalization(A0, system=system, c=c, l=fixed),
                    )

    def test_input_not_modified(self):
        for A in (self.A, self.A.astype(np.float32), self.A.round().astype(int)):
            with self.subTest(dtype=A.dtype):
                before = A.copy()
                matrix_normalization(A, system="continuous", zero_diagonal=True)
                np.testing.assert_array_equal(A, before)

    def test_keyword_only(self):
        with self.assertRaises(TypeError):
            matrix_normalization(self.A, "continuous", 1, None, True)


class TestDecay(unittest.TestCase):
    """matrix_normalization(..., decay=) implements Kim et al. (2025), Eq. 4."""

    def setUp(self):
        rng = np.random.default_rng(5)
        n = 10
        undirected = connectome(n)
        directed = rng.random((n, n))
        for A in (undirected, directed):
            np.fill_diagonal(A, rng.uniform(1, 5, size=n))  # self-connections, as in the PNC connectome
        self.connectomes = {"undirected": undirected, "directed": directed}
        self.decay = rng.uniform(0.5, 2.5, size=n)

    def test_default_is_unit_decay(self):
        for name, A in self.connectomes.items():
            with self.subTest(connectome=name):
                default = matrix_normalization(A, system="continuous")
                np.testing.assert_array_equal(matrix_normalization(A, system="continuous", decay=1), default)
                np.testing.assert_array_equal(matrix_normalization(A, system="continuous", decay=np.ones(10)), default)

    def test_decay_is_subtracted_and_self_connections_are_kept(self):
        for name, A in self.connectomes.items():
            for c, fixed in ((1, None), (0.5, 9.0)):
                with self.subTest(connectome=name, c=c, l=fixed):
                    scale = c + (spectral_radius(A) if fixed is None else fixed)
                    A_norm = matrix_normalization(A, system="continuous", c=c, l=fixed, decay=self.decay)
                    off = ~np.eye(10, dtype=bool)
                    np.testing.assert_array_equal(A_norm[off], (A / scale)[off])
                    np.testing.assert_allclose(np.diag(A_norm), np.diag(A) / scale - self.decay, rtol=1e-14)

    def test_composes_with_zero_diagonal(self):
        for name, A in self.connectomes.items():
            with self.subTest(connectome=name):
                A_norm = matrix_normalization(A, system="continuous", zero_diagonal=True, decay=self.decay)
                np.testing.assert_array_equal(np.diag(A_norm), -self.decay)
                A0 = A.copy()
                np.fill_diagonal(A0, 0)
                off = ~np.eye(10, dtype=bool)
                np.testing.assert_array_equal(A_norm[off], matrix_normalization(A0, system="continuous")[off])

    def test_uniform_decay_of_at_least_one_is_stable(self):
        for name, A in self.connectomes.items():
            for decay in (1, 1.5, 4):
                with self.subTest(connectome=name, decay=decay):
                    eigvals = np.linalg.eigvals(matrix_normalization(A, system="continuous", decay=decay))
                    self.assertLess(eigvals.real.max(), 0)

    def test_invalid_input_raises(self):
        A = self.connectomes["undirected"]
        with self.assertRaisesRegex(ValueError, "continuous-time systems only"):
            matrix_normalization(A, system="discrete", decay=self.decay)
        for bad in (np.ones(9), np.ones((10, 1)), np.ones((10, 10))):
            with self.subTest(shape=bad.shape), self.assertRaisesRegex(ValueError, "one value per node"):
                matrix_normalization(A, system="continuous", decay=bad)
        matrix_normalization(A, system="discrete")  # decay=None is fine in discrete time

    def test_keyword_only(self):
        with self.assertRaises(TypeError):
            matrix_normalization(self.connectomes["undirected"], "continuous", 1, None, False, self.decay)


class TestControlSets(unittest.TestCase):
    def test_random_control_set(self):
        for n, k in ((43, 5), (400, 124)):
            for seed in (0, 7):
                with self.subTest(n=n, k=k, seed=seed):
                    B = random_control_set(n, k, seed=seed)
                    self.assertEqual(B.shape, (n, n))
                    np.testing.assert_array_equal(B, np.diag(np.diag(B)))
                    self.assertEqual(np.count_nonzero(B), k)
                    np.testing.assert_array_equal(random_control_set(n, k, seed=seed), B)
        self.assertFalse(np.array_equal(random_control_set(400, 124, seed=0), random_control_set(400, 124, seed=1)))

    def test_random_control_set_matches_seeding_numpy_globally(self):
        # the same seed draws the same control nodes as np.random.seed(seed); np.random.choice(...), which is how
        # the control sets of Kim et al. (2025) were drawn, without touching numpy's global state
        saved = np.random.get_state()
        try:
            for seed in range(5):
                np.random.seed(seed)
                expected = np.sort(np.random.choice(np.arange(400), size=64, replace=False))
                np.random.seed(12345)
                before = np.random.get_state()[1].copy()
                got = np.flatnonzero(np.diag(random_control_set(400, 64, seed=seed)))
                np.testing.assert_array_equal(got, expected)
                np.testing.assert_array_equal(np.random.get_state()[1], before)
        finally:
            np.random.set_state(saved)

    def test_baseline(self):
        B = random_control_set(50, 10, seed=3, baseline=1e-5)
        weights = np.diag(B)
        self.assertEqual(np.sum(weights == 1), 10)
        self.assertTrue(np.all(weights[weights != 1] == 1e-5))

    def test_mask_control_set(self):
        mask = np.array([True, False, True, False])
        np.testing.assert_array_equal(mask_control_set(mask), np.diag([1.0, 0.0, 1.0, 0.0]))
        np.testing.assert_array_equal(mask_control_set(mask.astype(int)), mask_control_set(mask))
        np.testing.assert_array_equal(mask_control_set(mask, baseline=1e-3), np.diag([1.0, 1e-3, 1.0, 1e-3]))
        states, labels = convert_states_str2int(["Vis", "Vis", "Default", "Default"])
        np.testing.assert_array_equal(mask_control_set(states == labels.index("Vis")), np.diag([1.0, 1.0, 0.0, 0.0]))

    def test_invalid_input_raises(self):
        with self.assertRaisesRegex(ValueError, "one-dimensional"):
            mask_control_set(np.ones((4, 4), dtype=bool))
        with self.assertRaises(ValueError):
            random_control_set(10, 11)


class TestStates(unittest.TestCase):
    def test_expand_states(self):
        states = np.array([0, 0, 1, 1, 2, 2])
        x0_mat, xf_mat = expand_states(states)
        self.assertEqual((x0_mat.shape, xf_mat.shape), ((6, 9), (6, 9)))
        self.assertEqual((x0_mat.dtype, xf_mat.dtype), (bool, bool))
        for i in range(3):
            for j in range(3):
                with self.subTest(i=i, j=j):
                    np.testing.assert_array_equal(x0_mat[:, i * 3 + j], states == i)
                    np.testing.assert_array_equal(xf_mat[:, i * 3 + j], states == j)

    def test_convert_states_str2int(self):
        names = ["Vis", "SomMot", "Vis", "Default", "SomMot"]
        states, labels = convert_states_str2int(names)
        self.assertIsInstance(labels, list)
        self.assertEqual(labels, ["Default", "SomMot", "Vis"])
        np.testing.assert_array_equal(states, [2, 1, 2, 0, 1])
        self.assertEqual(states.dtype, int)
        self.assertEqual([labels[s] for s in states], names)


class TestStatistics(unittest.TestCase):
    def test_get_null_p_versions(self):
        null = np.array([1.0, 2.0, 3.0, 4.0])
        self.assertEqual(get_null_p(3.0, null), 0.5)  # null >= 3
        self.assertEqual(get_null_p(3.0, null, version="reverse"), 0.75)  # null <= 3
        self.assertEqual(get_null_p(3.0, null, version="smallest"), 0.5)
        self.assertEqual(get_null_p(-3.0, null, abs=True), 0.5)

    def test_get_null_p_unknown_version_raises(self):
        with self.assertRaises(ValueError):
            get_null_p(1.0, np.arange(5), version="two-sided")

    def test_get_fdr_p_any_shape(self):
        p = np.random.default_rng(1).random((3, 4, 5))
        flat = get_fdr_p(p.ravel())
        for shape in ((60,), (6, 10), (3, 4, 5)):
            with self.subTest(shape=shape):
                got = get_fdr_p(p.reshape(shape))
                self.assertEqual(got.shape, shape)
                np.testing.assert_array_equal(got.ravel(), flat)

    def test_get_p_val_string(self):
        self.assertEqual(get_p_val_string(0.0), r"-log10($\mathit{p}$)>25")
        self.assertEqual(get_p_val_string(0.0012), r"$\mathit{p}$ = 1e-03")
        self.assertEqual(get_p_val_string(0.25), r"$\mathit{p}$ = 0.250")

    def test_normalize_weights_default_range(self):
        w = normalize_weights(np.random.default_rng(2).normal(size=20))
        self.assertEqual((w.min(), w.max()), (1.0, 2.0))


class TestExpm(unittest.TestCase):
    def test_matches_scipy(self):
        import scipy.linalg as la

        A = matrix_normalization(connectome(), system="continuous")
        np.testing.assert_allclose(expm(A), la.expm(A), rtol=1e-10, atol=1e-14)


if __name__ == "__main__":
    unittest.main()
