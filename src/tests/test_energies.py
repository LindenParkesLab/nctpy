"""nctpy.energies: input handling and the Gramian's branches (Roadmap 2.2).

States may be Boolean, integer or float, as 1-D vectors or (N, 1) columns: all are converted to float64 and
give the same results. Array inputs are computed in float64 whatever their precision.
"""

import contextlib
import io
import unittest

import numpy as np
import scipy.linalg as la

from nctpy.energies import get_control_inputs, gramian, minimum_energy_fast, sim_state_eq
from nctpy.utils import matrix_normalization

N = 8
HORIZON = {"continuous": 1, "discrete": 4}


def connectome():
    rng = np.random.default_rng(0)
    A = rng.random((N, N))
    A = (A + A.T) / 2
    np.fill_diagonal(A, 0)
    return A


def state_variants(mask):
    """The same state as every accepted input form; the first is the reference (1-D float64)."""
    return {
        "float 1-D": mask.astype(float),
        "bool 1-D": mask,
        "bool (N, 1)": mask[:, None],
        "int 1-D": mask.astype(int),
        "float32 1-D": mask.astype(np.float32),
        "float (N, 1)": mask.astype(float)[:, None],
    }


class TestStateInputs(unittest.TestCase):
    def setUp(self):
        self.A = connectome()
        self.A_norm = {s: matrix_normalization(self.A, system=s) for s in HORIZON}
        self.x0 = state_variants(np.arange(N) < 3)
        self.xf = state_variants(np.arange(N) >= N - 3)

    def test_get_control_inputs(self):
        for system, T in HORIZON.items():
            for xr in ("zero", "midpoint"):
                args = (self.A_norm[system], T, np.eye(N))
                x_ref, u_ref, err_ref = get_control_inputs(
                    *args, self.x0["float 1-D"], self.xf["float 1-D"], system=system, xr=xr
                )
                for form in self.x0:
                    with self.subTest(system=system, xr=xr, form=form):
                        x, u, err = get_control_inputs(*args, self.x0[form], self.xf[form], system=system, xr=xr)
                        np.testing.assert_array_equal(x, x_ref)
                        np.testing.assert_array_equal(u, u_ref)
                        self.assertEqual(err, err_ref)
                        self.assertEqual((x.dtype, u.dtype), (np.float64, np.float64))

    def test_minimum_energy_fast(self):
        A_c = self.A_norm["continuous"]
        ref = minimum_energy_fast(A_c, 1, np.eye(N), self.x0["float 1-D"], self.xf["float 1-D"])
        for form in self.x0:
            with self.subTest(form=form):
                np.testing.assert_array_equal(minimum_energy_fast(A_c, 1, np.eye(N), self.x0[form], self.xf[form]), ref)

    def test_sim_state_eq(self):
        U = np.random.default_rng(1).normal(size=(N, 20))
        for system in HORIZON:
            ref = sim_state_eq(self.A_norm[system], np.eye(N), self.x0["float 1-D"], U, system=system)
            for form in self.x0:
                with self.subTest(system=system, form=form):
                    x = sim_state_eq(self.A_norm[system], np.eye(N), self.x0[form], U, system=system)
                    np.testing.assert_array_equal(x, ref)


class TestShapes(unittest.TestCase):
    """The trajectory and input lengths the docstring states: time x nodes."""

    def test_lengths(self):
        A = connectome()
        x0, xf = np.arange(N) < 3, np.arange(N) >= N - 3
        for system, T, n_x, n_u in (("continuous", 1, 1001, 1001), ("discrete", 4, 5, 4), ("discrete", 2, 3, 2)):
            with self.subTest(system=system, T=T):
                x, u, err = get_control_inputs(
                    matrix_normalization(A, system=system), T, np.eye(N), x0, xf, system=system
                )
                self.assertEqual((x.shape, u.shape, len(err)), ((n_x, N), (n_u, N), 2))

    def test_discrete_horizon_below_two_raises(self):
        A_d = matrix_normalization(connectome(), system="discrete")
        for T in (0, 1):
            with self.subTest(T=T):
                with self.assertRaises(Exception) as ctx:
                    get_control_inputs(A_d, T, np.eye(N), np.arange(N) < 3, np.arange(N) >= 5, system="discrete")
                self.assertIs(type(ctx.exception), Exception)
                self.assertEqual(str(ctx.exception), "Discrete time systems must have T >= 2")


class TestGramian(unittest.TestCase):
    def setUp(self):
        self.A = connectome()
        self.A_c = matrix_normalization(self.A, system="continuous")
        self.A_d = matrix_normalization(self.A, system="discrete")

    def test_continuous_infinite_horizon_solves_lyapunov(self):
        W = gramian(self.A_c, np.inf, system="continuous")
        residual = self.A_c @ W + W @ self.A_c.T + np.eye(N)
        np.testing.assert_allclose(residual, 0, atol=1e-10 * np.abs(W).max())

    def test_finite_horizon_matches_infinite_minus_tail(self):
        # for stable A, W(T) = W(inf) - e^{AT} W(inf) e^{A^T T}; checks the finite-horizon Simpson integration
        W_inf = gramian(self.A_c, np.inf, system="continuous")
        for T in (0.5, 1, 3):
            with self.subTest(T=T):
                E = la.expm(self.A_c * T)
                expected = W_inf - E @ W_inf @ E.T
                W = gramian(self.A_c, T, system="continuous")
                np.testing.assert_allclose(W, expected, rtol=0, atol=1e-10 * np.abs(expected).max())

    def test_unstable_infinite_horizon_returns_nan(self):
        # the raw connectome is unstable in both senses: positive entries give a positive, > 1 leading eigenvalue
        for system in ("continuous", "discrete"):
            with self.subTest(system=system):
                out = io.StringIO()
                with contextlib.redirect_stdout(out):
                    W = gramian(self.A, np.inf, system=system)
                self.assertTrue(np.isnan(W))
                self.assertEqual(out.getvalue(), "cannot compute infinite-time Gramian for an unstable system!\n")

    def test_lower_precision_input_computed_in_float64(self):
        for system, A_norm in (("continuous", self.A_c), ("discrete", self.A_d)):
            with self.subTest(system=system):
                np.testing.assert_array_equal(
                    gramian(A_norm.astype(np.float32).astype(float), 2, system=system),
                    gramian(A_norm.astype(np.float32), 2, system=system),
                )


class TestReuseAcrossTransitions(unittest.TestCase):
    """System-only work is reused between calls (Roadmap 2.5); results must not depend on call order."""

    def setUp(self):
        rng = np.random.default_rng(3)
        self.A = connectome()
        self.states = [(rng.normal(size=N), rng.normal(size=N)) for _ in range(3)]
        partial = np.diag((np.arange(N) % 3 > 0).astype(float))
        self.systems = [
            (system, T, B, S, rho)
            for system, T in HORIZON.items()
            for B in (np.eye(N), partial)
            for S, rho in ((np.eye(N), 1), (np.zeros((N, N)), 0.5))
        ]

    def solve(self, system, T, B, S, rho, x0, xf, A=None):
        A_norm = matrix_normalization(self.A if A is None else A, system=system)
        return get_control_inputs(A_norm, T, B, x0, xf, system=system, rho=rho, S=S)

    def assertSameResult(self, got, expected):
        np.testing.assert_array_equal(got[0], expected[0])
        np.testing.assert_array_equal(got[1], expected[1])
        self.assertEqual(got[2], expected[2])

    def test_call_order_does_not_matter(self):
        calls = [(s, x) for s in range(len(self.systems)) for x in range(len(self.states))]
        by_system = {c: self.solve(*self.systems[c[0]], *self.states[c[1]]) for c in calls}
        interleaved = calls[::2] + calls[1::2]  # alternate between systems on every call
        for c in reversed(interleaved):
            with self.subTest(system=c[0], states=c[1]):
                self.assertSameResult(self.solve(*self.systems[c[0]], *self.states[c[1]]), by_system[c])

    def test_matrix_changed_in_place_is_recomputed(self):
        for system, T in HORIZON.items():
            with self.subTest(system=system):
                A_norm = matrix_normalization(self.A, system=system)
                x0, xf = self.states[0]
                get_control_inputs(A_norm, T, np.eye(N), x0, xf, system=system)
                A_norm[0, 1] += 0.01  # same array object, new contents
                got = get_control_inputs(A_norm, T, np.eye(N), x0, xf, system=system)
                expected = get_control_inputs(A_norm.copy(), T, np.eye(N) * 1, x0, xf, system=system)
                self.assertSameResult(got, expected)

    def test_minimum_energy_fast_batches_states(self):
        # the numerics page's advice: one Gramian for many transitions, passed as columns
        A_c = matrix_normalization(self.A, system="continuous")
        X0 = np.column_stack([x0 for x0, _ in self.states])
        XF = np.column_stack([xf for _, xf in self.states])
        batched = minimum_energy_fast(A_c, 1, np.eye(N), X0, XF)
        for k, (x0, xf) in enumerate(self.states):
            with self.subTest(transition=k):
                # equal to rounding: BLAS multiplies a block of columns and a single column differently
                single = minimum_energy_fast(A_c, 1, np.eye(N), x0, xf)
                np.testing.assert_allclose(batched[:, [k]], single, rtol=1e-12, atol=1e-14 * np.abs(single).max())


if __name__ == "__main__":
    unittest.main()
