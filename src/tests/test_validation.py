"""Input validation.

Validation rejects input that does not define a problem, and nothing else: an ill-conditioned or
incomplete control problem must still return its energies and error terms.
"""

import contextlib
import io
import unittest

import numpy as np

from nctpy.energies import get_control_inputs, gramian, sim_state_eq
from nctpy.metrics import ave_control
from nctpy.pipelines import ComputeControlEnergy, ComputeOptimizedControlEnergy
from nctpy.utils import matrix_normalization, normalize_state

N = 8
MISSING = (
    "Time system not specified. "
    "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
    "(see matrix_normalization help)."
)
INVALID = "Incorrect system specification. Please specify either 'system=discrete' or 'system=continuous'."


def system_inputs():
    rng = np.random.default_rng(0)
    A = rng.random((N, N)) + 0.1
    A = (A + A.T) / 2
    np.fill_diagonal(A, 0)
    x0 = normalize_state(np.arange(N) < 3)
    xf = normalize_state(np.arange(N) >= N - 3)
    return A, x0, xf


class TestSystem(unittest.TestCase):
    """The system checks: same exception type (bare Exception) and messages as before."""

    def setUp(self):
        self.A, self.x0, self.xf = system_inputs()
        self.A_c = matrix_normalization(self.A, system="continuous")

    def assertRaisesExactly(self, message, fn, *args, **kwargs):
        with self.assertRaises(Exception) as ctx:
            fn(*args, **kwargs)
        self.assertIs(type(ctx.exception), Exception)
        self.assertEqual(str(ctx.exception), message)

    def test_get_control_inputs(self):
        for system, message in ((None, MISSING), ("cont", INVALID)):
            with self.subTest(system=system):
                self.assertRaisesExactly(
                    message, get_control_inputs, self.A_c, 1, np.eye(N), self.x0, self.xf, system=system
                )

    def test_sim_state_eq(self):
        U = np.zeros((N, 5))
        for system, message in ((None, MISSING), ("cont", INVALID)):
            with self.subTest(system=system):
                self.assertRaisesExactly(message, sim_state_eq, self.A_c, np.eye(N), self.x0, U, system=system)

    def test_gramian(self):
        # until 1.1.0 gramian returned None here
        for T in (1, np.inf):
            for system, message in ((None, MISSING), ("cont", INVALID)):
                with self.subTest(T=T, system=system):
                    self.assertRaisesExactly(message, gramian, self.A_c, T, system=system)

    def test_matrix_normalization(self):
        # utils points to its own docstring rather than to matrix_normalization's
        for system, message in (
            (None, MISSING.replace("matrix_normalization help", "function help")),
            ("cont", INVALID),
        ):
            with self.subTest(system=system):
                self.assertRaisesExactly(message, matrix_normalization, self.A, system=system)

    def test_ave_control(self):
        for system, message in ((None, MISSING), ("cont", INVALID)):
            with self.subTest(system=system):
                self.assertRaisesExactly(message, ave_control, self.A_c, system=system)


class TestRho(unittest.TestCase):
    """rho must be positive; rho <= 0 used to return NaN with numpy warnings."""

    def setUp(self):
        self.A, self.x0, self.xf = system_inputs()
        self.A_norm = {s: matrix_normalization(self.A, system=s) for s in ("continuous", "discrete")}
        self.T = {"continuous": 1, "discrete": 5}

    def solve(self, system, rho, S=None):
        S = np.eye(N) if S is None else S
        return get_control_inputs(
            self.A_norm[system], self.T[system], np.eye(N), self.x0, self.xf, system=system, rho=rho, S=S
        )

    def test_non_positive_rho_raises(self):
        for system in ("continuous", "discrete"):
            for rho in (0, 0.0, -1, np.nan):
                for S in (np.eye(N), np.zeros((N, N))):  # including minimum-energy control
                    with self.subTest(system=system, rho=rho, S_is_zero=not S.any()):
                        with self.assertRaises(ValueError) as ctx:
                            self.solve(system, rho, S)
                        self.assertIn("rho must be positive", str(ctx.exception))

    def test_positive_rho_returns(self):
        # a very small rho gives a poorly solved problem, not an invalid one: it returns
        for system in ("continuous", "discrete"):
            for rho in (1e-6, 1, 100):
                with self.subTest(system=system, rho=rho):
                    x, u, err = self.solve(system, rho)
                    self.assertEqual(x.shape[1], N)
                    self.assertEqual(len(err), 2)

    def test_pipelines_raise_for_rho_zero(self):
        task = {"x0": self.x0, "xf": self.xf, "B": np.eye(N), "S": np.eye(N), "rho": 0}
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(ValueError):
                ComputeControlEnergy(A=self.A, control_tasks=[task], system="continuous").run()
            with self.assertRaises(ValueError):
                ComputeOptimizedControlEnergy(A=self.A, control_task=task, system="continuous", n_steps=1).run()


if __name__ == "__main__":
    unittest.main()
