"""nctpy.optimize: fitting nodes' decay rates (Kim et al., Nat Commun 2025).

The cross-check is the important test: the PyTorch loss, which drives the optimisation, must agree with the same
quantities computed by nctpy's numpy implementation, so that the two implementations of the state-costate solution
cannot drift apart. Runs on CPU, in float64. Skipped where PyTorch is not installed.
"""

import os
import subprocess
import sys
import unittest

import numpy as np

from nctpy.energies import get_control_inputs
from nctpy.utils import matrix_normalization, normalize_state

N = 8


def setUpModule():
    try:
        import torch  # noqa: F401
    except ImportError as exc:
        raise unittest.SkipTest("PyTorch not installed: {0}".format(exc))


def connectome(directed=False, self_connections=False):
    rng = np.random.default_rng(0)
    A = rng.random((N, N))
    if not directed:
        A = (A + A.T) / 2
    np.fill_diagonal(A, rng.uniform(1, 3, size=N) if self_connections else 0)
    return A


def states(k=1, seed=1):
    rng = np.random.default_rng(seed)
    X0 = np.column_stack([normalize_state(rng.random(N)) for _ in range(k)])
    XF = np.column_stack([normalize_state(rng.random(N)) for _ in range(k)])
    return (X0[:, 0], XF[:, 0]) if k == 1 else (X0, XF)


def fit(*args, **kwargs):
    from nctpy.optimize import optimize_decay_rates

    kwargs.setdefault("device", "cpu")
    kwargs.setdefault("progress", False)
    return optimize_decay_rates(*args, **kwargs)


class TestCrossCheck(unittest.TestCase):
    """The PyTorch loss against numpy (get_control_inputs), at fixed decay parameters."""

    def numpy_loss(self, A_norm, w, x0, xf, eig_weight=1.0, reg_weight=1e-4):
        A_sub = A_norm - np.diag(w)
        x, _, _ = get_control_inputs(A_sub, 1, np.eye(N), x0, xf, system="continuous", rho=1, S=np.eye(N), xr="xf")
        x_mid = x[500]  # T/2 with T = 1 and steps of 0.001
        stability = np.mean(-eig_weight / np.linalg.eigvalsh(A_sub))
        return reg_weight * np.sum(np.diag(A_sub) ** 2) + stability + np.linalg.norm(x_mid - xf)

    def test_loss_matches_numpy(self):
        import torch

        from nctpy.optimize import _DecayRateModel

        A_norm = matrix_normalization(connectome(), system="continuous")
        x0, xf = states()
        for w in (np.ones(N), np.random.default_rng(2).uniform(0.5, 2, N)):
            with self.subTest(w=w[:2]):
                model = _DecayRateModel(
                    A_norm, x0[:, None], xf[:, None], xf[:, None], 1, np.eye(N), np.eye(N), 1, "one", 1.0, 1e-4, "l2",
                    "symmetric",
                )  # fmt: skip
                with torch.no_grad():
                    model.w.copy_(torch.from_numpy(w[:, None]))
                    loss = model().item()
                np.testing.assert_allclose(loss, self.numpy_loss(A_norm, w, x0, xf), rtol=1e-10)


class TestFit(unittest.TestCase):
    def test_result_is_consistent(self):
        A = connectome()
        x0, xf = states()
        result = fit(A, x0, xf, n_steps=40)
        np.testing.assert_array_equal(result.decay, 1 + result.w)
        A_norm = matrix_normalization(A, system="continuous", zero_diagonal=True, decay=result.decay)
        np.testing.assert_allclose(result.A_norm, A_norm, rtol=0, atol=1e-14)
        self.assertEqual(result.loss.shape, (result.n_steps_run,))
        self.assertEqual(result.max_eigenvalue.shape, (result.n_steps_run,))
        self.assertTrue(np.all(np.isfinite(result.loss)))
        self.assertLess(result.max_eigenvalue[-1], 0)

    def test_deterministic_on_cpu(self):
        x0, xf = states()
        first, second = fit(connectome(), x0, xf, n_steps=30), fit(connectome(), x0, xf, n_steps=30)
        np.testing.assert_array_equal(first.w, second.w)
        np.testing.assert_array_equal(first.loss, second.loss)

    def test_init(self):
        x0, xf = states()
        for init, start in (("one", 2.0), ("zero", 1.0)):
            with self.subTest(init=init):
                result = fit(connectome(), x0, xf, init=init, n_steps=1, lr=1e-12)
                np.testing.assert_allclose(result.decay, start, atol=1e-9)

    def test_early_stopping(self):
        x0, xf = states()
        result = fit(connectome(), x0, xf, n_steps=3000)
        self.assertTrue(result.stopped_early)
        self.assertLess(result.n_steps_run, 3000)
        result = fit(connectome(), x0, xf, n_steps=60, early_stopping=False)
        self.assertFalse(result.stopped_early)
        self.assertEqual(result.n_steps_run, 60)

    def test_self_connections_are_removed_by_default(self):
        A = connectome(self_connections=True)
        x0, xf = states()
        result = fit(A, x0, xf, n_steps=20)
        np.testing.assert_allclose(np.diag(result.A_norm), -result.decay, rtol=0, atol=1e-15)
        A0 = A.copy()
        np.fill_diagonal(A0, 0)
        off = ~np.eye(N, dtype=bool)
        np.testing.assert_array_equal(result.A_norm[off], matrix_normalization(A0, system="continuous")[off])

    def test_several_transitions_and_reference_states(self):
        X0, XF = states(k=3)
        for xr in ("xf", "midpoint", "zero", XF):
            with self.subTest(xr=xr if isinstance(xr, str) else "array"):
                result = fit(connectome(), X0, XF, xr=xr, n_steps=10)
                self.assertEqual(result.decay.shape, (N,))

    def test_general_eigenvalues(self):
        x0, xf = states()
        # identical in exact arithmetic for a symmetric connectome
        symmetric = fit(connectome(), x0, xf, n_steps=30)
        general = fit(connectome(), x0, xf, n_steps=30, eigenvalues="general")
        np.testing.assert_allclose(general.w, symmetric.w, rtol=1e-8)
        # for a directed connectome, 'general' tracks the true largest real part of the eigenvalues
        result = fit(connectome(directed=True), x0, xf, n_steps=30, eigenvalues="general")
        np.testing.assert_allclose(result.max_eigenvalue[-1], np.max(np.linalg.eigvals(result.A_norm).real), rtol=1e-12)

    def test_invalid_input_raises(self):
        x0, xf = states()
        for kwargs in ({"init": "half"}, {"reg_type": "l3"}, {"eigenvalues": "complex"}, {"rho": 0}, {"xr": "target"}):
            with self.subTest(**kwargs), self.assertRaises(ValueError):
                fit(connectome(), x0, xf, n_steps=1, **kwargs)


class TestImport(unittest.TestCase):
    def test_import_error_names_the_extra(self):
        # a fresh interpreter in which torch cannot be imported
        script = (
            "import sys\n"
            "class Block:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] == 'torch':\n"
            "            raise ModuleNotFoundError('blocked: ' + name)\n"
            "sys.meta_path.insert(0, Block())\n"
            "import nctpy.energies\n"
            "try:\n"
            "    import nctpy.optimize\n"
            "except ImportError as exc:\n"
            "    print('optimize raised ImportError:', exc)\n"
        )
        import nctpy

        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            filter(None, [os.path.dirname(os.path.dirname(nctpy.__file__)), os.environ.get("PYTHONPATH")])))  # fmt: skip
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("pip install 'nctpy[optimize]'", result.stdout)


if __name__ == "__main__":
    unittest.main()
