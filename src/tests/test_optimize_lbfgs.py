"""nctpy.optimize, method="lbfgs": the decay-rate loss in closed form, minimised with L-BFGS.

For a symmetric connectome with B and S multiples of the identity, the loss of Kim et al. (2025) and its gradient
have closed forms. These tests check the closed form against the loss computed from get_control_inputs (numpy, no
PyTorch) and the gradient against finite differences; with PyTorch installed, also against the PyTorch model and
its automatic differentiation, and the L-BFGS fit against the Adam fit. They run without PyTorch, except where
marked.
"""

import os
import subprocess
import sys
import unittest
from unittest import mock

import numpy as np
from scipy.optimize import minimize

from nctpy.energies import get_control_inputs, integrate_u
from nctpy.optimize import _closed_form_loss, optimize_decay_rates
from nctpy.utils import matrix_normalization, normalize_state

N = 12

try:
    import torch  # noqa: F401

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False


def connectome(directed=False):
    rng = np.random.default_rng(0)
    A = rng.random((N, N))
    if not directed:
        A = (A + A.T) / 2
    np.fill_diagonal(A, 0)
    return A


def states(k=1, seed=1):
    rng = np.random.default_rng(seed)
    X0 = np.column_stack([normalize_state(rng.random(N)) for _ in range(k)])
    XF = np.column_stack([normalize_state(rng.random(N)) for _ in range(k)])
    return X0, XF


def numpy_loss(A_norm, w, X0, XF, XR, T=1.0, rho=1.0, b=1.0, s=1.0, eig_weight=1.0, reg_weight=1e-4):
    """Kim et al.'s loss computed with get_control_inputs: the reference the closed form must match."""
    A_w = A_norm - np.diag(w)
    distances = []
    for j in range(X0.shape[1]):
        x, _, _ = get_control_inputs(A_w, T, np.sqrt(b) * np.eye(N), X0[:, j], XF[:, j], system="continuous",
                                     rho=rho, S=s * np.eye(N), xr=XR[:, j])  # fmt: skip
        distances.append(np.linalg.norm(x[x.shape[0] // 2] - XR[:, j]))  # x(T/2): T is a multiple of 2 steps here
    return reg_weight * np.sum(np.diag(A_w) ** 2) + np.mean(-eig_weight / np.linalg.eigvalsh(A_w)) + np.mean(distances)


def closed_loss(A_norm, w, X0, XF, XR, T=1.0, rho=1.0, b=1.0, s=1.0, eig_weight=1.0, reg_weight=1e-4, reg_type="l2"):
    return _closed_form_loss(w, A_norm, X0, XF, XR, T, b, s, rho, eig_weight, reg_weight, reg_type)


class TestClosedForm(unittest.TestCase):
    def setUp(self):
        self.A_norm = matrix_normalization(connectome(), system="continuous", zero_diagonal=True)
        self.rng = np.random.default_rng(2)

    def test_loss_matches_get_control_inputs(self):
        for case, (k, T, rho, b, s) in enumerate(((1, 1.0, 1.0, 1.0, 1.0), (3, 1.0, 1.0, 1.0, 1.0),
                                                  (1, 2.0, 0.5, 4.0, 2.0))):  # fmt: skip
            X0, XF = states(k, seed=case)
            for XR in (XF, X0, np.zeros_like(X0), (X0 + XF) / 2):
                w = self.rng.uniform(0.5, 3, N)
                with self.subTest(case=case):
                    got = closed_loss(self.A_norm, w, X0, XF, XR, T=T, rho=rho, b=b, s=s)[0]
                    expected = numpy_loss(self.A_norm, w, X0, XF, XR, T=T, rho=rho, b=b, s=s)
                    np.testing.assert_allclose(got, expected, rtol=1e-10)

    def test_gradient_matches_finite_differences(self):
        X0, XF = states(2)
        for reg_type in ("l2", "l1"):
            w = self.rng.uniform(0.5, 3, N)
            loss, grad, _ = closed_loss(self.A_norm, w, X0, XF, XF, reg_type=reg_type)
            step = 1e-6
            numeric = np.array([
                (closed_loss(self.A_norm, w + step * e, X0, XF, XF, reg_type=reg_type)[0]
                 - closed_loss(self.A_norm, w - step * e, X0, XF, XF, reg_type=reg_type)[0]) / (2 * step)
                for e in np.eye(N)
            ])  # fmt: skip
            with self.subTest(reg_type=reg_type):
                np.testing.assert_allclose(grad, numeric, rtol=1e-6, atol=1e-9)

    def test_largest_eigenvalue(self):
        w = self.rng.uniform(0.5, 3, N)
        X0, XF = states()
        max_eig = closed_loss(self.A_norm, w, X0, XF, XF)[2]
        self.assertAlmostEqual(max_eig, np.max(np.linalg.eigvalsh(self.A_norm - np.diag(w))), places=12)


@unittest.skipUnless(HAS_TORCH, "PyTorch not installed")
class TestAgainstPyTorch(unittest.TestCase):
    def test_loss_and_gradient_match_autograd(self):
        import torch

        from nctpy.optimize import _DecayRateModel

        A_norm = matrix_normalization(connectome(), system="continuous", zero_diagonal=True)
        X0, XF = states(3)
        w = np.random.default_rng(3).uniform(0.5, 3, N)
        model = _DecayRateModel(A_norm, X0, XF, XF, 1, np.eye(N), np.eye(N), 1, "one", 1.0, 1e-4, "l2", "symmetric")
        with torch.no_grad():
            model.w.copy_(torch.from_numpy(w[:, None]))
        loss = model()
        loss.backward()
        got_loss, got_grad, _ = closed_loss(A_norm, w, X0, XF, XF)
        np.testing.assert_allclose(got_loss, loss.item(), rtol=1e-12)
        np.testing.assert_allclose(got_grad, model.w.grad.numpy()[:, 0], rtol=1e-9, atol=1e-13)

    def test_lbfgs_reaches_the_minimum_adam_converges_to(self):
        # Adam run to convergence (no early stopping) and L-BFGS minimise the same loss to the same point. With its
        # default early stopping, Adam can stop short of the minimum, so it is not the comparison here.
        A = connectome()
        x0, xf = (X[:, 0] for X in states(seed=0))
        adam = optimize_decay_rates(A, x0, xf, device="cpu", progress=False, early_stopping=False, n_steps=4000)
        lbfgs = optimize_decay_rates(A, x0, xf, method="lbfgs")
        energies = []
        for fit in (adam, lbfgs):
            _, u, _ = get_control_inputs(fit.A_norm, 1, np.eye(N), x0, xf, "continuous", rho=1, S=np.eye(N), xr="xf")
            energies.append(np.sum(integrate_u(u)))
        self.assertGreaterEqual(np.corrcoef(adam.decay, lbfgs.decay)[0, 1], 0.9999)
        np.testing.assert_allclose(lbfgs.loss[-1], adam.loss[-1], rtol=1e-6)
        np.testing.assert_allclose(energies[1], energies[0], rtol=1e-3)


class TestFit(unittest.TestCase):
    def test_result_is_consistent(self):
        A = connectome()
        x0, xf = (X[:, 0] for X in states())
        fit = optimize_decay_rates(A, x0, xf, method="lbfgs")
        np.testing.assert_array_equal(fit.decay, 1 + fit.w)
        A_norm = matrix_normalization(A, system="continuous", zero_diagonal=True, decay=fit.decay)
        np.testing.assert_allclose(fit.A_norm, A_norm, rtol=0, atol=1e-14)
        self.assertEqual(fit.loss.shape, (fit.n_steps_run,))
        self.assertEqual(fit.max_eigenvalue.shape, (fit.n_steps_run,))
        self.assertLess(fit.max_eigenvalue[-1], 0)
        self.assertTrue(fit.stopped_early)

    def test_restarts_after_stopping_short(self):
        # L-BFGS-B can stop before the minimum, e.g. when a step lands in the unstable region (seen on real
        # connectomes). The fit then restarts from where it stopped. Here the first run is cut short artificially.
        A = connectome()
        x0, xf = (X[:, 0] for X in states())
        runs = []

        def first_run_short(*args, **kwargs):
            if not runs:
                kwargs["options"] = {**kwargs["options"], "maxiter": 2}
            runs.append(1)
            return minimize(*args, **kwargs)

        with mock.patch("nctpy.optimize.minimize", side_effect=first_run_short):
            restarted = optimize_decay_rates(A, x0, xf, method="lbfgs")
        direct = optimize_decay_rates(A, x0, xf, method="lbfgs")
        self.assertGreater(len(runs), 1)
        self.assertEqual(restarted.loss.shape, (restarted.n_steps_run,))
        np.testing.assert_allclose(restarted.loss[-1], direct.loss[-1], rtol=1e-10)
        np.testing.assert_allclose(restarted.decay, direct.decay, rtol=1e-4)

    def test_improves_on_the_start(self):
        A = connectome()
        X0, XF = states(2)
        A_norm = matrix_normalization(A, system="continuous", zero_diagonal=True)
        start = closed_loss(A_norm, np.ones(N), X0, XF, XF)[0]
        fit = optimize_decay_rates(A, X0, XF, method="lbfgs")
        self.assertLess(fit.loss[-1], start)
        self.assertEqual(fit.decay.shape, (N,))

    def test_options(self):
        A = connectome()
        X0, XF = states(2)
        for kwargs in ({"xr": "zero"}, {"xr": "midpoint"}, {"init": "zero"}, {"reg_type": "l1"}, {"rho": 2.0, "T": 2},
                       {"B": 2 * np.eye(N), "S": 0.5 * np.eye(N)}):  # fmt: skip
            with self.subTest(**{k: str(v)[:12] for k, v in kwargs.items()}):
                fit = optimize_decay_rates(A, X0, XF, method="lbfgs", **kwargs)
                self.assertTrue(np.all(np.isfinite(fit.decay)))
                self.assertLess(fit.max_eigenvalue[-1], 0)

    def test_invalid_input_raises(self):
        x0, xf = (X[:, 0] for X in states())
        with self.assertRaises(ValueError):
            optimize_decay_rates(connectome(directed=True), x0, xf, method="lbfgs")
        partial = np.eye(N)
        partial[0, 0] = 0
        for kwargs in ({"B": partial}, {"S": partial}, {"method": "newton"}):
            with self.subTest(**{k: str(v)[:12] for k, v in kwargs.items()}), self.assertRaises(ValueError):
                optimize_decay_rates(connectome(), x0, xf, **{"method": "lbfgs", **kwargs})


class TestWithoutPyTorch(unittest.TestCase):
    def test_imports_and_lbfgs_runs_without_torch(self):
        # a fresh interpreter in which torch cannot be imported
        script = (
            "import sys\n"
            "class Block:\n"
            "    def find_spec(self, name, path=None, target=None):\n"
            "        if name.split('.')[0] == 'torch':\n"
            "            raise ModuleNotFoundError('blocked: ' + name)\n"
            "sys.meta_path.insert(0, Block())\n"
            "import numpy as np\n"
            "from nctpy.optimize import optimize_decay_rates\n"
            "rng = np.random.default_rng(0)\n"
            "A = rng.random((8, 8)); A = (A + A.T) / 2; np.fill_diagonal(A, 0)\n"
            "fit = optimize_decay_rates(A, rng.random(8), rng.random(8), method='lbfgs')\n"
            "print('lbfgs ran:', fit.decay.shape)\n"
            "try:\n"
            "    optimize_decay_rates(A, rng.random(8), rng.random(8))\n"
            "except ImportError as exc:\n"
            "    print('adam raised ImportError:', exc)\n"
        )
        import nctpy

        env = dict(os.environ, PYTHONPATH=os.pathsep.join(
            filter(None, [os.path.dirname(os.path.dirname(nctpy.__file__)), os.environ.get("PYTHONPATH")])))  # fmt: skip
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=env)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("lbfgs ran: (8,)", result.stdout)
        self.assertIn("pip install 'nctpy[optimize]'", result.stdout)


if __name__ == "__main__":
    unittest.main()
