"""Reproduction of Kim et al., Nat Commun 16:11639 (2025) through nctpy.optimize: local only.

Refits the published models from nct_xr's local data with the published settings and compares them with nct_xr's
saved results. Nothing is read from, or stored in, this repository: the data and results stay in a local copy of
nct_xr, found through the NCT_XR_DIR environment variable (default: a sibling folder named nct_xr). The test skips
wherever that copy or PyTorch is absent, which includes CI, because real data may never be committed or released.

- Mouse connectome (43 nodes, directed, 36 transitions; CPU, minutes): each refit must stop at the published step,
  with the published decay parameters, and give the published uniform, fitted and static-decay energies.
- HCP-YA (400 nodes; 7 states, so 49 transitions including the 7 self-transitions; best on a GPU, about 20 minutes
  there): opt in with NCTPY_REPRODUCE_HCPYA=<number of transitions to refit>, e.g. 3 or 49.
"""

import contextlib
import io
import os
import unittest
from pathlib import Path

import numpy as np

from nctpy.energies import get_control_inputs, integrate_u
from nctpy.utils import matrix_normalization, normalize_state

REPO = Path(__file__).resolve().parents[2]
NCT_XR = Path(os.environ.get("NCT_XR_DIR", REPO.parent / "nct_xr"))
SETTINGS = "c-1_T-1_rho-1_refstate-xf_initweights-one_nsteps-1000_lr-0.01_eigweight-1.0_regweight-0.0001_regtype-l2"
MOUSE = {
    "A": NCT_XR / "data" / "mouse_cortex_Am.npy",
    "states": NCT_XR / "data" / "mouse_cortex_brain_states.npy",
    "results": NCT_XR / "results" / "mouse" / f"mouse-cortex-Am_optimal-optimized-energy_k-6_{SETTINGS}.npy",
}
HCPYA = {
    "A": NCT_XR / "data" / "HCPYA_Schaefer4007_Am.npy",
    "states": NCT_XR / "results" / "HCPYA" / "HCPYA_Schaefer4007_rsts_fmri_clusters_k-7.npy",
    "results": NCT_XR / "results" / "HCPYA" / f"HCPYA-Schaefer4007-Am_optimal-optimized-energy_k-7_{SETTINGS}.npy",
}


def setUpModule():
    try:
        import torch  # noqa: F401
    except ImportError as exc:
        raise unittest.SkipTest("PyTorch not installed: {0}".format(exc))
    absent = [str(p) for p in MOUSE.values() if not p.exists()]
    if absent:
        raise unittest.SkipTest("nct_xr's local data or results are not available: " + ", ".join(absent))


def load(paths):
    A = np.load(paths["A"])
    states = np.load(paths["states"], allow_pickle=True)
    centroids = states.item()["centroids"] if states.dtype == object else states
    return A, np.nan_to_num(centroids), np.load(paths["results"], allow_pickle=True).item()


def energy(A_norm, x0, xf):
    _, u, _ = get_control_inputs(
        A_norm, 1, np.eye(len(A_norm)), x0, xf, system="continuous", rho=1, S=np.eye(len(A_norm)), xr="xf"
    )
    return np.sum(integrate_u(u))


class Reproduction:
    paths: dict
    device = "cpu"

    def check(self, transitions):
        from nctpy.optimize import optimize_decay_rates

        A, centroids, published = load(self.paths)
        A_norm = matrix_normalization(A, system="continuous", c=1)
        for i, j in transitions:
            with self.subTest(transition=f"{i}->{j}"):
                x0, xf = normalize_state(centroids[i]), normalize_state(centroids[j])
                with contextlib.redirect_stderr(io.StringIO()):
                    fit = optimize_decay_rates(A, x0, xf, device=self.device, progress=False)

                n_steps = int(np.sum(~np.isnan(published["loss"][i, j])))
                self.assertEqual(fit.n_steps_run, n_steps)
                np.testing.assert_allclose(fit.w, published["optimized_weights"][i, j, n_steps - 1], rtol=0, atol=1e-10)
                np.testing.assert_allclose(fit.loss, published["loss"][i, j, :n_steps], rtol=1e-10)

                expected = {
                    "uniform": (A_norm, published["control_energy"][i, j]),
                    "fitted": (fit.A_norm, published["control_energy_variable_decay"][i, j]),
                    "static": (
                        matrix_normalization(A, system="continuous", c=1, decay=np.mean(fit.decay)),
                        published["control_energy_static_decay"][i, j],
                    ),
                }
                for model, (matrix, value) in expected.items():
                    np.testing.assert_allclose(energy(matrix, x0, xf), value, rtol=1e-6, err_msg=model)


class TestMouse(Reproduction, unittest.TestCase):
    paths = MOUSE

    def test_all_transitions(self):
        self.check([(i, j) for i in range(6) for j in range(6)])


class TestHCPYA(Reproduction, unittest.TestCase):
    paths = HCPYA

    def setUp(self):
        self.n_transitions = int(os.environ.get("NCTPY_REPRODUCE_HCPYA", "0"))
        if not self.n_transitions:
            self.skipTest("opt in with NCTPY_REPRODUCE_HCPYA=<number of transitions>")
        if not all(p.exists() for p in self.paths.values()):
            self.skipTest("nct_xr's local HCP-YA data or results are not available")
        import torch

        self.device = "cuda" if torch.cuda.is_available() else "cpu"

    def test_transitions(self):
        self.check([(i, j) for i in range(7) for j in range(7)][: self.n_transitions])


if __name__ == "__main__":
    unittest.main()
