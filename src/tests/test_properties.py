"""Property tests: mathematical facts the package's outputs must satisfy (Roadmap 0.6).

Fixtures catch drift from what nctpy used to return; these catch errors that were always there.
Each property below holds by theory. Where theory needs a condition, the test states it.

- Minimum-energy statements (energy rises as the control set shrinks, falls as a control weight
  grows, is independent of rho) are guaranteed for minimum-energy control, S = 0. With S != 0 the
  solver trades trajectory cost against input cost, so input energy alone is not guaranteed to be
  monotone; those cases are not asserted here.
- Discrete ave_control (the Schur-based formula) equals the diagonal of the infinite-horizon
  discrete Gramian exactly only for normal (e.g. symmetric) A; for directed A it approximates it,
  so that property is checked on symmetric A only.
- sim_state_eq integrates continuous time with forward Euler steps of 0.001, so it agrees with
  exact solutions to O(dt); discrete-time simulation is exact.

All systems are small, seeded and synthetic.
"""
import unittest

import numpy as np
import scipy.linalg as la

from nctpy.energies import get_control_inputs, integrate_u, gramian, minimum_energy_fast, sim_state_eq
from nctpy.metrics import ave_control, modal_control
from nctpy.utils import matrix_normalization, normalize_state, normalize_weights

N = 12
DT = 0.001  # continuous-time step used by get_control_inputs and sim_state_eq
THRESHOLD = 1e-8


def connectomes():
    rng = np.random.default_rng(0)
    A = rng.random((N, N)) + 0.1
    symmetric = (A + A.T) / 2
    np.fill_diagonal(symmetric, 0)
    directed = rng.random((N, N)) + 0.1
    np.fill_diagonal(directed, 0)
    return {'symmetric': symmetric, 'directed': directed}


X0 = normalize_state(np.r_[np.ones(3), np.zeros(N - 3)])
XF = normalize_state(np.r_[np.zeros(N - 3), np.ones(3)])


def solve(A_norm, B=None, S=None, system='continuous', T=1, rho=1, x0=X0, xf=XF):
    B = np.eye(N) if B is None else B
    S = np.eye(N) if S is None else S
    x, u, err = get_control_inputs(A_norm, T, B, x0, xf, system=system, rho=rho, S=S)
    return integrate_u(u), x, u, err


class PropertyTestCase(unittest.TestCase):
    def setUp(self):
        self.systems = {}
        for name, A in connectomes().items():
            self.systems[name] = (A, matrix_normalization(A, system='continuous'),
                                  matrix_normalization(A, system='discrete'))

    def assertCompletes(self, err):
        self.assertLess(err[0], THRESHOLD)
        self.assertLess(err[1], THRESHOLD)


class TestNormalization(PropertyTestCase):
    def test_stability(self):
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                self.assertLess(np.max(np.abs(np.linalg.eigvals(A_d))), 1)  # discrete: spectral radius < 1
                self.assertLess(np.max(np.linalg.eigvals(A_c).real), 0)  # continuous: all real parts < 0

    def test_eigenvalues_map_affinely(self):
        # A_norm = A / (|lambda|max + c), minus I in continuous time: every eigenvalue is mapped by the same
        # increasing affine map, so their order is preserved. Real and imaginary parts are compared as sorted
        # multisets, since complex-conjugate pairs share a real part and so have no unique order.
        for name, (A, A_c, A_d) in self.systems.items():
            ev = np.linalg.eigvals(A)
            scale = np.max(np.abs(ev)) + 1  # c = 1
            for system, A_norm, shift in (('continuous', A_c, -1), ('discrete', A_d, 0)):
                with self.subTest(connectome=name, system=system):
                    ev_norm = np.linalg.eigvals(A_norm)
                    np.testing.assert_allclose(np.sort(ev_norm.real), np.sort(ev.real) / scale + shift, atol=1e-12)
                    np.testing.assert_allclose(np.sort(ev_norm.imag), np.sort(ev.imag) / scale, atol=1e-12)


class TestMinimumEnergy(PropertyTestCase):
    """S = 0: get_control_inputs solves the minimum-energy problem."""

    def test_matches_minimum_energy_fast(self):
        # integrate_u integrates with unit sample spacing, so energy per unit time is its sum times dt
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                energy, _, _, err = solve(A_c, S=np.zeros((N, N)))
                self.assertCompletes(err)
                expected = np.sum(minimum_energy_fast(A_c, 1, np.eye(N), X0, XF))
                self.assertAlmostEqual(np.sum(energy) * DT / expected, 1.0, places=8)

    def test_independent_of_rho(self):
        for name, (A, A_c, A_d) in self.systems.items():
            for system, A_norm, T in (('continuous', A_c, 1), ('discrete', A_d, 5)):
                with self.subTest(connectome=name, system=system):
                    energies = [np.sum(solve(A_norm, S=np.zeros((N, N)), system=system, T=T, rho=rho)[0])
                                for rho in (0.01, 1, 100)]
                    np.testing.assert_allclose(energies, energies[1], rtol=1e-10)

    def test_energy_rises_as_control_set_shrinks(self):
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                energies = []
                for k in range(4):
                    B = np.eye(N)
                    B[np.arange(k), np.arange(k)] = 0  # remove control from nodes 0..k-1
                    energy, _, _, err = solve(A_c, B=B, S=np.zeros((N, N)))
                    self.assertCompletes(err)
                    energies.append(np.sum(energy))
                self.assertTrue(np.all(np.diff(energies) > 0), energies)

    def test_energy_falls_as_a_control_weight_grows(self):
        # the SI's premise for the gradient step on B (section C2, 3d)
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                base = np.sum(solve(A_c, S=np.zeros((N, N)))[0])
                for i in range(N):
                    B = np.eye(N)
                    B[i, i] += 0.1
                    self.assertLess(np.sum(solve(A_c, B=B, S=np.zeros((N, N)))[0]), base)


class TestRelabelling(PropertyTestCase):
    """Permuting the nodes permutes every node-level output and leaves totals unchanged."""

    def test_control_outputs_permute(self):
        p = np.random.default_rng(1).permutation(N)
        P = np.eye(N)[p]
        for name, (A, A_c, A_d) in self.systems.items():
            for system, A_norm, T in (('continuous', A_c, 1), ('discrete', A_d, 5)):
                with self.subTest(connectome=name, system=system):
                    energy, x, _, _ = solve(A_norm, system=system, T=T)
                    energy_p, x_p, _, _ = solve(P @ A_norm @ P.T, system=system, T=T, x0=X0[p], xf=XF[p])
                    np.testing.assert_allclose(energy_p, energy[p], rtol=1e-10)
                    np.testing.assert_allclose(x_p, x[:, p], atol=1e-12)
                    self.assertAlmostEqual(np.sum(energy_p) / np.sum(energy), 1.0, places=12)

    def test_metrics_permute(self):
        p = np.random.default_rng(2).permutation(N)
        P = np.eye(N)[p]
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                np.testing.assert_allclose(ave_control(P @ A_c @ P.T, 'continuous'), ave_control(A_c, 'continuous')[p],
                                           rtol=1e-10)
                np.testing.assert_allclose(matrix_normalization(P @ A @ P.T, system='continuous'), P @ A_c @ P.T,
                                           atol=1e-14)

    def test_schur_metrics_permute_for_symmetric_A(self):
        p = np.random.default_rng(2).permutation(N)
        P = np.eye(N)[p]
        A, A_c, A_d = self.systems['symmetric']
        np.testing.assert_allclose(ave_control(P @ A_d @ P.T, 'discrete'), ave_control(A_d, 'discrete')[p], rtol=1e-10)
        np.testing.assert_allclose(modal_control(P @ A_d @ P.T), modal_control(A_d)[p], rtol=1e-10)

    @unittest.expectedFailure
    def test_schur_metrics_permute_for_directed_A(self):
        # Known limitation (Roadmap D15): discrete ave_control and modal_control read only the diagonal of a
        # real Schur decomposition. That is exact for normal (e.g. symmetric) A, but for directed A the result
        # depends on node order: relabelling changes them by ~5e-4 relative here, while the exact quantity (the
        # infinite-horizon Gramian diagonal) does not change. Remove the decorator if the formulas are changed.
        p = np.random.default_rng(2).permutation(N)
        P = np.eye(N)[p]
        A, A_c, A_d = self.systems['directed']
        np.testing.assert_allclose(ave_control(P @ A_d @ P.T, 'discrete'), ave_control(A_d, 'discrete')[p], rtol=1e-10)
        np.testing.assert_allclose(modal_control(P @ A_d @ P.T), modal_control(A_d)[p], rtol=1e-10)


class TestAverageControllability(PropertyTestCase):
    def test_continuous_is_gramian_diagonal(self):
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                np.testing.assert_allclose(ave_control(A_c, 'continuous'), np.diag(gramian(A_c, 1, 'continuous')),
                                           rtol=1e-12)

    def test_discrete_is_infinite_gramian_diagonal_for_symmetric_A(self):
        A, A_c, A_d = self.systems['symmetric']
        np.testing.assert_allclose(ave_control(A_d, 'discrete'), np.diag(gramian(A_d, np.inf, 'discrete')), rtol=1e-10)


class TestSimulation(PropertyTestCase):
    """The optimal-control solver and the forward simulator describe the same system."""

    def test_discrete_solver_trajectory_is_reproduced(self):
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                _, x, u, _ = solve(A_d, system='discrete', T=5)
                # x has T + 1 states and u has T inputs; sim_state_eq records the state before each input
                x_sim = sim_state_eq(A_d, np.eye(N), X0, np.c_[u.T, np.zeros((N, 1))], system='discrete')
                np.testing.assert_allclose(x_sim.T, x, atol=1e-12)

    def test_continuous_solver_trajectory_is_reproduced_to_first_order(self):
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                _, x, u, _ = solve(A_c)
                x_sim = sim_state_eq(A_c, np.eye(N), X0, u.T, system='continuous')
                np.testing.assert_allclose(x_sim.T, x, atol=1e-3)  # forward Euler, dt = 0.001

    def test_uncontrolled_continuous_simulation_follows_the_matrix_exponential(self):
        for name, (A, A_c, A_d) in self.systems.items():
            with self.subTest(connectome=name):
                x_sim = sim_state_eq(A_c, np.eye(N), X0, np.zeros((N, 1001)), system='continuous')
                np.testing.assert_allclose(x_sim[:, -1], la.expm(A_c * 1.0) @ X0, atol=1e-3)


class TestSigns(PropertyTestCase):
    def test_energy_is_non_negative(self):
        rng = np.random.default_rng(3)
        weighted = np.diag(normalize_weights(rng.normal(size=N)))
        partial_S = np.diag(np.r_[np.ones(N // 2), np.zeros(N - N // 2)])
        for name, (A, A_c, A_d) in self.systems.items():
            for system, A_norm, T in (('continuous', A_c, 1), ('discrete', A_d, 5)):
                for B, S in ((np.eye(N), np.eye(N)), (weighted, partial_S), (np.eye(N), np.zeros((N, N)))):
                    with self.subTest(connectome=name, system=system):
                        self.assertTrue(np.all(solve(A_norm, B=B, S=S, system=system, T=T)[0] >= 0))


class TestUtils(unittest.TestCase):
    def test_normalize_state_has_unit_norm(self):
        rng = np.random.default_rng(4)
        for x in (rng.normal(size=N), rng.random(N) > 0.5, np.arange(1, N + 1)):
            self.assertAlmostEqual(np.linalg.norm(normalize_state(x)), 1.0, places=12)

    def test_normalize_weights_ranges(self):
        x = np.random.default_rng(5).normal(size=50)
        for rank in (True, False):
            with self.subTest(rank=rank):
                w = normalize_weights(x, rank=rank)
                self.assertAlmostEqual(np.min(w), 1.0)
                self.assertAlmostEqual(np.max(w), 2.0)
                w = normalize_weights(x, rank=rank, add_constant=False)
                self.assertAlmostEqual(np.min(w), 0.0)
                self.assertAlmostEqual(np.max(w), 1.0)
        # ranking preserves order
        np.testing.assert_array_equal(np.argsort(normalize_weights(x)), np.argsort(x))


if __name__ == '__main__':
    unittest.main()
