import unittest
import sys
import numpy as np
from network_control.utils import matrix_normalization, normalize_state, normalize_weights
from network_control.metrics import ave_control
from network_control.energies import integrate_u, get_control_inputs


class TestMatrixNormalization(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/A.npy', 'rb') as f:
            self.A = np.load(f)
        self.eps = 1e-10

    def test_matrix_normalization_stability(self):
        A = np.random.randn(20, 20)
        A = (A + A.T) / 2
        # discrete
        norm = matrix_normalization(A, system='discrete')
        w, _ = np.linalg.eig(norm)
        l = np.max(np.abs(w))
        self.assertLess(l, 1)
        # continuous
        norm = matrix_normalization(A, system='continuous')
        w, _ = np.linalg.eig(norm)
        self.assertTrue((l > 0).all())

    def test_matrix_normalization_success(self):
        # discrete time system default c=1
        with open('./fixtures/A_d_1.npy', 'rb') as f:
            A_d_1 = np.load(f)
        result = matrix_normalization(self.A, system='discrete')
        self.assertTrue((np.abs(A_d_1 - result) <= self.eps).all())
        # discrete time system c=2
        with open('./fixtures/A_d_2.npy', 'rb') as f:
            A_d_2 = np.load(f)
        result = matrix_normalization(self.A, system='discrete', c=2)
        self.assertTrue((np.abs(A_d_2 - result) <= self.eps).all())
        self.assertTrue((matrix_normalization(self.A, system='discrete', c=1) != result).any())
        # continuous time system c=1
        with open('./fixtures/A_c_1.npy', 'rb') as f:
            A_c_1 = np.load(f)
        result = matrix_normalization(self.A, system='continuous', c=1)
        self.assertTrue((np.abs(A_c_1 - result) <= self.eps).all())
        self.assertTrue((matrix_normalization(self.A, system='discrete', c=1) != result).any())
        # continuous time system c=2
        with open('./fixtures/A_c_2.npy', 'rb') as f:
            A_c_2 = np.load(f)
        result = matrix_normalization(self.A, system='continuous', c=2)
        self.assertTrue((np.abs(A_c_2 - result) <= self.eps).all())
        self.assertTrue((matrix_normalization(self.A, system='continuous', c=1) != result).all())

    def test_system_specification_error(self):
        # no default
        with self.assertRaises(Exception) as exception_context:
            matrix_normalization(self.A)
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                         "(see function help).")

        # typo
        with self.assertRaises(Exception) as exception_context:
            matrix_normalization(self.A, system="discrt")
        self.assertEqual(str(exception_context.exception),
                         "Incorrect system specification. "
                         "Please specify either 'system=discrete' or 'system=continuous'.")


class TestNormalizeState(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/x.npy', 'rb') as f:
            self.x = np.load(f)
        self.eps = 1e-10

    def test_normalize_state_success(self):
        with open('./fixtures/x_norm.npy', 'rb') as f:
            x_norm = np.load(f)
        self.assertTrue(((x_norm - normalize_state(self.x)) <= self.eps).all())


class TestNormalizeWeights(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/x.npy', 'rb') as f:
            self.weights = np.load(f)
        self.eps = 1e-10

    def test_normalize_weights_success(self):
        # defaults
        with open('./fixtures/weights_rank_scale_const.npy', 'rb') as f:
            weight_r_s_c = np.load(f)
        self.assertTrue((np.abs(weight_r_s_c - normalize_weights(self.weights)) <= self.eps).all())
        self.assertTrue((np.abs(weight_r_s_c - normalize_weights(self.weights,
                                                                 rank=True,
                                                                 add_constant=True)) <= self.eps).all())
        # no rank
        with open('./fixtures/weights_scale_const.npy', 'rb') as f:
            weight_s_c = np.load(f)
        self.assertTrue((np.abs(weight_s_c - normalize_weights(self.weights,
                                                               rank=False) <= self.eps)).all())
        self.assertTrue((normalize_weights(self.weights) != normalize_weights(self.weights,
                                                                              rank=False)).any())
        # no constant
        with open('./fixtures/weights_rank_scale.npy', 'rb') as f:
            weight_r_s = np.load(f)
        self.assertTrue((np.abs(weight_r_s - normalize_weights(self.weights,
                                                               add_constant=False)) <= self.eps).all())
        self.assertTrue((normalize_weights(self.weights) != normalize_weights(self.weights,
                                                                              add_constant=False)).any())
        # no rank or constant
        with open('./fixtures/weights_scale.npy', 'rb') as f:
            weight_s = np.load(f)
        self.assertTrue((np.abs(weight_s - normalize_weights(self.weights,
                                                             rank=False,
                                                             add_constant=False)) <= self.eps).all())
        self.assertTrue((normalize_weights(self.weights) != normalize_weights(self.weights,
                                                                              rank=False,
                                                                              add_constant=False)).any())


class TestGetControlInputs(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/A_d_1.npy', 'rb') as f:
            self.A_d = np.load(f)
        with open('./fixtures/A_c_1.npy', 'rb') as f:
            self.A_c = np.load(f)
        with open('./fixtures/x.npy', 'rb') as f:
            self.x0 = np.load(f)
        with open('./fixtures/B.npy', 'rb') as f:
            self.B = np.load(f)
        with open('./fixtures/xf.npy', 'rb') as f:
            self.xf = np.load(f)
        self.n = np.shape(self.A_d)[0]
        self.eps = 1e-10

    def test_get_control_inputs_success(self):
        # discrete
        with open('./fixtures/control_discrete.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, self.xf, system='discrete')
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())
        # T
        with open('./fixtures/control_T.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_d, 7, np.eye(self.n), self.x0, self.xf, system='discrete')
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())
        # B
        with open('./fixtures/control_B.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_c, 2, self.B, self.x0, self.xf, system='continuous')
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())
        # rho
        with open('./fixtures/control_rho.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n),
                                                      self.x0, self.xf, system='discrete',
                                                      rho=100)
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())
        # S
        with open('./fixtures/control_S.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0,
                                                      self.xf, system='discrete', S=self.B)
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())
        # system
        with open('./fixtures/control_continuous.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_c, 2, np.eye(self.n), self.x0, self.xf,
                                                      system='continuous')
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())
        # reference state
        with open('./fixtures/control_ref.npz', 'rb') as f:
            data = np.load(f)
            x = data['x']
            u = data['u']
            err = data['err']
        x_test, u_test, err_test = get_control_inputs(self.A_c, 2, np.eye(self.n), self.x0, self.xf,
                                                      system='continuous', xr='x0')
        self.assertTrue((np.abs(x - x_test) <= self.eps).all())
        self.assertTrue((np.abs(u - u_test) <= self.eps).all())
        self.assertTrue((np.abs(err - err_test) <= self.eps).all())

    def test_get_control_inputs_reaching_xf(self):
        for i in range(10):
            # for states with increasingly large values, are we always getting to the final?
            xf = np.random.rand(self.n, ) * ((i + 1) * 10)
            x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, xf, system='discrete')
            for j, state in enumerate(x_test[-1, :]):
                self.assertAlmostEqual(state, xf[j], places=3)

    def test_get_control_inputs_consistency(self):
        # boolean states
        with open('./fixtures/x_bin.npy', 'rb') as f:
            x_bin = np.load(f)
        x_bool = x_bin == 1
        x_bi, u_bi, err_bi = get_control_inputs(self.A_d, 2, np.eye(self.n), x_bin, x_bin, system='discrete')
        x_bo, u_bo, err_bo = get_control_inputs(self.A_d, 2, np.eye(self.n),
                                                x_bool, x_bool, system='discrete')
        self.assertTrue((x_bi == x_bo).all())
        self.assertTrue((u_bi == u_bo).all())
        self.assertEqual(err_bi, err_bo)
        # state dimensions
        x0 = np.random.rand(self.n, )
        xf = np.random.rand(self.n, )
        x_1d, u_1d, err_1d = get_control_inputs(self.A_d, 2, np.eye(self.n), x0, xf, system='discrete')
        x_2d, u_2d, err_2d = get_control_inputs(self.A_d, 2, np.eye(self.n),
                                                x0.reshape(-1, 1), xf.reshape(-1, 1), system='discrete')
        self.assertTrue((x_1d == x_2d).all())
        self.assertTrue((u_1d == u_2d).all())
        self.assertEqual(err_1d, err_2d)
        # specifying reference state
        x_d, u_d, err_d = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, self.xf,
                                             system='discrete', xr='x0')
        x_s, u_s, err_s = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, self.xf,
                                             system='discrete', xr=self.x0)
        self.assertTrue((x_d == x_s).all())
        self.assertTrue((u_d == u_s).all())
        self.assertEqual(err_d, err_s)
        x_d, u_d, err_d = get_control_inputs(self.A_c, 2, np.eye(self.n), self.x0, self.xf,
                                             system='continuous', xr='xf')
        x_s, u_s, err_s = get_control_inputs(self.A_c, 2, np.eye(self.n), self.x0, self.xf,
                                             system='continuous', xr=self.xf)
        self.assertTrue((x_d == x_s).all())
        self.assertTrue((u_d == u_s).all())
        self.assertEqual(err_d, err_s)

    def test_get_control_inputs_bounds(self):
        # AC id a lower bound on input
        with open('./fixtures/ac_d.npy', 'rb') as f:
            ac = np.load(f)
        for i in range(10):
            xf = np.random.rand(self.n, )
            xf = xf / np.linalg.norm(xf, ord=2)
            _, u_1d, _ = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, xf, system='discrete')
            self.assertTrue((u_1d <= ac).all())

    def test_get_control_inputs_error(self):
        # no system
        with self.assertRaises(Exception) as exception_context:
            get_control_inputs(self.A_d, 5, np.eye(self.n), self.x0, self.xf)
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                         "(see matrix_normalization help).")

        # typo
        with self.assertRaises(Exception) as exception_context:
            get_control_inputs(self.A_c, 1, self.B, self.x0, self.xf, system='cont')
        self.assertEqual(str(exception_context.exception),
                         "Incorrect system specification. "
                         "Please specify either 'system=discrete' or 'system=continuous'.")


class TestIntegrateU(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/u.npy', 'rb') as f:
            self.u = np.load(f)
        self.eps = 2 * np.finfo(float).eps

    def test_integrate_u_bounds(self):
        for i in range(10):
            # for increasingly large U, test that we always get positive energy
            u = np.random.randn(100, 10000) * ((i + 1) * 10)
            self.assertTrue((integrate_u(u) > 0).all())

    def test_integrate_u_success(self):
        with open('./fixtures/u_int.npy', 'rb') as f:
            energy = np.load(f)
        self.assertTrue((np.abs(energy - integrate_u(self.u)) <= self.eps).all())
        # TODO different scipy version (how do I do this???)


class TestAveControl(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/A_d_1.npy', 'rb') as f:
            self.A_d = np.load(f)
        with open('./fixtures/A_c_1.npy', 'rb') as f:
            self.A_c = np.load(f)
        self.eps = 1e-13

    def test_ave_control_success(self):
        # discrete
        with open('./fixtures/ac_d.npy', 'rb') as f:
            ac = np.load(f)
        self.assertTrue((np.abs(ac - ave_control(self.A_d, 'discrete')) <= self.eps).all())
        # continuous
        with open('./fixtures/ac_c.npy', 'rb') as f:
            ac = np.load(f)
        self.assertTrue((np.abs(ac - ave_control(self.A_c, 'continuous')) <= self.eps).all())
        self.assertTrue((ave_control(self.A_d, 'discrete') != ave_control(self.A_c, 'continuous')).any())

    def test_ave_control_error(self):
        # no system
        with self.assertRaises(Exception) as exception_context:
            ave_control(self.A_d)
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                         "(see matrix_normalization help).")

        # typo
        with self.assertRaises(Exception) as exception_context:
            ave_control(self.A_c, 'contigiuous')
        self.assertEqual(str(exception_context.exception),
                         "Incorrect system specification. "
                         "Please specify either 'system=discrete' or 'system=continuous'.")


if __name__ == '__main__':
    # setting path
    sys.path.append('../network_control/')

    unittest.main()
