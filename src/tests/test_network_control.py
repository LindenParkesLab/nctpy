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
        self.assertTrue((A_d_1 == matrix_normalization(self.A, system='discrete')).all())
        # discrete time system c=2
        with open('./fixtures/A_d_2.npy', 'rb') as f:
            A_d_2 = np.load(f)
        self.assertTrue((A_d_2 == matrix_normalization(self.A, system='discrete', c=2)).all())
        self.assertTrue((matrix_normalization(self.A, system='discrete', c=1) != matrix_normalization(self.A,
                                                                                                      system='discrete',
                                                                                                      c=2)).any())
        # continuous time system c=1
        with open('./fixtures/A_c_1.npy', 'rb') as f:
            A_c_1 = np.load(f)
        self.assertTrue((A_c_1 == matrix_normalization(self.A, system='continuous', c=1)).all())
        self.assertTrue((matrix_normalization(self.A, system='discrete', c=1) != matrix_normalization(self.A,
                                                                                                      system='continuous',
                                                                                                      c=1)).any())
        # continuous time system c=2
        with open('./fixtures/A_c_2.npy', 'rb') as f:
            A_c_2 = np.load(f)
        self.assertTrue((A_c_2 == matrix_normalization(self.A, system='continuous', c=2)).all())
        self.assertTrue((matrix_normalization(self.A, system='continuous', c=1) != matrix_normalization(self.A,
                                                                                                        system='continuous',
                                                                                                        c=2)).any())

    def test_system_specification_error(self):
        # no default
        with self.assertRaises(Exception) as exception_context:
            matrix_normalization(self.A)
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please nominate whether you are normalizing A for a continuous-time or a discrete-time "
                         "system "
                         "(see function help).")

        # typo
        with self.assertRaises(Exception) as exception_context:
            matrix_normalization(self.A, system="discrt")
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please nominate whether you are normalizing A for a continuous-time or a discrete-time "
                         "system "
                         "(see function help).")

    # TODO add test for dimensions of A


class TestNormalizeState(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/x.npy', 'rb') as f:
            self.x = np.load(f)

    def test_normalize_state_success(self):
        with open('./fixtures/x_norm.npy', 'rb') as f:
            x_norm = np.load(f)
        self.assertTrue((x_norm == normalize_state(self.x)).all())

        # TODO add type and size checking for inputs


class TestNormalizeWeights(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/x.npy', 'rb') as f:
            self.weights = np.load(f)

    def test_normalize_weights_success(self):
        # defaults
        with open('./fixtures/weights_rank_scale_const.npy', 'rb') as f:
            weight_r_s_c = np.load(f)
        self.assertTrue((weight_r_s_c == normalize_weights(self.weights)).all())
        self.assertTrue((weight_r_s_c == normalize_weights(self.weights,
                                                           rank=True,
                                                           add_constant=True)).all())
        # no rank
        with open('./fixtures/weights_scale_const.npy', 'rb') as f:
            weight_s_c = np.load(f)
        self.assertTrue((weight_s_c == normalize_weights(self.weights,
                                                         rank=False)).all())
        self.assertTrue((normalize_weights(self.weights) != normalize_weights(self.weights,
                                                                              rank=False)).any())
        # no constant
        with open('./fixtures/weights_rank_scale.npy', 'rb') as f:
            weight_r_s = np.load(f)
        self.assertTrue((weight_r_s == normalize_weights(self.weights,
                                                         add_constant=False)).all())
        self.assertTrue((normalize_weights(self.weights) != normalize_weights(self.weights,
                                                                              add_constant=False)).any())
        # no rank or constant
        with open('./fixtures/weights_scale.npy', 'rb') as f:
            weight_s = np.load(f)
        self.assertTrue((weight_s == normalize_weights(self.weights,
                                                       rank=False,
                                                       add_constant=False)).all())
        self.assertTrue((normalize_weights(self.weights) != normalize_weights(self.weights,
                                                                              rank=False,
                                                                              add_constant=False)).any())

        # TODO add type and size checking for inputs


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

    def test_get_control_inputs_success(self):
        # discrete
        with open('./fixtures/control_discrete.npy', 'rb') as f:
            x, u, err = np.load(f)
        x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, self.xf, system='discrete')
        self.assertTrue((x == x_test).all())
        self.assertTrue((u == u_test).all())
        self.assertTrue((err == err_test).all())
        # T
        with open('./fixtures/control_T.npy', 'rb') as f:
            x, u, err = np.load(f)
        x_test, u_test, err_test = get_control_inputs(self.A_d, 7, np.eye(self.n), self.x0, self.xf, system='discrete')
        self.assertTrue((x == x_test).all())
        self.assertTrue((u == u_test).all())
        self.assertTrue((err == err_test).all())
        # TODO test for T=1
        # B
        with open('./fixtures/control_B.npy', 'rb') as f:
            x, u, err = np.load(f)
        x_test, u_test, err_test = get_control_inputs(self.A_c, 2, self.B, self.x0, self.xf, system='continuous')
        self.assertTrue((x == x_test).all())
        self.assertTrue((u == u_test).all())
        self.assertTrue((err == err_test).all())
        # rho
        with open('./fixtures/control_rho.npy', 'rb') as f:
            x, u, err = np.load(f)
        x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n),
                                                      self.x0, self.xf, system='discrete',
                                                      rhp=100)
        self.assertTrue((x == x_test).all())
        self.assertTrue((u == u_test).all())
        self.assertTrue((err == err_test).all())
        # S
        with open('./fixtures/control_T.npy', 'rb') as f:
            x, u, err = np.load(f)
        x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0,
                                                      self.xf, system='discrete', S=self.B)
        self.assertTrue((x == x_test).all())
        self.assertTrue((u == u_test).all())
        self.assertTrue((err == err_test).all())
        # system
        with open('./fixtures/control_continuous.npy', 'rb') as f:
            x, u, err = np.load(f)
        x_test, u_test, err_test = get_control_inputs(self.A_c, 2, np.eye(self.n), self.x0, self.xf,
                                                      system='continuous')
        self.assertTrue((x == x_test).all())
        self.assertTrue((u == u_test).all())
        self.assertTrue((err == err_test).all())
        # TODO test reference state

    def test_get_control_inputs_reaching_xf(self):
        for i in range(10):
            # for states with increasingly large values, are we always getting to the final?
            xf = np.random.rand(self.n, ) * ((i+1)*10)
            x_test, u_test, err_test = get_control_inputs(self.A_d, 2, np.eye(self.n), self.x0, xf, system='discrete')
            for j, state in enumerate(x_test):
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
        self.assertTrue((err_bi == err_bo).all())
        # state dimensions
        x0 = np.random.rand(self.n, )
        xf = np.random.rand(self.n, )
        x_1d, u_1d, err_1d = get_control_inputs(self.A_d, 2, np.eye(self.n), x0, xf, system='discrete')
        x_2d, u_2d, err_2d = get_control_inputs(self.A_d, 2, np.eye(self.n),
                                                x0.reshape(-1, 1), xf.reshape(-1, 1), system='discrete')
        self.assertTrue((x_1d == x_2d).all())
        self.assertTrue((u_1d == u_2d).all())
        self.assertTrue((err_1d == err_2d).all())

    # TODO test that energy required to get to 1st eig is sim to IR
    def test_get_control_inputs_bounds(self):
        self.assertEqual(True, False)

    def test_get_control_inputs_error(self):
        # no system
        with self.assertRaises(Exception) as exception_context:
            get_control_inputs(self.A_d, 5, np.eye(self.n), self.x0, self.xf)
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please indicate whether you are simulating a continuous-time or a discrete-time system "
                         "(see matrix_normalization for help).")
        # typo
        with self.assertRaises(Exception) as exception_context:
            get_control_inputs(self.A_c, 1, self.B, self.x0, self.xf, system='cont')
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please indicate whether you are simulating a continuous-time or a discrete-time system "
                         "(see matrix_normalization for help).")


class TestIntegrateU(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/u.npy', 'rb') as f:
            self.u = np.load(f)

    def test_integrate_u_bounds(self):
        for i in range(10):
            # for increasingly large U, test that we always get positive energy
            u = np.random.randn(100, 10000) * ((i + 1) * 10)
            self.assertTrue((integrate_u(u) > 0).all())

    def test_integrate_u_success(self):
        with open('./fixtures/u_int.npy', 'rb') as f:
            energy = np.load(f)
        self.assertTrue((energy == integrate_u(self.u)).all())
        # TODO different scipy version (how do I do this???)


class TestAveControl(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/A_d_1.npy', 'rb') as f:
            self.A_d = np.load(f)
        with open('./fixtures/A_c_1.npy', 'rb') as f:
            self.A_c = np.load(f)

    def test_ave_control_success(self):
        # discrete
        with open('./fixtures/ac_d.npy', 'rb') as f:
            ac = np.load(f)
        self.assertTrue((ac == ave_control(self.A_d, 'discrete')).all())
        # continuous
        with open('./fixtures/ac_c.npy', 'rb') as f:
            ac = np.load(f)
        self.assertTrue((ac == ave_control(self.A_c, 'continuous')).all())
        self.assertTrue((ave_control(self.A_d, 'discrete') != ave_control(self.A_c, 'continuous')).any())

    def test_ave_control_error(self):
        # no system
        with self.assertRaises(Exception) as exception_context:
            ave_control(self.A_d)
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please indicate whether you are simulating a continuous-time or a discrete-time system "
                         "(see matrix_normalization for help).")

        # typo
        with self.assertRaises(Exception) as exception_context:
            ave_control(self.A_c, 'contigiuous')
        self.assertEqual(str(exception_context.exception),
                         "Time system not specified. "
                         "Please indicate whether you are simulating a continuous-time or a discrete-time system "
                         "(see matrix_normalization for help).")


if __name__ == '__main__':
    # setting path
    sys.path.append('../network_control/')

    unittest.main()
