import unittest
import sys
import numpy as np

# setting path
sys.path.append('../network_control/')
from network_control.utils import matrix_normalization, normalize_state, normalize_weights
from network_control.metrics import ave_control
from network_control.energies import integrate_u, get_control_inputs


class TestMatrixNormalization(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/A.npy', 'rb') as f:
            self.A = np.load(f)

    def test_matrix_normalization_stability(self):
        A = np.random.randn(20, 20)
        A = (A + A.T)/2
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
        A_d_1_test = matrix_normalization(self.A, system='discrete')
        self.assertTrue((A_d_1 == A_d_1_test).all())
        # discrete time system c=2
        with open('./fixtures/A_d_2.npy', 'rb') as f:
            A_d_2 = np.load(f)
        self.assertTrue((A_d_2 == matrix_normalization(self.A, system='discrete', c=2)).all())
        # continuous time system c=1
        with open('./fixtures/A_c_1.npy', 'rb') as f:
            A_c_1 = np.load(f)
        self.assertTrue((A_c_1 == matrix_normalization(self.A, system='continuous', c=1)).all())
        # continuous time system c=2
        with open('./fixtures/A_c_2.npy', 'rb') as f:
            A_c_2 = np.load(f)
        self.assertTrue((A_c_2 == matrix_normalization(self.A, system='continuous', c=2)).all())

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
        # no constant
        with open('./fixtures/weights_rank_scale.npy', 'rb') as f:
            weight_r_s = np.load(f)
        self.assertTrue((weight_r_s == normalize_weights(self.weights,
                                                         add_constant=False)).all())
        # no rank or constant
        with open('./fixtures/weights_scale.npy', 'rb') as f:
            weight_s = np.load(f)
        self.assertTrue((weight_s == normalize_weights(self.weights,
                                                       rank=False,
                                                       add_constant=False)).all())

        # TODO add type and size checking for inputs


class TestGetControlInputs(unittest.TestCase):
    def setUp(self):
        pass

    def test_get_control_inputs_success(self):
        # TODO defaults

        # TODO T

        # TODO B

        # TODO rho

        # TODO S

        # TODO system
        self.assertEqual(True, False)

    # TODO test that output dimensions are correct
    def test_get_control_inputs_dimensions(self):
        pass

    def test_get_control_inputs_error(self):
        # TODO system
        self.assertEqual(True, False)


class TestIntegrateU(unittest.TestCase):
    def setUp(self):
        with open('./fixtures/u.npy', 'rb') as f:
            self.u = np.load(f)

    def test_integrate_u_success(self):
        with open('./fixtures/u_int.npy', 'rb') as f:
            energy = np.load(f)
        self.assertTrue((energy == integrate_u(self.u)).all())

        # TODO different scipy version


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

        # TODO, check with different T

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
    unittest.main()
