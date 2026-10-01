"""nctpy.metrics: input handling (Roadmap 2.3a).

The formulas themselves are checked against theory in test_properties.py and against stored values in
test_regression.py.
"""

import unittest

import numpy as np

from nctpy.metrics import ave_control, modal_control
from nctpy.utils import matrix_normalization

N = 10


def connectome():
    rng = np.random.default_rng(0)
    A = rng.random((N, N))
    A = (A + A.T) / 2
    np.fill_diagonal(A, 0)
    return A


class TestInputs(unittest.TestCase):
    """Array inputs are computed in float64 whatever their precision."""

    def test_lower_precision_input_computed_in_float64(self):
        for system in ("continuous", "discrete"):
            A32 = matrix_normalization(connectome(), system=system).astype(np.float32)
            metrics = [("ave_control", lambda A, s=system: ave_control(A, system=s))]
            if system == "discrete":
                metrics.append(("modal_control", modal_control))
            for name, fn in metrics:
                with self.subTest(system=system, metric=name):
                    got = fn(A32)
                    self.assertEqual(got.dtype, np.float64)
                    np.testing.assert_array_equal(got, fn(A32.astype(np.float64)))

    def test_integer_input(self):
        A_int = np.round(10 * matrix_normalization(connectome(), system="discrete")).astype(int)
        np.testing.assert_array_equal(
            ave_control(A_int, system="discrete"), ave_control(A_int.astype(float), "discrete")
        )
        np.testing.assert_array_equal(modal_control(A_int), modal_control(A_int.astype(float)))

    def test_one_value_per_node(self):
        for system in ("continuous", "discrete"):
            with self.subTest(system=system):
                self.assertEqual(ave_control(matrix_normalization(connectome(), system=system), system).shape, (N,))
        self.assertEqual(modal_control(matrix_normalization(connectome(), system="discrete")).shape, (N,))


if __name__ == "__main__":
    unittest.main()
