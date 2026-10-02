"""nctpy.pipelines: the wrapper classes do what the direct calls do, in task order."""

import contextlib
import io
import unittest
from unittest import mock

import numpy as np

from nctpy.energies import get_control_inputs, integrate_u
from nctpy import pipelines
from nctpy.pipelines import ComputeControlEnergy, ComputeOptimizedControlEnergy
from nctpy.utils import matrix_normalization, normalize_state

N = 9
MODULES = np.repeat(np.arange(3), 3)


def connectome():
    rng = np.random.default_rng(0)
    A = rng.random((N, N))
    return (A + A.T) / 2


def make_tasks(**extra):
    tasks = []
    for i in range(3):
        for j in range(3):
            task = dict(
                x0=normalize_state(MODULES == i), xf=normalize_state(MODULES == j), B=np.eye(N), S=np.eye(N), rho=1
            )
            task.update(extra)
            tasks.append(task)
    return tasks


def direct_energy(A_norm, task, system, T, B=None):
    _, u, _ = get_control_inputs(
        A_norm, T, task["B"] if B is None else B, task["x0"], task["xf"], system=system, rho=task["rho"],
        S=task["S"], xr=task.get("xr", "zero"),
    )  # fmt: skip
    return np.sum(integrate_u(u))


def quietly(fn):
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        fn()
    return out.getvalue()


class TestComputeControlEnergy(unittest.TestCase):
    """Tasks sharing a system are solved together.

    A transition that completes agrees with a direct call to rounding; one that does not complete is solved again
    on its own, so it equals the direct call exactly.
    """

    def assertMatchesDirectCalls(self, pipeline, A_norm, tasks, system, T):
        self.assertEqual(pipeline.E.shape, (len(tasks),))
        for i, task in enumerate(tasks):
            _, _, err = get_control_inputs(
                A_norm, T, task["B"], task["x0"], task["xf"], system=system, rho=task["rho"], S=task["S"],
                xr=task.get("xr", "zero"),
            )  # fmt: skip
            expected = direct_energy(A_norm, task, system, T)
            if max(err) < 1e-8:
                np.testing.assert_allclose(pipeline.E[i], expected, rtol=1e-10)
            else:
                self.assertEqual(pipeline.E[i], expected)

    def test_matches_direct_calls_in_task_order(self):
        A = connectome()
        partial = np.diag((np.arange(N) % 3 > 0).astype(float))
        for system, T in (("continuous", 1), ("discrete", 3)):
            for extra in ({}, {"xr": "midpoint"}, {"rho": 1e-3, "B": partial}):  # the last does not complete
                with self.subTest(system=system, **{k: v for k, v in extra.items() if k != "B"}):
                    tasks = make_tasks(**extra)
                    pipeline = ComputeControlEnergy(A=A, control_tasks=tasks, system=system, c=1, T=T)
                    quietly(pipeline.run)
                    A_norm = matrix_normalization(A, system=system, c=1)
                    np.testing.assert_array_equal(pipeline.A_norm, A_norm)
                    self.assertMatchesDirectCalls(pipeline, A_norm, tasks, system, T)

    def test_mixed_systems_and_small_batches(self):
        # tasks alternate between two control sets; a tiny memory budget forces batches of one or two
        A = connectome()
        partial = np.diag((np.arange(N) % 3 > 0).astype(float))
        tasks = make_tasks()
        for i in range(0, len(tasks), 2):
            tasks[i]["B"] = partial
        A_norm = matrix_normalization(A, system="continuous")
        for budget in (pipelines.BATCH_BYTES, 3 * 8 * N * 1001 * 2):
            with self.subTest(budget=budget), mock.patch.object(pipelines, "BATCH_BYTES", budget):
                pipeline = ComputeControlEnergy(A=A, control_tasks=tasks, system="continuous")
                quietly(pipeline.run)
                self.assertMatchesDirectCalls(pipeline, A_norm, tasks, "continuous", 1)

    def test_unusual_tasks_are_solved_on_their_own(self):
        # S given as 'identity', or a task the batch cannot take, goes through get_control_inputs as before
        A = connectome()
        A_norm = matrix_normalization(A, system="continuous")
        tasks = make_tasks(S="identity")
        pipeline = ComputeControlEnergy(A=A, control_tasks=tasks, system="continuous")
        quietly(pipeline.run)
        expected = make_tasks()
        self.assertMatchesDirectCalls(pipeline, A_norm, expected, "continuous", 1)

        tasks[4]["xr"] = "target"  # not a reference state: fails as a direct call does
        with self.assertRaisesRegex(ValueError, "xr must be"):
            quietly(ComputeControlEnergy(A=A, control_tasks=tasks, system="continuous").run)

    def test_xr_key_is_used(self):
        A = connectome()
        energies = []
        for extra in ({}, {"xr": "midpoint"}):
            pipeline = ComputeControlEnergy(A=A, control_tasks=make_tasks(**extra), system="continuous")
            quietly(pipeline.run)
            energies.append(pipeline.E)
        self.assertFalse(np.allclose(*energies))

    def test_preset_A_norm_is_reused(self):
        A = connectome()
        pipeline = ComputeControlEnergy(A=A, control_tasks=make_tasks(), system="continuous")
        pipeline.A_norm = matrix_normalization(A, system="continuous", c=3)
        quietly(pipeline.run)
        np.testing.assert_array_equal(pipeline.A_norm, matrix_normalization(A, system="continuous", c=3))

    def test_non_square_A_raises(self):
        pipeline = ComputeControlEnergy(A=connectome()[:, :5], control_tasks=make_tasks(), system="continuous")
        with self.assertRaises(Exception) as ctx:
            pipeline.run()
        self.assertEqual(str(ctx.exception), "A matrix is not square. This routine requires A.shape[0] == A.shape[1]")


class TestComputeOptimizedControlEnergy(unittest.TestCase):
    def test_runs_without_B_and_reports_each_step(self):
        # the SI's task dict has no 'B': the class optimises B and never reads it
        task = {k: v for k, v in make_tasks()[1].items() if k != "B"}
        pipeline = ComputeOptimizedControlEnergy(A=connectome(), control_task=task, system="continuous", n_steps=3)
        printed = quietly(pipeline.run)
        self.assertEqual(printed, "".join(f"Running gradient step {i}\n" for i in range(3)))
        self.assertEqual((pipeline.E_opt.shape, pipeline.B_opt.shape), ((3,), (3, N)))

    def test_first_step(self):
        # one step from B = I: the gradient is estimated by adding 0.1 to each weight, then rescaled to ||I||.
        # Every energy uses the task's xr, if it has one (until 1.1 the class always used 'zero').
        A = connectome()
        A_norm = matrix_normalization(A, system="continuous")
        for extra in ({}, {"xr": "xf"}, {"xr": "midpoint"}):
            with self.subTest(**extra):
                task = make_tasks(**extra)[1]
                pipeline = ComputeOptimizedControlEnergy(A=A, control_task=task, system="continuous", n_steps=1)
                quietly(pipeline.run)

                E = direct_energy(A_norm, task, "continuous", 1, B=np.eye(N))
                bumps = [np.eye(N) + 0.1 * np.diag(np.arange(N) == i) for i in range(N)]
                E_d = np.array([direct_energy(A_norm, task, "continuous", 1, B=B) for B in bumps]) - E
                weights = 1 - 0.01 * E_d
                weights = weights / np.linalg.norm(weights) * np.sqrt(N)
                np.testing.assert_allclose(pipeline.B_opt[0], weights, rtol=1e-12)
                expected = direct_energy(A_norm, task, "continuous", 1, B=np.diag(weights))
                np.testing.assert_allclose(pipeline.E_opt[0], expected, rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
