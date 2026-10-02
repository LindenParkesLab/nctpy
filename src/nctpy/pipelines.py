"""Wrapper classes that normalise a connectome and compute control energy for one or many control tasks.

Both take the raw adjacency matrix A, normalise it with :func:`nctpy.utils.matrix_normalization` when ``run()``
is called, and compute energy as the sum over nodes of :func:`nctpy.energies.integrate_u` applied to the inputs
from :func:`nctpy.energies.get_control_inputs`.
"""

from typing import Any, cast

import numpy as np
import numpy.typing as npt
import scipy.linalg as la
from tqdm import tqdm

from nctpy.energies import DT, _as_float, _column, _control_inputs, _reference, get_control_inputs, integrate_u
from nctpy.utils import matrix_normalization

COMPLETE = 1e-8  # a transition completed if both its error terms are below this (the protocol paper's threshold)
BATCH_BYTES = 256e6  # memory budget for one batch of trajectories in ComputeControlEnergy


class _ControlEnergyPipeline:
    """The normalisation and energy computation the pipeline classes share."""

    A: Any
    system: str | None
    c: float
    T: float
    n_nodes: int
    A_norm: npt.NDArray[np.float64]

    def _check_inputs(self) -> None:
        if self.A.shape[0] != self.A.shape[1]:
            raise Exception("A matrix is not square. This routine requires A.shape[0] == A.shape[1]")
        self.n_nodes = self.A.shape[0]
        # normalise once; an A_norm already set on the instance (e.g. by an earlier run()) is reused
        if not hasattr(self, "A_norm"):
            self.A_norm = matrix_normalization(self.A, system=self.system, c=self.c)

    def _energy(self, task: dict[str, Any], B: npt.ArrayLike, xr: npt.ArrayLike | str = "zero") -> np.float64:
        _, u, _ = get_control_inputs(
            A_norm=self.A_norm,
            T=self.T,
            B=B,
            x0=task["x0"],
            xf=task["xf"],
            system=self.system,
            rho=task["rho"],
            S=task["S"],
            xr=xr,
        )
        return np.sum(integrate_u(u))


class ComputeControlEnergy(_ControlEnergyPipeline):
    """Compute control energy for a set of independent control tasks.

    See /scripts/control_energy_wrapper.ipynb for an example.

    Parameters
    ----------
    A : (N, N) ndarray
        Adjacency matrix representing a structural connectome (not normalised).
    control_tasks : list of dict
        One dict per control task, with the keys ``'x0'``, ``'xf'``, ``'B'``, ``'S'`` and ``'rho'``, and
        optionally ``'xr'`` (default ``'zero'``); see :func:`nctpy.energies.get_control_inputs`. For example::

            control_tasks = []  # list of tasks

            control_task = dict()  # initialize dict
            control_task['x0'] = x0  # store initial state
            control_task['xf'] = xf  # store target state
            control_task['B'] = B  # store control nodes
            control_task['S'] = S  # store state trajectory constraints
            control_task['rho'] = rho  # store rho
            control_tasks.append(control_task)

    system : {'continuous', 'discrete'}
        Time system to normalise A for. Required.
    c : float, default 1
        Normalisation constant.
    T : float, default 1
        Time horizon.

    Attributes
    ----------
    E : (n_tasks,) ndarray
        Set by ``run()``: the total control energy of each task, in task order.
    A_norm : (N, N) ndarray
        Set by ``run()``: the normalised adjacency matrix.
    """

    def __init__(
        self, A: Any, control_tasks: list[dict[str, Any]], system: str | None = None, c: float = 1, T: float = 1
    ) -> None:
        self.A = A
        self.control_tasks = control_tasks
        self.system = system
        self.c = c
        self.T = T

    def run(self) -> None:
        """Compute the energy of every task and store it in ``self.E``.

        Consecutive tasks that share a system (B, S and rho) are solved together. A transition that does not
        complete (an error term of 1e-8 or more) is solved again on its own, so its energy is exactly what
        :func:`nctpy.energies.get_control_inputs` gives; completed transitions agree with it to rounding.
        """
        self._check_inputs()
        tasks = [self._prepare(task) for task in self.control_tasks]
        n_t = int(np.round(self.T / DT)) + 1 if self.system == "continuous" else int(self.T) + 1
        max_batch = max(1, int(BATCH_BYTES // (3 * 8 * self.n_nodes * n_t)))

        E = np.empty(len(tasks))
        with tqdm(total=len(tasks)) as progress:
            start = 0
            while start < len(tasks):
                stop = start + 1
                while stop < len(tasks) and stop - start < max_batch and _same_system(tasks[start], tasks[stop]):
                    stop += 1
                E[start:stop] = self._batch_energy(start, stop, tasks)
                progress.update(stop - start)
                start = stop
        self.E = E

    def _prepare(self, task: dict[str, Any]) -> tuple[Any, ...] | None:
        """A task's system key, B, S and states as columns; None if it should be solved on its own."""
        B, S = task["B"], task["S"]  # read in the order a single solve reads them, so errors are the same
        x0, xf = _column(_as_float(task["x0"])), _column(_as_float(task["xf"]))
        xr = _reference(task.get("xr", "zero"), x0, xf, self.n_nodes)
        if isinstance(S, str) and S == "identity":
            S = np.eye(self.n_nodes)
        if isinstance(S, str) or any(isinstance(a, str) or a.shape != (self.n_nodes, 1) for a in (x0, xf, xr)):
            return None
        B, S = _as_float(B), _as_float(S)
        key = (B.dtype.str, B.shape, B.tobytes(), S.dtype.str, S.shape, S.tobytes(), type(task["rho"]), task["rho"])
        return key, B, S, x0, xf, xr

    def _batch_energy(self, start: int, stop: int, tasks: list[tuple[Any, ...] | None]) -> list[np.float64]:
        def alone(i: int) -> np.float64:
            task = self.control_tasks[i]
            return self._energy(task, B=task["B"], xr=task.get("xr", "zero"))

        if stop - start == 1:
            return [alone(start)]
        _, B, S, *_ = cast(tuple[Any, ...], tasks[start])
        X0, XF, XR = (
            np.concatenate([cast(tuple[Any, ...], tasks[i])[c] for i in range(start, stop)], axis=1) for c in (3, 4, 5)
        )
        _, u, err = _control_inputs(
            self.A_norm, self.T, B, X0, XF, XR, S, self.system, self.control_tasks[start]["rho"], "scipy"
        )
        return [
            np.sum(integrate_u(u[:, :, j])) if all(e < COMPLETE for e in err[j]) else alone(start + j)
            for j in range(stop - start)
        ]


def _same_system(a: tuple[Any, ...] | None, b: tuple[Any, ...] | None) -> bool:
    return a is not None and b is not None and a[0] == b[0]


class ComputeOptimizedControlEnergy(_ControlEnergyPipeline):
    """Compute `optimized` control energy for a single control task, by gradient descent on the control weights.

    Optimisation starts from a uniform full control set (B = I). At each step, the energy's sensitivity to each
    node's weight is estimated by adding 0.1 to that weight; the weights take a step of size `lr` down that
    gradient and are rescaled to the Frobenius norm of the identity. See
    /scripts/path_a_control_energy_binary.ipynb for an example.

    Parameters
    ----------
    A : (N, N) ndarray
        Adjacency matrix representing a structural connectome (not normalised).
    control_task : dict
        Control task with the keys ``'x0'``, ``'xf'``, ``'S'`` and ``'rho'``, and optionally ``'xr'`` (default
        ``'zero'``); see :func:`nctpy.energies.get_control_inputs`. A ``'B'`` key is not needed and is ignored, since
        the control weights are what is optimised. For example::

            control_task = dict()  # initialize dict
            control_task['x0'] = x0  # store initial state
            control_task['xf'] = xf  # store target state
            control_task['S'] = S  # store state trajectory constraints
            control_task['rho'] = rho  # store rho

    system : {'continuous', 'discrete'}
        Time system to normalise A for. Required.
    c : float, default 1
        Normalisation constant.
    T : float, default 1
        Time horizon.
    n_steps : int, default 2
        Number of gradient steps.
    lr : float, default 0.01
        Learning rate.

    Attributes
    ----------
    E_opt : (n_steps,) ndarray
        Set by ``run()``: the energy with the optimised weights after each step.
    B_opt : (n_steps, N) ndarray
        Set by ``run()``: the optimised control weights (the diagonal of B) after each step.
    A_norm : (N, N) ndarray
        Set by ``run()``: the normalised adjacency matrix.
    """

    def __init__(
        self,
        A: Any,
        control_task: dict[str, Any],
        system: str | None = None,
        c: float = 1,
        T: float = 1,
        n_steps: int = 2,
        lr: float = 0.01,
    ) -> None:
        self.A = A
        self.control_task = control_task
        self.system = system
        self.c = c
        self.T = T
        self.n_steps = n_steps
        self.lr = lr

    def _get_energy(self, B: npt.ArrayLike) -> np.float64:
        return self._energy(self.control_task, B=B, xr=self.control_task.get("xr", "zero"))

    def _get_energy_perturbed(self, B: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Energy with 0.1 added to each node's control weight in turn."""
        E_p = np.zeros(self.n_nodes)
        for i in tqdm(range(self.n_nodes)):
            B_p = B.copy()
            B_p[i, i] += 0.1
            E_p[i] = self._get_energy(B=B_p)
        return E_p

    def run(self) -> None:
        """Run the gradient descent and store the results in ``self.E_opt`` and ``self.B_opt``."""
        self._check_inputs()

        E_opt = np.zeros(self.n_steps)
        B_opt = np.zeros((self.n_steps, self.n_nodes))
        B_I = np.eye(self.n_nodes)

        for i in range(self.n_steps):
            print("Running gradient step {0}".format(i))

            B = B_I if i == 0 else np.diag(B_opt[i - 1, :])  # start from identity, then the previous step
            E = self._get_energy(B=B)
            E_d = self._get_energy_perturbed(B=B) - E  # energy delta per node
            B_o = np.diag(B.diagonal() - (E_d * self.lr))  # step down gradient
            B_o = B_o / la.norm(B_o) * la.norm(B_I)  # normalize

            E_opt[i] = self._get_energy(B=B_o)
            B_opt[i, :] = B_o.diagonal()

        self.E_opt = E_opt
        self.B_opt = B_opt
