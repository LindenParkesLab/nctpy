"""Control inputs, trajectories and energies for linear network dynamics.

The models are, in continuous and discrete time,

    dx/dt = A x(t) + B u(t)        and        x(t + 1) = A x(t) + B u(t),

with A normalised by :func:`nctpy.utils.matrix_normalization` for the matching time system. Equation numbers
in comments refer to Kim et al., Nat Commun 16:11639 (2025), Methods.
"""

from collections.abc import Callable
from typing import Any, TypeVar, cast

import numpy as np
import numpy.typing as npt
import scipy as sp
import scipy.integrate  # noqa: F401  (makes sp.integrate available)
import scipy.linalg as la
import scipy.sparse.linalg
from scipy import sparse

from nctpy._validation import _check_rho, _check_system
from nctpy.utils import expm

FloatArray = npt.NDArray[np.float64]
_R = TypeVar("_R")

DT = 0.001  # continuous-time integration step


class _LastCall:
    """Wrap a function so that it remembers its most recent result.

    Transitions are usually computed one after another for the same system (A_norm, B, S, rho, T), and the
    expensive parts of each depend only on the system. Wrapped, those parts are computed once for a run of
    calls on one system. Arguments are compared byte for byte, so a result is reused only for identical input.
    One result is kept per wrapped function.
    """

    def __init__(self, fn: Callable[..., Any]) -> None:
        self._fn = fn
        self._entry: tuple[tuple[Any, ...], Any] | None = None

    def __call__(self, *args: Any) -> Any:
        key = tuple((a.dtype.str, a.shape, a.tobytes()) if isinstance(a, np.ndarray) else (type(a), a) for a in args)
        entry = self._entry  # read once: another thread may replace it
        if entry is None or entry[0] != key:
            entry = (key, self._fn(*args))
            self._entry = entry
        return entry[1]


def _last_call(fn: Callable[..., _R]) -> Callable[..., _R]:
    return _LastCall(fn)


def _as_float(a: npt.ArrayLike) -> npt.NDArray[Any]:
    """`a` as an array of at least float64 precision.

    Boolean and integer arrays become float64 (True -> 1.0), lower-precision floats are promoted, and float64 or
    complex input is returned unchanged, without a copy.
    """
    a = np.asarray(a)
    return a.astype(np.result_type(a.dtype, np.float64), copy=False)


def _column(a: npt.NDArray[Any]) -> npt.NDArray[Any]:
    """A 1-D state vector as an (N, 1) column; anything else unchanged."""
    return a.reshape(-1, 1) if a.ndim == 1 else a


def sim_state_eq(
    A_norm: npt.ArrayLike, B: npt.ArrayLike, x0: npt.ArrayLike, U: npt.ArrayLike, system: str | None = None
) -> FloatArray:
    """Simulate the state trajectory driven by given inputs, with no target state and no constraints.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix.
    B : (N, N) array_like
        Control node matrix. Diagonal entries designate which nodes are control nodes and how much influence those
        nodes have on system dynamics. For example, ``B=np.eye(A_norm.shape[0])`` sets all nodes as controllers with
        equal weight (1); this is referred to as a uniform full control set.
    x0 : (N,) or (N, 1) array_like
        Initial state. Boolean states are converted to 0/1 floats.
    U : (N, T) array_like
        System inputs, one column per time point. For example, to simulate the trajectory resulting from
        stimulation, every element of U could be ``log(StimFreq)*StimAmp*StimDur``. U may also vary over time.
    system : {'continuous', 'discrete'}
        Time system that A_norm was normalised for. Continuous time is integrated with forward Euler steps of
        0.001.

    Returns
    -------
    x : (N, T) ndarray
        State trajectory: the neural activity that results from simulating the system with the above parameters.
        Note the orientation, nodes x time, which is the transpose of :func:`get_control_inputs`' output.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    """
    A_norm, B, U = _as_float(A_norm), _as_float(B), _as_float(U)
    x0 = _column(_as_float(x0))

    n_steps = np.size(U, 1)
    n_nodes = np.size(A_norm, 0)
    x = np.zeros((n_nodes, n_steps))
    xt = x0

    _check_system(system)
    if system == "continuous":
        for t in range(n_steps):
            x[:, t] = xt[:, 0]
            dxdt = A_norm @ xt + B @ np.reshape(U[:, t], (n_nodes, 1))  # state equation
            xt = dxdt * 0.001 + xt
    else:
        for t in range(n_steps):
            x[:, t] = xt[:, 0]
            xt = A_norm @ xt + B @ np.reshape(U[:, t], (n_nodes, 1))  # state equation

    return x


def get_control_inputs(
    A_norm: npt.ArrayLike,
    T: float,
    B: npt.ArrayLike,
    x0: npt.ArrayLike,
    xf: npt.ArrayLike,
    system: str | None = None,
    rho: float = 1,
    S: npt.ArrayLike | str = "identity",
    xr: npt.ArrayLike | str = "zero",
    expm_version: str = "scipy",
) -> tuple[FloatArray, FloatArray, list[np.floating[Any]]]:
    """Compute the optimal control inputs and state trajectory that drive the system from x0 to xf.

    The inputs u(t) minimise the cost ``integral of (x - xr)^T S (x - xr) + rho u^T u`` subject to the dynamics
    and to the boundary conditions x(0) = x0 and x(T) = xf.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix.
    T : float
        Time horizon: the amount of time the model runs for. Too long yields a large error; too short does not
        give enough time for control (i.e., the state transition may not complete). For discrete-time systems T
        is an integer number of steps, at least 2.
    B : (N, N) array_like
        Control node matrix. Diagonal entries designate which nodes are control nodes and how much influence those
        nodes have on system dynamics. For example, ``B=np.eye(A_norm.shape[0])`` sets all nodes as controllers with
        equal weight (1); this is referred to as a uniform full control set.
    x0 : (N,) or (N, 1) array_like
        Initial state: the initial condition of the system. Boolean states are converted to 0/1 floats.
    xf : (N,) or (N, 1) array_like
        Target state: the state the system is controlled toward and should arrive at at time T. Boolean states are
        converted to 0/1 floats.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system. Required.
    rho : float, default 1
        Mixing parameter: the weight of the control signals' cost relative to the state trajectory's, in the cost
        above. rho=1 weights them equally; a smaller rho constrains the state trajectory more, a larger one less.
        Must be > 0, and has no effect if S is all zeros (``S=np.zeros((N, N))``, minimum-energy control).
    S : (N, N) array_like or 'identity', default 'identity'
        Constraint matrix for the state trajectory. Determines which nodes in the state trajectory are
        constrained. By default, all nodes' neural activity is constrained.
    xr : (N,) or (N, 1) array_like, or {'zero', 'x0', 'xf', 'midpoint'}, default 'zero'
        Reference state. This state governs the constraints placed on the state trajectory. By default it is a
        vector of zeros, so the state trajectory is penalised for departing too far from 0 activity. 'midpoint'
        is ``x0 + (xf - x0) / 2``. This only applies if S contains nodes to constrain.
    expm_version : {'scipy', 'eig'}, default 'scipy'
        Matrix exponential used in continuous time: :func:`scipy.linalg.expm`, or :func:`nctpy.utils.expm`
        (eigendecomposition).

    Returns
    -------
    x : (n_t, N) ndarray
        State trajectory (neural activity), time x nodes. In continuous time n_t = T/0.001 + 1; in discrete time
        n_t = T + 1.
    u : (n_t, N) ndarray
        Control signals, time x nodes. In continuous time n_t = T/0.001 + 1; in discrete time n_t = T.
    err : list of two floats
        Numerical error: ``[inversion error, reconstruction error]``. The inversion error is the residual of the
        linear solve for the costate (continuous time) or for the whole trajectory (discrete time); the
        reconstruction error is how far the trajectory misses xf (continuous time) or departs from the dynamics
        (discrete time). A transition that did not complete or is ill-conditioned still returns, with large
        errors, for the caller to inspect.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options, or if T < 2 in discrete time.
    ValueError
        If rho is not positive, if xr is a string other than the four options, or if x0, xf or xr hold more
        than one state.

    Notes
    -----
    Everything that depends only on the system (A_norm, T, B, S, rho), not on the states, is computed once and
    reused while consecutive calls share that system, e.g. in a loop over transitions. Results are identical either
    way: reuse happens only when those inputs are identical byte for byte.
    """
    A_norm, B = _as_float(A_norm), _as_float(B)
    n_nodes = A_norm.shape[0]

    X0 = _column(_as_float(x0))
    XF = _column(_as_float(xf))
    XR = _reference(xr, X0, XF, n_nodes)
    if isinstance(S, str):
        if S == "identity":
            S = np.eye(n_nodes)
    else:
        S = _as_float(S)

    _check_system(system)
    _check_rho(rho)
    for name, state in (("x0", X0), ("xf", XF), ("xr", XR)):
        if not isinstance(state, str) and state.ndim == 2 and state.shape[1] != 1:
            raise ValueError(f"{name} must be a single state, of shape (N,) or (N, 1); got shape {state.shape}")

    x, u, err = _control_inputs(A_norm, T, B, X0, XF, XR, S, system, rho, expm_version)
    return x[:, :, 0], u[:, :, 0], err[0]


def _reference(xr: npt.ArrayLike | str, x0: FloatArray, xf: FloatArray, n_nodes: int) -> Any:
    """The reference state as a column (or columns, for columns of states). An unknown string is returned as is."""
    if not isinstance(xr, str):
        return _column(_as_float(xr))
    if xr == "x0":
        return x0
    if xr == "xf":
        return xf
    if xr == "zero":
        return np.zeros((n_nodes, x0.shape[1]))
    if xr == "midpoint":
        return x0 + ((xf - x0) * 0.5)
    return xr


def _check_reference(XR: Any) -> None:
    """Raise if the reference state is a string that _reference did not recognise."""
    if isinstance(XR, str):
        raise ValueError(f"xr must be 'zero', 'x0', 'xf', 'midpoint' or a state of shape (N,) or (N, 1); got {XR!r}")


def _control_inputs(
    A_norm: FloatArray,
    T: float,
    B: FloatArray,
    X0: FloatArray,
    XF: FloatArray,
    XR: Any,
    S: Any,
    system: str | None,
    rho: float,
    expm_version: str,
) -> tuple[FloatArray, FloatArray, list[list[np.floating[Any]]]]:
    """Solve k transitions of one system at once: column j of X0, XF and XR is transition j.

    Returns x and u as (n_t, N, k) arrays, time x nodes x transitions, and err as k pairs
    [inversion error, reconstruction error]. With k = 1 this is exactly get_control_inputs.
    """
    _check_system(system)
    _check_rho(rho)
    n_nodes, k = X0.shape
    if system == "continuous":
        M, E, E_dt = _continuous_system(A_norm, T, B, S, rho, expm_version)  # Eq. 6, e^{MT}, e^{M DT}

        # Eq. 8: [x(t); p(t)] = e^{Mt} [x0; p0] + (e^{Mt} - I) ref_offset
        _check_reference(XR)
        ref_input = np.concatenate((np.zeros((n_nodes, k)), 2 * np.dot(S, XR)), axis=0)
        ref_offset = np.linalg.solve(M, ref_input)

        # Eq. 9: the top block row of E = e^{MT} maps [x0; p0] to x(T). Solve it for the initial costate p0.
        r = np.arange(n_nodes)
        E11 = E[r, :][:, r]
        E12 = E[r, :][:, r + n_nodes]
        b1 = np.dot(np.concatenate((E11 - np.eye(n_nodes), E12), axis=1), ref_offset)
        p0_rhs = XF - np.dot(E11, X0) - b1  # E12 p0 = xf - E11 x0 - b1
        P0 = np.linalg.solve(E12, p0_rhs)

        # Integrate the state-costate system exactly over steps of DT, with E_dt = e^{M DT}
        n_steps = int(np.round(T / DT))
        z = np.zeros((2 * n_nodes, n_steps + 1, k))  # [x; p] x time x transition, as before batching for k = 1
        z[:, 0, :] = np.concatenate((X0, P0), axis=0)
        offset_dt = np.dot((E_dt - np.eye(2 * n_nodes)), ref_offset)
        # a single transition steps a (strided) vector, not a one-column matrix: BLAS rounds the two differently
        steps, offset = (z[:, :, 0], offset_dt[:, 0]) if k == 1 else (z, offset_dt)
        for i in range(1, n_steps + 1):
            steps[:, i] = np.dot(E_dt, steps[:, i - 1]) + offset

        # Extract state and input from the joint state-costate trajectory, as time x nodes x transitions
        x = z[:n_nodes].transpose(1, 0, 2)
        u = np.dot(-B.T, z[n_nodes:].reshape(n_nodes, -1)) / (2 * rho)
        u = u.reshape(n_nodes, n_steps + 1, k).transpose(1, 0, 2)

        # Collect error
        costate_residual = np.dot(E12, P0) - p0_rhs
        err = [
            [np.linalg.norm(costate_residual[:, [j]]), np.linalg.norm(x[-1, :, j].reshape(-1, 1) - XF[:, [j]])]
            for j in range(k)
        ]
        return x, u, err

    if T <= 1:
        raise Exception("Discrete time systems must have T >= 2")
    T = cast(int, T)
    M_sparse, M_lu = _discrete_system(A_norm, T, B, S, rho)  # the system matrix and its LU

    # Right-hand side: A x0 in the first state row, -xf in the last, -state_cost xr in every costate row
    boundary = np.concatenate(
        (np.dot(A_norm, X0), np.zeros((n_nodes * (T - 2), k)), -XF, np.zeros((n_nodes * (T - 1), k))), axis=0
    )
    _check_reference(XR)
    reference = np.concatenate((np.zeros((n_nodes * T, k)), np.tile(2 * np.dot(S, XR), (T - 1, 1))), axis=0)
    b = boundary - reference

    # Solve the simultaneous state and costate equations, one transition at a time with the shared LU
    xs, us, err = [], [], []
    for j in range(k):
        v = M_lu.solve(b[:, j])
        V = v.reshape((n_nodes, int(len(v) / n_nodes)), order="F")
        x = np.concatenate((X0[:, [j]], V[:, : T - 1], XF[:, [j]]), axis=1)
        u = np.dot(-B.T, V[:, T - 1 :]) / (2 * rho)
        residual = np.dot(M_sparse, sparse.csc_matrix(np.expand_dims(v, axis=1))) - sparse.csc_matrix(b[:, [j]])
        err_traj = np.linalg.norm(x[:, 1:] - (np.dot(A_norm, x[:, 0:-1]) + np.dot(B, u)))
        xs.append(x.T)
        us.append(u.T)
        err.append([np.linalg.norm(residual.todense()), err_traj])
    return np.stack(xs, axis=2), np.stack(us, axis=2), err


def integrate_u(u: npt.ArrayLike) -> FloatArray:
    """Integrate squared control inputs over time, with Simpson's rule, to give the energy at each node.

    If the control set (B) is the identity this gives energies nearly identical to a Riemann sum. When control
    sets are sparse the inputs can be very curved, and Simpson's rule is more accurate.

    Samples are taken to be one unit apart: the result is not scaled by the time step.

    Parameters
    ----------
    u : (n_t, N) array_like
        Control signals input to the system, time x nodes, as returned by :func:`get_control_inputs`.

    Returns
    -------
    energy : (N,) ndarray
        Energy input into each node.
    """
    u = np.asarray(u)
    return sp.integrate.simpson(u.T**2)


def gramian(A_norm: npt.ArrayLike, T: float, system: str | None = None) -> FloatArray | float:
    """Compute the controllability Gramian of (A_norm, I).

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix.
    T : float
        Time horizon. ``np.inf`` gives the infinite-horizon Gramian. For discrete-time systems a finite T is an
        integer number of steps.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system. Required.

    Returns
    -------
    Wc : (N, N) ndarray, or float
        The Gramian. With finite T, in continuous time it is the integral of ``e^{At} e^{A^T t}`` over [0, T]
        by Simpson's rule over steps of 0.001; in discrete time it is ``I + sum over k = 1..T of A^k (A^k)^T``.
        With T = np.inf it solves the Lyapunov equation if the system is stable; if it is not, it prints a
        message and returns ``np.nan``.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    """
    _check_system(system)
    A_norm = _as_float(A_norm)
    n_nodes = A_norm.shape[0]
    eye = np.eye(n_nodes)

    # Only a stable system has an infinite-horizon Gramian: it solves a Lyapunov equation
    if T == np.inf:
        eigvals = np.linalg.eig(A_norm)[0]
        stable = np.max(np.real(eigvals)) < 0 if system == "continuous" else np.max(np.abs(eigvals)) < 1
        if not stable:
            print("cannot compute infinite-time Gramian for an unstable system!")
            return np.nan
        if system == "continuous":
            return la.solve_continuous_lyapunov(A_norm, -eye)
        return la.solve_discrete_lyapunov(A_norm, eye)

    if system == "continuous":
        step = 0.001
        t = np.arange(0, (T + step / 2), step)
        e_step = sp.linalg.expm(A_norm * step)
        e_At = eye  # e^{A t}, accumulated one step at a time
        integrand = np.zeros((n_nodes, n_nodes, len(t)))
        integrand[:, :, 0] = eye
        for i in range(1, len(t)):
            e_At = e_At @ e_step
            integrand[:, :, i] = e_At @ e_At.T
        return sp.integrate.simpson(integrand, x=t, axis=2)

    A_k = eye  # A^k
    Wc = eye
    for _ in range(cast(int, T)):
        A_k = A_k @ A_norm
        Wc = Wc + A_k @ A_k.T
    return Wc


def minimum_energy_fast(
    A_norm: npt.ArrayLike, T: float, B: npt.ArrayLike, x0: npt.ArrayLike, xf: npt.ArrayLike
) -> FloatArray:
    """Compute the minimum control energy from x0 to xf in continuous time, from the controllability Gramian.

    The Gramian of (A_norm, B) over [0, T] is integrated with Simpson's rule over 1000 steps, and the energy is
    split across nodes as ``(pinv(Wc) d) * d``, where ``d = xf - e^{A T} x0``. It equals minimum-energy control
    (``get_control_inputs`` with S all zeros) without simulating the trajectory.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix, for a continuous-time system.
    T : float
        Time horizon.
    B : (N, N) array_like
        Control node matrix.
    x0 : (N,) or (N, k) array_like
        Initial state, or k initial states as columns. Boolean states are converted to 0/1 floats.
    xf : (N,) or (N, k) array_like
        Target state, or k target states as columns, paired with those of x0.

    Returns
    -------
    energy : (N, 1) or (N, k) ndarray
        Energy at each node, one column per transition; sum a column for the total.

    Notes
    -----
    The Gramian and e^{AT} depend only on (A_norm, T, B). They are computed once and reused while consecutive calls
    share them, so looping over transitions is fast; passing all transitions at once as columns is equivalent.
    """
    A_norm, B = _as_float(A_norm), _as_float(B)
    x0 = _column(_as_float(x0))
    xf = _column(_as_float(xf))

    G_pinv, eAT = _minimum_energy_system(A_norm, T, B)

    delx = xf - eAT @ x0
    return np.multiply(G_pinv @ delx, delx)


def minimum_energy_infinite(
    A_norm: npt.ArrayLike, B: npt.ArrayLike, xf: npt.ArrayLike, system: str | None = None
) -> tuple[FloatArray, FloatArray]:
    """Compute the minimum control energy to reach xf over an infinite time horizon.

    Uses the infinite-horizon controllability Gramian Wc of (A_norm, B), the solution of a Lyapunov equation, so it
    needs no numerical integration or matrix exponential (Kim et al., Nat Commun 2025). The minimum energy is
    ``xf^T Wc^-1 xf``. Over an infinite horizon the initial state's contribution decays away, so the energy depends
    only on the target.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix. It must be stable: in continuous time every eigenvalue has a
        negative real part, in discrete time every eigenvalue lies inside the unit circle.
    B : (N, N) array_like
        Control node matrix.
    xf : (N,) or (N, k) array_like
        Target state, or k target states as columns. Boolean states are converted to 0/1 floats.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system. Required.

    Returns
    -------
    energy : (N, 1) or (N, k) ndarray
        Energy at each node, ``(Wc^-1 xf) * xf``, one column per target; sum a column for the total. All NaN if
        A_norm is not stable, since the infinite-horizon Gramian does not exist.
    xf_reached : (N, 1) or (N, k) ndarray
        ``Wc Wc^-1 xf``: the target as reconstructed through the inverse Gramian. Its difference from xf measures
        how accurately Wc could be inverted, which degrades as control sets get sparser (Kim et al., Fig. 3A-B).

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.

    Notes
    -----
    With few control nodes Wc is close to singular, and ``xf^T Wc^-1 xf`` can lose all accuracy (even its sign);
    check ``xf_reached`` against xf. See also :func:`average_energy_infinite`, which does not depend on a target.
    """
    _check_system(system)
    A_norm, B = _as_float(A_norm), _as_float(B)
    xf = _column(_as_float(xf))
    Wc = _infinite_gramian(A_norm, B, system)
    if Wc is None:
        nan = np.full(xf.shape, np.nan)
        return nan, nan.copy()
    Wc_inv_xf = np.linalg.solve(Wc, xf)
    return np.multiply(Wc_inv_xf, xf), Wc @ Wc_inv_xf


def average_energy_infinite(A_norm: npt.ArrayLike, B: npt.ArrayLike, system: str | None = None) -> float:
    """Compute the average minimum control energy over an infinite time horizon: the trace of the inverse Gramian.

    ``trace(Wc^-1)``, where Wc is the infinite-horizon controllability Gramian of (A_norm, B). It equals the sum of the
    minimum energies to reach each unit vector, and is proportional to the average minimum energy over all
    unit-norm target states, so it does not depend on any particular transition. This is the control energy shown in
    Kim et al. (Nat Commun 2025), Fig. 3C, there multiplied by dt = 0.001.

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix. It must be stable (see :func:`minimum_energy_infinite`).
    B : (N, N) array_like
        Control node matrix.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system. Required.

    Returns
    -------
    energy : float
        ``trace(Wc^-1)``, computed as the sum of the reciprocal singular values of Wc. NaN if A_norm is not stable.

    Raises
    ------
    Exception
        If `system` is missing or not one of the two options.
    """
    _check_system(system)
    Wc = _infinite_gramian(_as_float(A_norm), _as_float(B), system)
    if Wc is None:
        return np.nan
    return np.sum(1 / la.svdvals(Wc))


# Everything below depends only on the system, not on the states, and is reused across transitions (_last_call).


@_last_call
def _continuous_system(
    A_norm: FloatArray, T: float, B: FloatArray, S: FloatArray, rho: float, expm_version: str
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """M (Eq. 6), e^{MT} and e^{M DT}.

    The state x and the costate p evolve jointly, d/dt [x; p] = M [x; p] + [0; 2 S xr], with
    M = [[A, -B B^T / (2 rho)], [-2 S, -A^T]].
    """
    costate_to_state = np.dot(-B, B.T) / (2 * rho)  # B u(t) = costate_to_state p(t)
    M = np.concatenate(
        (np.concatenate((A_norm, costate_to_state), axis=1), np.concatenate((-2 * S, -A_norm.T), axis=1)),
        axis=0,
    )
    if expm_version == "scipy":
        E, E_dt = sp.linalg.expm(M * T), sp.linalg.expm(M * DT)
    elif expm_version == "eig":
        E, E_dt = expm(M * T), expm(M * DT)
    return M, E, E_dt


@_last_call
def _discrete_system(A_norm: FloatArray, T: int, B: FloatArray, S: FloatArray, rho: float) -> tuple[Any, Any]:
    """The discrete-time system matrix (sparse) and its LU factorisation."""
    n_nodes = A_norm.shape[0]

    # Solve for every unknown at once. The unknowns are the free states x(1)..x(T-1) and the costates
    # p(0)..p(T-1), stacked in blocks of n_nodes, and the inputs are u(t) = -B^T p(t) / (2 rho). The rows are
    #   state equations,   t = 0..T-1:  x(t+1) = A x(t) + costate_to_state p(t)
    #   costate equations, t = 1..T-1:  p(t-1) = A^T p(t) + state_cost (x(t) - xr)
    # with the known x(0) = x0 and x(T) = xf moved to the right-hand side.
    costate_to_state = np.dot(-B, B.T) / (2 * rho)  # B u(t) = costate_to_state p(t)
    state_cost = 2 * S
    eye = np.eye(n_nodes)
    block = np.arange(n_nodes)

    def x_col(t: int) -> npt.NDArray[np.int_]:
        return block + (t - 1) * n_nodes  # x(1) is the first block

    def p_col(t: int) -> npt.NDArray[np.int_]:
        return block + (T - 1 + t) * n_nodes

    M = np.zeros(((2 * T - 1) * n_nodes, (2 * T - 1) * n_nodes))
    for t in range(T):
        row = block + t * n_nodes
        M[np.ix_(row, p_col(t))] = -costate_to_state
        if t + 1 < T:
            M[np.ix_(row, x_col(t + 1))] = eye
        if t > 0:
            M[np.ix_(row, x_col(t))] = -A_norm
    for t in range(1, T):
        row = block + (T - 1 + t) * n_nodes
        M[np.ix_(row, x_col(t))] = -state_cost
        M[np.ix_(row, p_col(t - 1))] = eye
        M[np.ix_(row, p_col(t))] = -A_norm.T

    M_sparse = sparse.csc_matrix(M)
    return M_sparse, sparse.linalg.splu(M_sparse)


@_last_call
def _minimum_energy_system(A_norm: FloatArray, T: float, B: FloatArray) -> tuple[FloatArray, FloatArray]:
    """The pseudo-inverse of the controllability Gramian of (A_norm, B) over [0, T], and e^{AT}."""
    n_nodes = A_norm.shape[0]

    # Number of integration steps
    nt = 1000
    dt = T / nt

    # Numerical integration with Simpson's 1/3 rule
    dE = sp.linalg.expm(A_norm * dt)  # integration step
    dEA = np.eye(n_nodes)  # accumulates expm(A * dt)
    G = np.zeros((n_nodes, n_nodes))  # Gramian

    for _ in range(1, nt // 2):
        # Add odd terms
        dEA = dEA @ dE
        p1 = dEA @ B
        # Add even terms
        dEA = dEA @ dE
        p2 = dEA @ B
        G = G + 4 * (p1 @ p1.T) + 2 * (p2 @ p2.T)

    # Add final odd term
    dEA = dEA @ dE
    p1 = dEA @ B
    G = G + 4 * (p1 @ p1.T)

    # Add the end points and scale by the step
    eAT = sp.linalg.expm(A_norm * T)
    G = (G + B @ B.T + (eAT @ B) @ (eAT @ B).T) * dt / 3
    return np.linalg.pinv(G), eAT


@_last_call
def _infinite_gramian(A_norm: FloatArray, B: FloatArray, system: str) -> FloatArray | None:
    """The infinite-horizon controllability Gramian of (A_norm, B); None if A_norm is not stable."""
    eigvals = np.linalg.eigvals(A_norm)
    if system == "continuous":
        if np.max(eigvals.real) >= 0:
            return None
        return la.solve_continuous_lyapunov(A_norm, -(B @ B.T))
    if np.max(np.abs(eigvals)) >= 1:
        return None
    return la.solve_discrete_lyapunov(A_norm, B @ B.T)
