"""Control inputs, trajectories and energies for linear network dynamics.

The models are, in continuous and discrete time,

    dx/dt = A x(t) + B u(t)        and        x(t + 1) = A x(t) + B u(t),

with A normalised by :func:`nctpy.utils.matrix_normalization` for the matching time system. Equation numbers
in comments refer to Kim et al., Nat Commun 16:11639 (2025), Methods.
"""

from typing import Any, cast

import numpy as np
import numpy.typing as npt
import scipy as sp
import scipy.integrate  # noqa: F401  (makes sp.integrate available)
import scipy.linalg as la
from scipy import sparse

from nctpy._validation import _check_rho, _check_system
from nctpy.utils import expm

FloatArray = npt.NDArray[np.float64]


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
    elif system == "discrete":
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
        Mixing parameter. Determines the extent to which the state trajectory is constrained alongside the
        control signals. rho=1 equals maximum constraint. Must be > 0, and has no effect if S is all zeros
        (``S=np.zeros((N, N))``, minimum-energy control).
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
        If rho is not positive.
    """
    A_norm, B = _as_float(A_norm), _as_float(B)
    n_nodes = A_norm.shape[0]

    x0 = _column(_as_float(x0))
    xf = _column(_as_float(xf))

    if isinstance(xr, str):
        if xr == "x0":
            xr = x0
        elif xr == "xf":
            xr = xf
        elif xr == "zero":
            xr = np.zeros((n_nodes, 1))
        elif xr == "midpoint":
            xr = x0 + ((xf - x0) * 0.5)
    else:
        xr = _column(_as_float(xr))

    if isinstance(S, str):
        if S == "identity":
            S = np.eye(n_nodes)
    else:
        S = _as_float(S)

    _check_system(system)
    _check_rho(rho)
    if system == "continuous":
        dt = 0.001
        xr, S = cast(FloatArray, xr), cast(FloatArray, S)

        # Eq. 6: the state x and the costate p evolve jointly, d/dt [x; p] = M [x; p] + ref_input
        costate_to_state = np.dot(-B, B.T) / (2 * rho)  # B u(t) = costate_to_state p(t)
        M = np.concatenate(
            (np.concatenate((A_norm, costate_to_state), axis=1), np.concatenate((-2 * S, -A_norm.T), axis=1)),
            axis=0,
        )
        ref_input = np.concatenate((np.zeros((n_nodes, 1)), 2 * np.dot(S, xr)), axis=0)

        # Eq. 8: [x(t); p(t)] = e^{Mt} [x0; p0] + (e^{Mt} - I) ref_offset
        ref_offset = np.linalg.solve(M, ref_input)

        # Eq. 9: the top block row of e^{MT} maps [x0; p0] to x(T). Solve it for the initial costate p0.
        if expm_version == "scipy":
            E = sp.linalg.expm(M * T)
        elif expm_version == "eig":
            E = expm(M * T)
        r = np.arange(n_nodes)
        E11 = E[r, :][:, r]
        E12 = E[r, :][:, r + n_nodes]
        b1 = np.dot(np.concatenate((E11 - np.eye(n_nodes), E12), axis=1), ref_offset)
        p0_rhs = xf - np.dot(E11, x0) - b1  # E12 p0 = xf - E11 x0 - b1
        p0 = np.linalg.solve(E12, p0_rhs)

        # Integrate the state-costate system exactly over steps of dt
        n_steps = int(np.round(T / dt))
        z = np.zeros((2 * n_nodes, n_steps + 1))
        z[:, 0] = np.concatenate((x0, p0), axis=0).flatten()
        if expm_version == "scipy":
            E_dt = sp.linalg.expm(M * dt)
        elif expm_version == "eig":
            E_dt = expm(M * dt)
        offset_dt = np.dot((E_dt - np.eye(2 * n_nodes)), ref_offset).flatten()
        for i in np.arange(1, n_steps + 1):
            z[:, i] = np.dot(E_dt, z[:, i - 1]) + offset_dt

        # Extract state and input from the joint state-costate trajectory
        x = z[r, :]
        u = np.dot(-B.T, z[r + n_nodes, :]) / (2 * rho)

        # Collect error
        err_costate = np.linalg.norm(np.dot(E12, p0) - p0_rhs)
        err_xf = np.linalg.norm(x[:, -1].reshape(-1, 1) - xf)
        err = [err_costate, err_xf]

        return x.T, u.T, err
    elif system == "discrete":
        if T <= 1:
            raise Exception("Discrete time systems must have T >= 2")
        T = cast(int, T)
        xr, S = cast(FloatArray, xr), cast(FloatArray, S)

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

        # Right-hand side: A x0 in the first state row, -xf in the last, -state_cost xr in every costate row
        boundary = np.concatenate(
            (np.dot(A_norm, x0), np.zeros((n_nodes * (T - 2), 1)), -xf, np.zeros((n_nodes * (T - 1), 1))), axis=0
        )
        reference = np.concatenate((np.zeros((n_nodes * T, 1)), np.tile(2 * np.dot(S, xr), (T - 1, 1))), axis=0)
        b = boundary - reference

        # Solve simultaneous state and costate equations
        M_sparse = sparse.csc_matrix(M)
        b_sparse = sparse.csc_matrix(b)
        v = sparse.linalg.spsolve(M_sparse, b_sparse)
        V = v.reshape((n_nodes, int(len(v) / n_nodes)), order="F")
        x = np.concatenate((x0, V[:, : T - 1], xf), axis=1)
        u = np.dot(-B.T, V[:, T - 1 :]) / (2 * rho)

        # Collect error
        residual = np.dot(M_sparse, sparse.csc_matrix(np.expand_dims(v, axis=1))) - b_sparse
        err_system = np.linalg.norm(residual.todense())
        err_traj = np.linalg.norm(x[:, 1:] - (np.dot(A_norm, x[:, 0:-1]) + np.dot(B, u)))
        err = [err_system, err_traj]

        return x.T, u.T, err
    raise AssertionError("unreachable: _check_system accepts only 'continuous' and 'discrete'")


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


def gramian(A_norm: npt.ArrayLike, T: float, system: str | None = None) -> FloatArray | float | None:
    """Compute the controllability Gramian of (A_norm, I).

    Parameters
    ----------
    A_norm : (N, N) array_like
        Normalised structural connectivity matrix.
    T : float
        Time horizon. ``np.inf`` gives the infinite-horizon Gramian. For discrete-time systems a finite T is an
        integer number of steps.
    system : {'continuous', 'discrete'}
        Whether A_norm was normalised for a continuous-time or a discrete-time system.

    Returns
    -------
    Wc : (N, N) ndarray, or float
        The Gramian. In continuous time with finite T it is integrated with Simpson's rule over steps of 0.001.
        With T = np.inf it solves the Lyapunov equation if the system is stable; if it is not, it prints a
        message and returns ``np.nan``. For any other `system`, including None, it returns None.
    """
    A_norm = _as_float(A_norm)
    n_nodes = A_norm.shape[0]
    B = np.eye(n_nodes)

    eigvals, _ = np.linalg.eig(A_norm)
    BB = B @ B.T

    # If time horizon is infinite, can only compute the Gramian when stable
    if T == np.inf:
        if system == "continuous":
            # If stable: solve using Lyapunov equation
            if np.max(np.real(eigvals)) < 0:
                return la.solve_continuous_lyapunov(A_norm, -BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.nan
        elif system == "discrete":
            # If stable: solve using Lyapunov equation
            if np.max(np.abs(eigvals)) < 1:
                return la.solve_discrete_lyapunov(A_norm, BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.nan
    # If time horizon is finite, perform numerical integration
    else:
        if system == "continuous":
            STEP = 0.001
            t = np.arange(0, (T + STEP / 2), STEP)
            # Accumulate e^{A t} over the steps, and the integrand e^{A t} B B^T e^{A^T t}
            dE = sp.linalg.expm(A_norm * STEP)
            dEa = np.zeros((n_nodes, n_nodes, len(t)))
            dEa[:, :, 0] = np.eye(n_nodes)
            dG = np.zeros((n_nodes, n_nodes, len(t)))
            dG[:, :, 0] = B @ B.T
            for i in np.arange(1, len(t)):
                dEa[:, :, i] = dEa[:, :, i - 1] @ dE
                dEab = dEa[:, :, i] @ B
                dG[:, :, i] = dEab @ dEab.T

            return sp.integrate.simpson(dG, x=t, dx=STEP, axis=2)
        elif system == "discrete":
            Ap = np.eye(n_nodes)
            Wc = np.eye(n_nodes)
            for _ in range(cast(int, T)):
                Ap = Ap @ A_norm
                Wc = Wc + Ap @ Ap.T

            return Wc
    return None  # system not recognised: returned None since 1.0, kept for compatibility


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
    x0 : (N,) or (N, 1) array_like
        Initial state. Boolean states are converted to 0/1 floats.
    xf : (N,) or (N, 1) array_like
        Target state. Boolean states are converted to 0/1 floats.

    Returns
    -------
    energy : (N, 1) ndarray
        Energy at each node; sum it for the total.
    """
    A_norm, B = _as_float(A_norm), _as_float(B)
    n_nodes = A_norm.shape[0]
    x0 = _column(_as_float(x0))
    xf = _column(_as_float(xf))

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

    delx = xf - eAT @ x0
    return np.multiply(np.linalg.pinv(G) @ delx, delx)
