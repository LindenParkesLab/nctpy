import numpy as np 
import scipy as sp
from scipy import sparse
import scipy.integrate
import scipy.linalg as la
from numpy.linalg import eig
from numpy import matmul as mm
from scipy.linalg import expm as expm
from numpy import transpose as tp

from nctpy.utils import expm

def sim_state_eq(A_norm, B, x0, U, system=None):
    """This function calculate the trajectory for the system given our model if there are no constraints,
    and the target state is unknown.

    Args:
        A_norm (NxN numpy array): Normalized structural connectivity matrix
        B (NxN, numpy array): Control node matrix. Diagonal entries designate which nodes are control nodes and how much
            influence those nodes have of system dynamics. For example, B=np.eye(A_norm.shape[0]) would set all nodes as
            controllers with equal weight (1). This is referred to as a uniform full control set.
        x0 (Nx1, numpy array): Initial state. The initial condition of the system.
        U (NxT numpy array): System inputs. T is the number of time points.
            For example, if you want to simulate the trajectory resulting from stimulation, U could have
            log(StimFreq)*StimAmp*StimDur as every element. You can also enter U's that vary with time
        system (str): Time system. options, 'continuous' or 'discrete'.

    Returns:
        x (NxT, numpy array): Trajectory. This is trajectory of neural activity that results from simulating the
            system with the above parameters

    """

    # state vectors to float if they're bools
    if type(x0[0]) == np.bool_:
        x0 = x0.astype(float)

    # check dimensions of states
    if x0.ndim == 1:
        x0 = x0.reshape(-1, 1)

    # Simulate trajectory
    T = np.size(U, 1)
    N = np.size(A_norm, 0)

    # initialize x
    x = np.zeros((N, T))
    xt = x0

    if system is None:
        raise Exception("Time system not specified. "
                        "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                        "(see matrix_normalization help).")
    elif system != 'continuous' and system != 'discrete':
        raise Exception("Incorrect system specification. "
                        "Please specify either 'system=discrete' or 'system=continuous'.")
    elif system == 'continuous':
        for t in range(T):
            x[:, t] = xt[:, 0]
            dt = np.matmul(A_norm, xt) + np.matmul(B, np.reshape(U[:, t], (N, 1)))  # state equation
            xt = dt + xt
    elif system == 'discrete':
        for t in range(T):
            x[:, t] = xt[:, 0]
            xt_1 = np.matmul(A_norm, xt) + np.matmul(B, np.reshape(U[:, t], (N, 1)))  # state equation
            xt = xt_1

    return x


def get_control_inputs(A_norm, T, B, x0, xf, system=None, xr='zero', rho=1, S='identity', expm_version='scipy'):
    """This function extracts the state trajectory (x) and the control signals (u) associated with a control task.

    Args:
        A_norm (NxN, numpy array): normalized structural connectivity matrix.
        T (float): time horizon. The amount of time the model will run for. Too long will yield a large error,
            too short will not give enough time for control (i.e., the state transition may not complete).
        B (NxN, numpy array): control node matrix. Diagonal entries designate which nodes are control nodes and how much
            influence those nodes have of system dynamics. For example, B=np.eye(A_norm.shape[0]) would set all nodes as
            controllers with equal weight (1). This is referred to as a uniform full control set.
        x0 (Nx1, numpy array): initial state. The initial condition of the system.
        xf (Nx1, numpy array): target state. The state that the system will be controlled toward and should arrive at
            at time T.
        system (str): string variable that designates whether A was normalized for a continuous-time system or a
            discrete-time system. options: 'continuous' or 'discrete'. default=None.
        xr (Nx1, numpy array): reference state. This state governs the constraints placed on the state trajectory.
            By default this will be a vector of zeros of length N. This will result in the state trajectory being
            penalized for departing too far from 0 activity. Note, this only applies is S contains nodes to constrain.
        rho (float): mixing parameter. Determines the extent to which the state trajectory is constrained alongside the
            control signals. rho=1 equals maximum constraint. Note, rho must be >0 and will be ignored if
            S=np.zeros((A_norm.shape[0], A_norm.shape[0])).
        S (NxN, numpy array): constraint matrix for state trajectory. Determines which nodes in the state trajectory
            will be constrained. By default, all nodes' neural activity will be constrained.

    Returns:
        x (txN, numpy array): state trajectory (neural activity). t will be equal to (T/0.001)+1.
        u (txN, numpy array): control signals. t will be equal to (T/0.001)+1.
        n_err (list): numerical error.

    """

    n_nodes = A_norm.shape[0]

    # state vectors to float if they're bools
    if type(x0[0]) == np.bool_:
        x0 = x0.astype(float)
    if type(xf[0]) == np.bool_:
        xf = xf.astype(float)

    # check dimensions of states
    if x0.ndim == 1:
        x0 = x0.reshape(-1, 1)
    if xf.ndim == 1:
        xf = xf.reshape(-1, 1)

    if type(xr) == str:
        if xr == 'x0':
            xr = x0
        elif xr == 'xf':
            xr = xf
        elif xr == 'zero':
            xr = np.zeros((n_nodes, 1))
    else:
        if xr.ndim == 1:
            xr = xr.reshape(-1, 1)

    if type(S) == str and S == 'identity':
        S = np.eye(n_nodes)

    if system is None:
        raise Exception("Time system not specified. "
                        "Please nominate whether you are normalizing A for a continuous-time or a discrete-time system "
                        "(see matrix_normalization help).")
    elif system != 'continuous' and system != 'discrete':
        raise Exception("Incorrect system specification. "
                        "Please specify either 'system=discrete' or 'system=continuous'.")
    elif system == 'continuous':
        # Set parameters
        dt = 0.001

        # Define joint state-costate matrix
        M = np.concatenate((np.concatenate((A_norm, np.dot(-B, B.T) / (2 * rho)), axis=1),
                            np.concatenate((-2 * S, -A_norm.T), axis=1)), axis=0)

        # Define constant vector due to cost deviation from reference state
        c = np.concatenate((np.zeros((n_nodes, 1)), 2 * np.dot(S, xr)), axis=0)
        c = np.linalg.solve(M, c)

        # Compute matrix exponential and decompose into NxN blocks
        if expm_version == 'scipy':
            E = sp.linalg.expm(M * T)
        elif expm_version == 'eig':
            E = expm(M * T)
        r = np.arange(n_nodes)
        E11 = E[r, :][:, r]
        E12 = E[r, :][:, r + n_nodes]

        # Solve for initial costate as a function of initial and final states
        l0 = np.linalg.solve(E12,
                             (xf - np.dot(E11, x0) - np.dot(np.concatenate((E11 - np.eye(n_nodes), E12), axis=1), c)))

        # Construct discretized matrices to numerically integrate
        z = np.zeros((2 * n_nodes, int(np.round(T / dt) + 1)))
        z[:, 0] = np.concatenate((x0, l0), axis=0).flatten()
        I = np.eye(2 * n_nodes)
        if expm_version == 'scipy':
            Ad = sp.linalg.expm(M * dt)
        elif expm_version == 'eig':
            Ad = expm(M * dt)
        Bd = np.dot((Ad - I), c)

        # Simulate the state-costate trajectory
        for i in np.arange(1, int(np.round(T / dt)) + 1):
            z[:, i] = np.dot(Ad, z[:, i - 1]) + Bd.flatten()

        # Extract state and input from the joint state-costate equation
        x = z[r, :]
        u = np.dot(-B.T, z[r + n_nodes, :]) / (2 * rho)

        # Collect error
        err_costate = np.linalg.norm(np.dot(E12, l0) -
                                     (xf - np.dot(E11, x0) -
                                      np.dot(np.concatenate((E11 - np.eye(n_nodes), E12), axis=1), c)))
        err_xf = np.linalg.norm(x[:, -1].reshape(-1, 1) - xf)
        err = [err_costate, err_xf]

        return x.T, u.T, err
    elif system == 'discrete':
        # Define joint state - costate matrix
        C = np.dot(-B, B.T) / (2 * rho)
        D = 2 * S
        I = np.eye(n_nodes)
        z = np.arange(n_nodes)

        # Construct system matrix
        M = np.zeros(((2 * T - 1) * n_nodes, (2 * T - 1) * n_nodes))

        # Constraints for the state equations
        for i in np.arange(T + 1):
            M[np.ix_(z + (i - 1) * n_nodes, z + (T - 2 + i) * n_nodes)] = -C
            if i != T:
                M[np.ix_(z + (i - 1) * n_nodes, z + (i - 1) * n_nodes)] = I
            if i != 0:
                M[np.ix_(z + (i - 0) * n_nodes, z + (i - 1) * n_nodes)] = -A_norm

        # Constraints for the costate equations
        for i in np.arange(1, T):
            M[np.ix_(z + (i - 1 + T) * n_nodes, z + (i - 1) * n_nodes)] = -D
            M[np.ix_(z + (i - 1 + T) * n_nodes, z + (T - 2 + i) * n_nodes)] = I
            M[np.ix_(z + (i - 1 + T) * n_nodes, z + (T - 1 + i) * n_nodes)] = -A_norm.T

        # Construct boundary condition vector
        b = np.concatenate((np.dot(A_norm, x0), np.zeros((n_nodes * (T - 2), 1)), -xf, np.zeros((n_nodes * (T - 1), 1))),
                           axis=0) - \
            np.concatenate((np.zeros((n_nodes * T, 1)), np.tile(2 * np.dot(S, xr), (T - 1, 1))), axis=0)

        # Solve simultaneous state and costate equations
        M = sparse.csc_matrix(M)
        b = sparse.csc_matrix(b)
        v = sparse.linalg.spsolve(M, b)
        # v = np.linalg.solve(M, b)
        V = v.reshape((n_nodes, int(len(v) / n_nodes)), order='F')
        x = np.concatenate((x0, V[:, :T - 1], xf), axis=1)
        u = np.dot(-B.T, V[:, T - 1:]) / (2 * rho)

        # Collect error
        v = sparse.csc_matrix(np.expand_dims(v, axis=1))
        err_system = np.dot(M, v) - b
        err_system = np.linalg.norm(err_system.todense())
        err_traj = np.linalg.norm(x[:, 1:] - (np.dot(A_norm, x[:, 0:-1]) + np.dot(B, u)))
        err = [err_system, err_traj]

        return x.T, u.T, err


def integrate_u(u):
    """This function integrates over some input squared to calculate energy using Simpson's integration.

    If your control set (B) is the identity this will likely give energies that are nearly identical to those calculated
    using a Reimann sum. However, when control sets are sparse inputs can be super curvy, so this method will be a bit
    more accurate.

    Args:
        u (txN, numpy array): Control signals input to the system.
            These can be obtained using get_control_inputs.
      
    Returns:
        energy ((N,) numpy array): Energy input into each node.

    """

    if sp.__version__ < '1.6.0':
        energy = sp.integrate.simps(u.T**2)
    else:
        energy = sp.integrate.simpson(u.T**2)
    return energy


def gramian(A_norm, T, system=None):
    """This function computes the controllability Gramian.

    Args:
        A_norm (NxN, numpy array): normalized structural connectivity matrix.
        T (float): time horizon.
        system (str): string variable that designates whether A was normalized for a continuous-time system or a
            discrete-time system. options: 'continuous' or 'discrete'. default=None.

    Returns:
        Wc (NxN, numpy array): gramian matrix.

    """

    # System Size
    n_nodes = A_norm.shape[0]
    B = np.eye(n_nodes)

    u, v = eig(A_norm)
    BB = mm(B, np.transpose(B))

    # If time horizon is infinite, can only compute the Gramian when stable
    if T == np.inf:
        # check system
        if system == 'continuous':
            # If stable: solve using Lyapunov equation
            if np.max(np.real(u)) < 0:
                return la.solve_continuous_lyapunov(A_norm, -BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.nan
        elif system == 'discrete':
            # If stable: solve using Lyapunov equation
            if np.max(np.abs(u)) < 1:
                return la.solve_discrete_lyapunov(A_norm, BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.nan
    # If time horizon is finite, perform numerical integration
    else:
        # check system
        if system == 'continuous':
            # Number of integration steps
            STEP = 0.001
            t = np.arange(0, (T+STEP/2), STEP)
            # Collect exponential difference
            dE = sp.linalg.expm(A_norm * STEP)
            dEa = np.zeros((n_nodes, n_nodes, len(t)))
            dEa[:, :, 0] = np.eye(n_nodes)
            # Collect Gramian difference
            dG = np.zeros((n_nodes, n_nodes, len(t)))
            dG[:, :, 0] = mm(B, B.T)
            for i in np.arange(1, len(t)):
                dEa[:, :, i] = mm(dEa[:, :, i-1], dE)
                dEab = mm(dEa[:, :, i], B)
                dG[:, :, i] = mm(dEab, dEab.T)

            # Integrate
            if sp.__version__ < '1.6.0':
                G = sp.integrate.simps(dG, t, STEP, 2)
            else:
                G = sp.integrate.simpson(dG, t, STEP, 2)

            return G
        elif system == 'discrete':
            Ap = np.eye(n_nodes)
            Wc = np.eye(n_nodes)
            for i in range(T):
                Ap = mm(Ap, A_norm)
                Wc = Wc + mm(Ap, tp(Ap))

            return Wc


def minimum_energy_fast(A_norm, T, B, x0, xf):
    # System Size
    n_nodes = A_norm.shape[0]

    try:
        if type(x0[0][0]) == np.bool_:
            x0 = x0.astype(float)
        if type(xf[0][0]) == np.bool_:
            xf = xf.astype(float)
    except:
        if type(x0[0]) == np.bool_:
            x0 = x0.astype(float)
        if type(xf[0]) == np.bool_:
            xf = xf.astype(float)

    if x0.ndim == 1:
        x0 = x0.reshape(-1, 1)
    if xf.ndim == 1:
        xf = xf.reshape(-1, 1)

    # Number of integration steps
    nt = 1000
    dt = T/nt

    # Numerical integration with Simpson's 1/3 rule
    # Integration step
    dE = sp.linalg.expm(A_norm * dt)
    # Accumulation of expm(A * dt)
    dEA = np.eye(n_nodes)
    # Gramian
    G = np.zeros((n_nodes, n_nodes))

    for i in np.arange(1, nt/2):
        # Add odd terms
        dEA = np.matmul(dEA, dE)
        p1 = np.matmul(dEA, B)
        # Add even terms
        dEA = np.matmul(dEA, dE)
        p2 = np.matmul(dEA, B)
        G = G + 4 * (np.matmul(p1, p1.transpose())) + 2 * (np.matmul(p2, p2.transpose()))

    # Add final odd term
    dEA = np.matmul(dEA, dE)
    p1 = np.matmul(dEA, B)
    G = G + 4 * (np.matmul(p1, p1.transpose()))

    # Divide by integration step
    E = sp.linalg.expm(A_norm * T)
    G = (G + np.matmul(B, B.transpose()) + np.matmul(np.matmul(E, B), np.matmul(E, B).transpose())) * dt / 3

    delx = xf - np.matmul(E, x0)
    E = np.multiply(np.matmul(np.linalg.pinv(G), delx), delx)

    return E
