import numpy as np 
import scipy as sp
from scipy import sparse
import scipy.integrate
import scipy.linalg as la
from scipy.linalg import eig
from numpy import matmul as mm
from scipy.linalg import expm as expm
from numpy import transpose as tp


def get_control_inputs(A_norm, T, B, x0, xf, system, xr='zero', rho=1, S='identity'):
    """

    Args:
        A_norm: (NxN numpy array) Normalized structural connectivity matrix
        T: (float) Time horizon: how long you want to control for. Too large will give large error,
            too short will not give enough time for control
        B: (NxN numpy array) Input matrix: selects which nodes to put input into. Define
            so there is a 1 on the diagonal of elements you want to add input to, and 0 otherwise
        x0: (Nx1, numpy array) Initial state
        xf: (Nx1, numpy array) Target state
        system: (str) Time system: options, 'continuous' or 'discrete'
        xr: (Nx1, numpy array) Reference state
        rho: (float) Mixing parameter: determines the extent to which x is constrained alongside u
        S: (NxN numpy array) Constraint matrix: determines which nodes for which u (and x) will be constrained

    Returns:
        x: (t x N numpy array) State trajectory (neural activity)
        u: (t x N numpy array) Control inputs
        err: (1 x 2 numpy array) Numerical error
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

    if xr == 'x0':
        xr = x0
    elif xr == 'xf':
        xr = xf
    elif xr == 'zero':
        xr = np.zeros((n_nodes, 1))

    if S == 'identity':
        S = np.eye(n_nodes)

    if system == 'continuous':
        # Set parameters
        dt = 0.001

        # Define joint state-costate matrix
        M = np.concatenate((np.concatenate((A_norm, np.dot(-B, B.T) / (2 * rho)), axis=1),
                            np.concatenate((-2 * S, -A_norm.T), axis=1)), axis=0)

        # Define constant vector due to cost deviation from reference state
        c = np.concatenate((np.zeros((n_nodes, 1)), 2 * np.dot(S, xr)), axis=0)
        c = np.linalg.solve(M, c)

        # Compute matrix exponential and decompose into NxN blocks
        E = sp.linalg.expm(M * T)
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
        Ad = sp.linalg.expm(M * dt)
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
        err_xf = np.linalg.norm(x[:, -1] - xf)
        err = [err_costate, err_xf]
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


def minimum_energy_fast(A, T, B, x0_mat, xf_mat):
    """ This function computes the minimum energy required to transition between all pairs of brain states
    encoded in (x0_mat,xf_mat)

     Args:
      A: numpy array (N x N)
            System adjacency matrix
      B: numpy array (N x N)
            Control input matrix
      x0_mat: numpy array (N x n_transitions)
             Initial states (see expand_states)
      xf_mat: numpy array (N x n_transitions)
            Final states (see expand_states)
      T: float (1 x 1)
           Control horizon

    Returns:
      E: numpy array (N x n_transitions)
            Regional energy for all state transition pairs.
            Notes,
                np.sum(E, axis=0)
                    collapse over regions to yield energy associated with all transitions.
                np.sum(E, axis=0).reshape(n_states, n_states)
                    collapse over regions and reshape into a state by state transition matrix.
    """
    if type(x0_mat[0][0]) == np.bool_:
        x0_mat = x0_mat.astype(float)
    if type(xf_mat[0][0]) == np.bool_:
        xf_mat = xf_mat.astype(float)

    G = gramian(A, B, T, version='continuous')
    delx = xf_mat - np.matmul(expm(A*T), x0_mat)
    E = np.multiply(np.linalg.solve(G, delx), delx)

    return E


def integrate_u(U):
    """ This function integrates over some input squared to calculate energy using Simpson's integration.

    If your control set (B) is the identity this will likely give energies that are nearly identical to those calculated using a Reimann sum.
    However, when control sets are sparse inputs can be super curvy, so this method will be a bit more accurate.
     Args:
      U: numpy array (N x T)
            Input to the system (likely the output from minimum_input or optimal_input)
      
    Returns:
      energy: numpy array (N x 1)
            energy input into each node
    """

    if sp.__version__ < '1.6.0':
        energy = sp.integrate.simps(U.T**2)
    else:
        energy = sp.integrate.simpson(U.T**2)
    return energy


def gramian(A, B, T, system=None):
    """
    This function computes the controllability Gramian.
    Args:
        A:             np.array (n x n)
        B:             np.array (n x k)
        T:             np.array (1 x 1)
        system:       str
            options: 'continuous' or 'discrete'. default=None
    Returns:
        Wc:            np.array (n x n)
    """

    # System Size
    n_parcels = A.shape[0]

    u, v = eig(A)
    BB = mm(B, np.transpose(B))
    n = A.shape[0]

    # If time horizon is infinite, can only compute the Gramian when stable
    if T == np.inf:
        # check system
        if system == 'continuous':
            # If stable: solve using Lyapunov equation
            if np.max(np.real(u)) < 0:
                return la.solve_continuous_lyapunov(A,-BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.nan
        elif system == 'discrete':
            # If stable: solve using Lyapunov equation
            if np.max(np.abs(u)) < 1:
                return la.solve_discrete_lyapunov(A,BB)
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
            dE = sp.linalg.expm(A * STEP)
            dEa = np.zeros((n_parcels, n_parcels, len(t)))
            dEa[:, :, 0] = np.eye(n_parcels)
            # Collect Gramian difference
            dG = np.zeros((n_parcels, n_parcels, len(t)))
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
            Ap = np.eye(n)
            Wc = np.eye(n)
            for i in range(T):
                Ap = mm(Ap, A)
                Wc = Wc + mm(Ap, tp(Ap))

            return Wc


class ComputeControlEnergy():
    def __init__(self, A, control_tasks, system='continuous', c=1, cost='optimal', T=1):
        self.A = A
        self.control_tasks = control_tasks

        self.system = system
        self.c = c
        self.cost = cost
        self.T = T

    def _check_inputs(self):
        if self.A.shape[0] == self.A.shape[1]:
            self.n_nodes = self.A.shape[0]
        elif self.A.shape[0] != self.A.shape[1]:
            raise Exception("A matrix is not square. This routine requires A.shape[0] == A.shape[1]")

        try:
            A_norm = self.A_norm
        except AttributeError:
            self.A_norm = matrix_normalization(self.A, system=self.system, c=self.c)

    def run(self):
        self._check_inputs()

        E = []

        for control_task in tqdm(self.control_tasks):
            if self.cost == 'optimal':
                _, u, _ = optimal_input(A=self.A_norm, T=self.T, B=control_task['B'],
                                            x0=control_task['x0'], xf=control_task['xf'],
                                            rho=control_task['rho'], S=control_task['S'])
            elif self.cost == 'minimum':
                _, u, _ = minimum_input(A=self.A_norm, T=self.T, B=control_task['B'],
                                            x0=control_task['x0'], xf=control_task['xf'])

            # get energy
            e = integrate_u(u)
            E.append(np.sum(e))

        # store outputs as array
        self.E = np.array(E)


class ComputeOptimizedControlEnergy():
    def __init__(self, A, control_task, system='continuous', c=1, cost='optimal', T=1, n_steps=2, lr=0.01):
        self.A = A
        self.control_task = control_task

        self.system = system
        self.c = c
        self.cost = cost
        self.T = T

        self.n_steps = n_steps
        self.lr = lr

    def _check_inputs(self):
        if self.A.shape[0] == self.A.shape[1]:
            self.n_nodes = self.A.shape[0]
        elif self.A.shape[0] != self.A.shape[1]:
            raise Exception("A matrix is not square. This routine requires A.shape[0] == A.shape[1]")

        try:
            A_norm = self.A_norm
        except AttributeError:
            self.A_norm = matrix_normalization(self.A, system=self.system, c=self.c)

    def _get_energy(self, B):
        if self.cost == 'optimal':
            _, u, _ = optimal_input(A=self.A_norm, T=self.T, B=B,
                                        x0=self.control_task['x0'], xf=self.control_task['xf'],
                                        rho=self.control_task['rho'], S=self.control_task['S'])
        elif self.cost == 'minimum':
            _, u, _ = minimum_input(A=self.A_norm, T=self.T, B=B,
                                        x0=self.control_task['x0'], xf=self.control_task['xf'])

        E = integrate_u(u)
        E = np.sum(E)

        return E

    def _get_energy_perturbed(self, B):
        # container for perturbed energies
        E_p = np.zeros(self.n_nodes)

        # backup B
        B_bak = B.copy()

        for i in tqdm(np.arange(self.n_nodes)):
            # get B
            B_p = B_bak.copy()

            # add arbitrary amount of additional control to node i
            B_p[i, i] += 0.1

            # get the state trajectory (x_p) and the control inputs (u_p)
            E_p[i] = self._get_energy(B=B_p)

        return E_p

    def run(self):
        self._check_inputs()

        E_opt = np.zeros(self.n_steps)
        B_opt = np.zeros((self.n_steps, self.n_nodes))

        B_I = np.eye(self.n_nodes)  # B = identity

        for i in np.arange(self.n_steps):
            print('Running gradient step {0}'.format(i))

            if i == 0:
                B = B_I  # use identity
            else:
                B = np.zeros((self.n_nodes, self.n_nodes))
                B[np.diag_indices(self.n_nodes)] = B_opt[i-1, :]  # get optimized from previous step

            E = self._get_energy(B=B)  # get energy
            E_p = self._get_energy_perturbed(B=B)  # get perturbed energy
            E_d = E_p - E  # calculate energy delta

            B_o = np.zeros((self.n_nodes, self.n_nodes))  # initialize container for optimized weights
            B_o[np.diag_indices(self.n_nodes)] = B.diagonal() - (E_d * self.lr)  # step down gradient
            B_o = B_o / sp.linalg.norm(B_o) * sp.linalg.norm(B_I)  # normalize

            E_o = self._get_energy(B=B_o)  # get optimized energy

            E_opt[i] = E_o  # store optimized energy
            B_opt[i, :] = B_o.diagonal()  # retain optimized

        self.E_opt = E_opt
        self.B_opt = B_opt