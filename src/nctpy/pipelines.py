import numpy as np
import scipy as sp
import scipy.linalg as la
from tqdm import tqdm

from nctpy.energies import get_control_inputs, integrate_u
from nctpy.utils import matrix_normalization

class ComputeControlEnergy():
    """This class will compute control energy for a set of independent control tasks stored in control_tasks.
    See /scripts/control_energy_wrapper.ipynb for an example.

    """
    def __init__(self, A, control_tasks, system=None, c=1, T=1):
        """Initializes inputs to ComputeControlEnergy

        Args:
            A (NxN, numpy array): adjacency matrix representing a structural connectome.
            control_tasks (list): list of control tasks wherein each entry is a dictionary variable.
                For example,
                    control_tasks = []  # list of tasks

                    control_task = dict()  # initialize dict
                    control_task['x0'] = x0  # store initial state
                    control_task['xf'] = xf  # store target state
                    control_task['B'] = B  # store control nodes
                    control_task['S'] = S  # store state trajectory constraints
                    control_task['rho'] = rho  # store rho
                    control_tasks.append(control_task)
                See get_control_inputs for information on above inputs.
            system (str): string variable that determines whether A is normalized for a continuous-time system or a
                discrete-time system. options: 'continuous' or 'discrete'. default=None.
            c (int): normalization constant, default=1.
            T (float): time horizon.

        """

        self.A = A
        self.control_tasks = control_tasks

        self.system = system
        self.c = c
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
            if control_task['rho'] == 0:
                S = np.eye(self.n_nodes)
            else:
                S = control_task['S']

            _, u, _ = get_control_inputs(A_norm=self.A_norm, T=self.T, B=control_task['B'],
                                         x0=control_task['x0'], xf=control_task['xf'], system=self.system,
                                         rho=control_task['rho'], S=S)

            # get energy
            e = integrate_u(u)
            E.append(np.sum(e))

        # store outputs as array
        self.E = np.array(E)


class ComputeOptimizedControlEnergy():
    """This class will compute `optimized` control energy for a single control task stored in control_task.
    See /scripts/path_a_control_energy_binary.ipynb for an example.

    """

    def __init__(self, A, control_task, system=None, c=1, T=1, n_steps=2, lr=0.01):
        """Initializes inputs to ComputeOptimizedControlEnergy

        Args:
            A (NxN, numpy array): adjacency matrix representing a structural connectome.
            control_task (dict): control task.
                For example,
                    control_task = dict()  # initialize dict
                    control_task['x0'] = x0  # store initial state
                    control_task['xf'] = xf  # store target state
                    control_task['B'] = B  # store control nodes
                    control_task['S'] = S  # store state trajectory constraints
                    control_task['rho'] = rho  # store rho
                See get_control_inputs for information on above inputs.
            system (str): string variable that determines whether A is normalized for a continuous-time system or a
                discrete-time system. options: 'continuous' or 'discrete'. default=None.
            c (int): normalization constant, default=1.
            T (float): time horizon.
            n_steps (int): number of steps down gradient that will be performed.
            lr (float): learning rate for gradient descent.

        """

        self.A = A
        self.control_task = control_task

        self.system = system
        self.c = c
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
        if self.control_task['rho'] == 0:
            S = np.eye(self.n_nodes)
        else:
            S = self.control_task['S']

        _, u, _ = get_control_inputs(A_norm=self.A_norm, T=self.T, B=B,
                                     x0=self.control_task['x0'], xf=self.control_task['xf'], system=self.system,
                                     rho=self.control_task['rho'], S=S)

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
