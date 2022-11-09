import numpy as np
from numpy.linalg import eig
import scipy as sp
from scipy.linalg import schur
from numpy import matmul as mm
from scipy import sparse

np.random.seed(1)
# matrix norm
N = 50
b = np.random.randn(N, N)
b_symm = (b + b.T)/2
np.save('./fixtures/A.npy', b_symm)

# eigenvalue decomposition
A = b_symm
w, _ = eig(A)
l = np.abs(w).max()
# Matrix normalization for discrete-time systems
A_norm = A / (1 + l)
np.save('./fixtures/A_d_1.npy', A_norm)
A_norm = A_norm - np.eye(A.shape[0])
np.save('./fixtures/A_c_1.npy', A_norm)

# Matrix normalization for discrete-time systems
A_norm = A / (2 + l)
np.save('./fixtures/A_d_2.npy', A_norm)
A_norm = A_norm - np.eye(A.shape[0])
np.save('./fixtures/A_c_2.npy', A_norm)


# state and weight norm
x = np.random.randn(N, )
np.save('./fixtures/x.npy', x)
x_norm = x / np.linalg.norm(x, ord=2)
np.save('./fixtures/x_norm.npy', x_norm)
x_rank = sp.stats.rankdata(x)
# rescale to unit interval
x_scale = (x_rank - min(x_rank)) / (max(x_rank) - min(x_rank))
np.save('./fixtures/weights_rank_scale.npy', x_scale)
x_scale_const = x_scale + 1
np.save('./fixtures/weights_rank_scale_const.npy', x_scale_const)
x_scale = (x - min(x)) / (max(x) - min(x))
np.save('./fixtures/weights_scale.npy', x_scale)
x_scale_const = x_scale + 1
np.save('./fixtures/weights_scale_const.npy', x_scale_const)

# integrate U
T = 1000
U = np.random.randn(T, N)
np.save('./fixtures/u.npy', U)
energy = sp.integrate.simps(U.T ** 2)
np.save('./fixtures/u_int.npy', energy)

# avg control
# discrete
A_norm = A / (1 + l)
T, U = schur(A_norm, 'real')  # Schur stability
midMat = np.multiply(U, U).transpose()
v = np.matrix(np.diag(T)).transpose()
N = A_norm.shape[0]
P = np.diag(1 - np.matmul(v, v.transpose()))
P = np.tile(P.reshape([N, 1]), (1, N))
ac = sum(np.divide(midMat, P))
np.save('./fixtures/ac_d.npy', ac)
# continuous
A_norm = A_norm - np.eye(A.shape[0])
# System Size
n_nodes = A_norm.shape[0]
B = np.eye(n_nodes)
u, v = eig(A_norm)
BB = mm(B, np.transpose(B))
# If time horizon is infinite, can only compute the Gramian when stable
T=1
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
ac = G.diagonal()
np.save('./fixtures/ac_c.npy', ac)

# get control inputs
xf = np.random.rand(N,)
np.save('./fixtures/xf.npy', xf)
B = np.eye(N)
B[2, 2] = 0
B[12, 12] = 0
B[20, 20] = 0
B[4, 4] = 0
B[15, 15] = 0
np.save('./fixtures/B.npy', B)
def control(Anorm, T, B, x0, xf, system, xr='zero', rho=1, S='identity'):
    n_nodes = Anorm.shape[0]

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

    if system == 'continuous':
        # Set parameters
        dt = 0.001

        # Define joint state-costate matrix
        M = np.concatenate((np.concatenate((Anorm, np.dot(-B, B.T) / (2 * rho)), axis=1),
                            np.concatenate((-2 * S, -Anorm.T), axis=1)), axis=0)

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
                M[np.ix_(z + (i - 0) * n_nodes, z + (i - 1) * n_nodes)] = -Anorm

        # Constraints for the costate equations
        for i in np.arange(1, T):
            M[np.ix_(z + (i - 1 + T) * n_nodes, z + (i - 1) * n_nodes)] = -D
            M[np.ix_(z + (i - 1 + T) * n_nodes, z + (T - 2 + i) * n_nodes)] = I
            M[np.ix_(z + (i - 1 + T) * n_nodes, z + (T - 1 + i) * n_nodes)] = -Anorm.T

        # Construct boundary condition vector
        b = np.concatenate((np.dot(Anorm, x0), np.zeros((n_nodes * (T - 2), 1)), -xf, np.zeros((n_nodes * (T - 1), 1))),
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
        err_traj = np.linalg.norm(x[:, 1:] - (np.dot(Anorm, x[:, 0:-1]) + np.dot(B, u)))
        err = [err_system, err_traj]

    return x.T, u.T, err


# defaults (discrete)
x0=x
x,u,err = control(A / (1 + l), 2, np.eye(N), x0, xf, system='discrete', xr='zero', rho=1, S='identity')
np.save('./fixtures/control_discrete.npy', x,u,err)
# T
x,u,err = control(A / (1 + l), 7, np.eye(N), x0, xf, system='discrete', xr='zero', rho=1, S='identity')
np.save('./fixtures/control_T.npy', x,u,err)
# B
x,u,err = control(A_norm, 2, B, x0, xf, system='continuous', xr='zero', rho=1, S='identity')
np.save('./fixtures/control_B.npy', x,u,err)
# rho
x,u,err = control(A / (1 + l), 2, np.eye(N), x0, xf, system='discrete', xr='zero', rho=100, S='identity')
np.save('./fixtures/control_rho.npy', x,u,err)
# S
x,u,err = control(A / (1 + l), 2, np.eye(N), x0, xf, system='discrete', xr='zero', rho=1, S=B)
np.save('./fixtures/control_S.npy', x,u,err)
# system
x,u,err = control(A_norm, 2, np.eye(N), x0, xf, system='continuous', xr='zero', rho=100, S='identity')
np.save('./fixtures/control_continuous.npy', x,u,err)
# binary x
x_bin = np.random.randint(0,1,size=(N,1))
np.save('./fixtures/x_bin.npy', x,u,err)

