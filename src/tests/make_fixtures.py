import numpy as np
from numpy.linalg import eig
import scipy as sp
from scipy.linalg import schur
from numpy import matmul as mm

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
