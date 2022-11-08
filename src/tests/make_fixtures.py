import numpy as np
from numpy.linalg import eig
import scipy as sp

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
