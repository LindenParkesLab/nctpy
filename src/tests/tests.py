# Control code
import os, sys
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)
from network_control.energies import minimum_input, optimal_input, optimal_input_gen, sim_state_eq, integrate_u, gramian, minimum_energy_fast
from network_control.utils import matrix_normalization

# Make it like MATLAB
# Numpy
import numpy as np 
from numpy.linalg import inv as inv
from numpy.linalg import solve as lasolve
from numpy import transpose as tp
from numpy import concatenate as cat
from numpy.random import random as rand
from numpy import identity as eye
from numpy import zeros as zeros
from numpy import matmul as mm
np.set_printoptions(threshold=sys.maxsize)
# Scipy
import scipy as sp
import scipy.integrate
import scipy.linalg as la
from scipy.linalg import eig
from scipy.linalg import svd
from scipy.linalg import expm as expm
# Performance and visualization
from time import perf_counter as tpf
from matplotlib.pyplot import plot as plot


# Define test cases
n = 100
T = 1
A = rand((n,n)) - 0.5
A = A
B = eye(n)
B[:,n-5:n] = B[:,n-5:n]*0
S = zeros((n,n))
x0 = rand((n,1)) - 0.5
xf = rand((n,1)) - 0.5

A = matrix_normalization(A,version='continuous')


# Optimal Hamiltonian
tic = tpf()
X_opt, U_opt, err_opt = optimal_input(A,T,B,x0,xf,100000,B)
Eopt = sum(integrate_u(U_opt))*.001

# Optimal Hamiltonian Gen
tic_opt = tpf()
X_optG, U_optG, err_optG = optimal_input_gen(A,T,B,x0,xf,100000,B)
EoptG = sum(integrate_u(U_optG))*.001

# Minimum Hamiltonian
tic_optG = tpf()
X_min, U_min, err_min = minimum_input(A,T,B,x0,xf)
Emin = sum(integrate_u(U_min))*.001
tic_min = tpf()

# Minimum Gramian
EminF = minimum_energy_fast(A,T,B,x0,xf)
EminF = sum(EminF).item()
tic_minf = tpf()

Eall = np.array([Eopt, EoptG, Emin, EminF])
print(np.c_[Eall] - Eall)
print(err_opt)
print(err_optG)
print(err_min)
# print([tic_opt-tic, tic_optG-tic_opt, tic_min-tic_opt, tic_minf-tic_min])
print(Eall)
# print(la.norm(X_optG[1000,:] - xf.T))
# print(U_min)