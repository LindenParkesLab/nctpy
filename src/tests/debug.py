# Control code
import os, sys
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)
import numpy as np
from network_control.utils import matrix_normalization
from network_control.energies import optimal_input

# initialize inputs
np.random.seed(28)
A = np.random.rand(5,5)
A = matrix_normalization(A,c=1,version='continuous')
x0 = np.random.rand(5,1)
xf = np.random.rand(5,1)
B = np.eye(5)
S = np.eye(5)
T = 4
rho = 1

# compute optimal energy
x_opt,u_opt,n_err_opt = optimal_input(A,T,B,x0,xf,rho,S)

print(x_opt.shape)