import numpy as np 
import scipy as sp
import scipy.linalg as la
from scipy.linalg import eig
from scipy.linalg import svd
from numpy import matmul as mm
from scipy.linalg import expm as expm
from numpy import transpose as tp
from numpy import concatenate as cat

def sim_state_eq( A, B, xi, U, version=None):
    """This function caclulates the trajectory for the network given our model
     if there are no constraints, and the target state is unknown, using the
     control equation precess x(t+1) = Ax(t) + BU(t). x(t) is the state vector, A is
     the adjacency matrix, U(t) is the time varying input as specified by the
     user, and B selects the control set (stimulating electrodes). This code assumes a DISCRETE system
    
    Args:
     A             : NxN state matrix (numpy array), where N is the number of nodes in your
                   network (for example, a structural connectivity matrix 
                   constructed from DTI). A should be stable to prevent
                   uncontrolled trajectories.
    
     B             : NxN input matrix (numpy array), where N is the number of nodes. B
                   selects where you want your input energy to be applied to.
                   For example, if B is the Identity matrix, then input energy
                   will be applied to all nodes in the network. If B is a
                   matrix of zeros, but B(1,1) = 1. then energy will only be
                   applied at the first node.
     
     xi            : Nx1 initial state (numpy array) of your system where N is the number of
                   nodes. xi MUST have N rows. 
    
     U             : NxT matrix of input (numpy array), where N is the number of nodes
                   and T is the number of
                   time points. For example, if you want to simulate the
                   trajectory resulting from stimulation, U could have
                   log(StimFreq)*StimAmp*StimDur as every element. You can
                   also enter U's that vary with time

     version        :options: 'continuous' or 'discrete' (str). default=None
                    string variable that determines whether A is a continuous-time system or a discrete-time
                    system
    
      Returns:
     x             : x is the NxT trajectory (numpy array) that results from simulating
                   x(t+1) = Ax(t) + Bu(t) the equation with the parameters
                   above.
    
     @author JStiso 
     June 2017
    """
    # check inouts
    if version == 'continuous':
            print("Simulating for a continuous-time system")
    elif version == 'discrete':
        print("Simulating for a discrete-time system")
    elif version == None:
        raise Exception("Time system not specified. "
                            "Please indicate whether you are simulating a continuous-time or a discrete-time system "
                            "(see matrix_normalization for help).")
    # state vectors to float if they're bools
    if type(xi[0]) == np.bool_:
        x0 = xi.astype(float)
    # check dimensions of states
    if xi.ndim == 1:
        xi = xi.reshape(-1, 1)

    # Simulate trajectory
    T = np.size(U,1)
    N = np.size(A,0)

    # initialize x
    x = np.zeros((N, T))
    xt = xi
    if version == 'discrete':
        for t in range(T):
            x[:,t] = np.reshape(xt, N) # annoying python 1d array thing
            xt_1 = np.matmul(A,xt) + np.matmul(B,np.reshape(U[:,t],(N,1) ))# state equation
            xt = xt_1
    elif version == 'continuous':
        for t in range(T):
            x[:,t] = np.reshape(xt, N) # annoying python 1d array thing
            dt = np.matmul(A,xt) + np.matmul(B,np.reshape(U[:,t],(N,1) ))# state equation
            xt = dt + xt
    return x

def optimal_input_gen(A, T, B, x0, xf, rho, S):
    """This is a python adaptation of matlab code originally written by Tomaso Menara and Jason Kim
     compute optimal inputs/trajectories for a system to transition between two states
     Fabio, Tommy September 2017. This code assumes a CONTINUOUS system

     Args:
     A: (NxN numpy array) Structural connectivity matrix
     B: (NxN numpy array) Input matrix: selects which nodes to put input into. Define
           so there is a 1 on the diagonal of elements you want to add input to, 
           and 0 otherwise 
     S: (NxN numpy array) Selects nodes whose distance you want to constrain, Define so
           that there is a 1 on the diagonal of elements you want to
           constrain, and a zero otherwise
     T: (float) Time horizon: how long you want to control for. Too large will give
           large error, too short will not give enough time for control
    rho: (float) weights energy and distance constraints. Small rho leads to larger
           energy
    
     Returns:
     X_opt: (TxN numpy array) 
          The optimal trajectory through state space
     U_opt: (TxN numpy array) 
          The optimal energy
     n_err: (float) 
          the error associated with this calculation. Errors will be larger when B is not identity, 
              and when A is large. Large T and rho will also tend to increase the error
    
    -------------- Change Log -------------
        JStiso April 2018
          Changed S to be an input, rather than something defined internally
        
        Jason Kim January 2021
          Changed the forward propagation of states to matrix exponential to
          avoid reliance on MATLAB toolboxes. Also changed definition of expanded
          input U to save time by avoiding having to resize the matrix.
          Also changed the initialization of U_opt for the same reason.
        
        JStiso 2021
            Translated to Python

        Jason Kim August 2021
          Different calculation of initial condition of costate, p(0), due to some convergence
          questions in the original implementation

    """

    n = np.shape(A)[1]

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

    Sbar = np.eye(n) - S
    np.shape(np.dot(-B,B.T)/(2*rho))

    Atilde = np.concatenate((np.concatenate((A, np.dot(-B,B.T)/(2*rho)), axis=1), 
                            np.concatenate((-2*S, -A.T), axis=1)), axis=0)

    M = sp.linalg.expm(Atilde*T)
    M11 = M[0:n,0:n]
    M12 = M[0:n,n:]
    M21 = M[n:,0:n]
    M22 = M[n:,n:]

    N = np.linalg.solve(Atilde,(M-np.eye(np.shape(Atilde)[0])))
    c = np.dot(np.dot(N,np.concatenate((np.zeros((n,n)),S),axis = 0)),2*xf)
    c1 = c[0:n]
    c2 = c[n:]

    p0 = np.linalg.solve(M12, xf - mm(M11,x0) - c1)
    n_err = np.linalg.norm(mm(M12,p0) - (xf - mm(M11,x0) - c1)) # norm(error)

    STEP = 0.001
    t = np.arange(0,(T+STEP/2),STEP)

    U = np.dot(np.ones((np.size(t),1)),2*xf.T)

    # Discretize continuous-time input for convolution
    Atilde_d = sp.linalg.expm(Atilde*STEP)
    Btilde_d = np.linalg.solve(Atilde,
                               np.dot((Atilde_d-np.eye(2*n)),np.concatenate((np.zeros((n,n)),S), axis=0)))

    # Propagate forward discretized model
    xp = np.zeros((2*n,np.size(t)))
    xp[:,0:1] = np.concatenate((x0,p0), axis=0)
    for i in np.arange(1,np.size(t)):
        xp[:,i] = np.dot(Atilde_d,xp[:,i-1]) + np.dot(Btilde_d,U[i-1,:].T)

    xp = xp.T

    U_opt = np.zeros((np.size(t),np.shape(B)[1]))
    for i in range(np.size(t)):
        U_opt[i,:] = -(1/(2*rho))*np.dot(B.T,xp[i,n:].T)

    X_opt = xp[:,0:n]
    
    return X_opt, U_opt, n_err

def optimal_input(A, T, B, x0, xf, rho, S):
    """This is a python adaptation of matlab code originally written by Tomaso Menara and Jason Kim
     compute optimal inputs/trajectories for a system to transition between two states
     Fabio, Tommy September 2017. This code assumes a CONTINUOUS system

     Args:
     A: (NxN numpy array) Structural connectivity matrix
     B: (NxN numpy array) Input matrix: selects which nodes to put input into. Define
           so there is a 1 on the diagonal of elements you want to add input to, 
           and 0 otherwise 
     S: (NxN numpy array) Selects nodes whose distance you want to constrain, Define so
           that there is a 1 on the diagonal of elements you want to
           constrain, and a zero otherwise
     T: (float) Time horizon: how long you want to control for. Too large will give
           large error, too short will not give enough time for control
    rho: (float) weights energy and distance constraints. Small rho leads to larger
           energy
    
     Returns:
     X_opt: (TxN numpy array) 
          The optimal trajectory through state space
     U_opt: (TxN numpy array) 
          The optimal energy
     n_err: (float) 
          the error associated with this calculation. Errors will be larger when B is not identity, 
              and when A is large. Large T and rho will also tend to increase the error
    
    -------------- Change Log -------------
        JStiso April 2018
          Changed S to be an input, rather than something defined internally
        
        Jason Kim January 2021
          Changed the forward propagation of states to matrix exponential to
          avoid reliance on MATLAB toolboxes. Also changed definition of expanded
          input U to save time by avoiding having to resize the matrix.
          Also changed the initialization of U_opt for the same reason.
        
        JStiso 2021
            Translated to Python

    """

    n = np.shape(A)[1]

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

    Sbar = np.eye(n) - S
    np.shape(np.dot(-B,B.T)/(2*rho))

    Atilde = np.concatenate((np.concatenate((A, np.dot(-B,B.T)/(2*rho)), axis=1), 
                            np.concatenate((-2*S, -A.T), axis=1)), axis=0)

    M = sp.linalg.expm(Atilde*T)
    M11 = M[0:n,0:n]
    M12 = M[0:n,n:]
    M21 = M[n:,0:n]
    M22 = M[n:,n:]

    N = np.linalg.solve(Atilde,(M-np.eye(np.shape(Atilde)[0])))
    c = np.dot(np.dot(N,np.concatenate((np.zeros((n,n)),S),axis = 0)),2*xf)
    c1 = c[0:n]
    c2 = c[n:]

    p0 = np.dot(np.linalg.pinv(np.concatenate((np.dot(S,M12),np.dot(Sbar,M22)), axis = 0)),
                        (-np.dot(np.concatenate((np.dot(S,M11),np.dot(Sbar,M21)),axis=0),x0) - 
                         np.concatenate((np.dot(S,c1),np.dot(Sbar,c2)), axis=0) + 
                         np.concatenate((np.dot(S,xf),np.zeros((n,1))), axis=0)))
    
    n_err = np.linalg.norm(np.dot(np.concatenate((np.dot(S,M12),np.dot(Sbar,M22)), axis = 0),p0) - 
                           (-np.dot(np.concatenate((np.dot(S,M11),np.dot(Sbar,M21)),axis=0),x0) - 
                            np.concatenate((np.dot(S,c1),np.dot(Sbar,c2)), axis=0) + 
                            np.concatenate((np.dot(S,xf),np.zeros((n,1))), axis=0))) # norm(error)

    STEP = 0.001
    t = np.arange(0,(T+STEP/2),STEP)

    U = np.dot(np.ones((np.size(t),1)),2*xf.T)

    # Discretize continuous-time input for convolution
    Atilde_d = sp.linalg.expm(Atilde*STEP)
    Btilde_d = np.linalg.solve(Atilde,
                               np.dot((Atilde_d-np.eye(2*n)),np.concatenate((np.zeros((n,n)),S), axis=0)))

    # Propagate forward discretized model
    xp = np.zeros((2*n,np.size(t)))
    xp[:,0:1] = np.concatenate((x0,p0), axis=0)
    for i in np.arange(1,np.size(t)):
        xp[:,i] = np.dot(Atilde_d,xp[:,i-1]) + np.dot(Btilde_d,U[i-1,:].T)

    xp = xp.T

    U_opt = np.zeros((np.size(t),np.shape(B)[1]))
    for i in range(np.size(t)):
        U_opt[i,:] = -(1/(2*rho))*np.dot(B.T,xp[i,n:].T)

    X_opt = xp[:,0:n]
    
    return X_opt, U_opt, n_err

def minimum_input(A, T, B, x0, xf):
    """ This function computes the minimum input required to transition between two states
     
     This is a python adaptation of code originally written by Jason Kim
     
     Computes minimum control input for state transition.

     This code assumes a CONTINUOUS system
     Args:
      A: numpy array (N x N)
            System adjacency matrix 
      B: numpy array (N x N)
            Control input matrix 
      x0: numpy array (N x t)
             Initial state 
      xf: numpy array (N x t)
            Final state 
      T: float (1 x 1) 
           Control horizon 
      
    Returns:
      x: numpy array (N x t)
            State Trajectory 
      u: numpy array (N x t)
           Control Input 
    """

    # System Size
    n = np.shape(A)[0]

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
    
    # Compute Matrix Exponential
    AT = np.concatenate((np.concatenate((A, -.5*(B.dot(B.T))), axis=1), 
                         np.concatenate((np.zeros(np.shape(A)), -A.T), axis=1)), axis=0)

    E = sp.linalg.expm(AT*T)

    # Compute Costate Initial Condition
    E12 = E[0:n,n:]
    E11 = E[0:n,0:n]
    p0 = la.solve(E12,xf - E11.dot(x0))

    # Compute Costate Initial Condition Error Induced by Inverse
    n_err = np.linalg.norm(E12.dot(p0) - (xf - E11.dot(x0)))

    # Prepare Simulation
    STEP = 0.001
    t = np.arange(0,(T+STEP/2),STEP)

    v0 = np.concatenate((x0, p0), axis=0)          # Initial Condition
    v = np.zeros((2*n,len(t)))          # Trajectory
    Et = sp.linalg.expm(AT*STEP)
    v[:,0] = v0.T

    # Simulate State and Costate Trajectories
    for i in np.arange(1,len(t)):
        v[:,i] = Et.dot(v[:,i-1])

    x = v[0:n,:]
    u = -0.5*B.T.dot(v[np.arange(0,n)+n,:])

    # transpose to be similar to opt_eng_cont
    u = u.T
    x = x.T

    return x, u, n_err


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

    G = gramian(A,B,T,version='continuous')
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



def gramian(A,B,T,version=None,tol=1e-12):
    """
    This function computes the controllability Gramian.
    Args:
        A:             np.array (n x n)
        B:             np.array (n x k)
        T:             np.array (1 x 1)
        version:       str
            options: 'continuous' or 'discrete'. default=None
    Returns:
        Wc:            np.array (n x n)
    """

    # System Size
    n_parcels = A.shape[0]

    u,v = eig(A)
    BB = mm(B,np.transpose(B))
    n = A.shape[0]
    
    
    # If time horizon is infinite, can only compute the Gramian when stable
    if T == np.inf:
        # check version
        if version=='continuous':
            # If stable: solve using Lyapunov equation
            if(np.max(np.real(u)) < 0):
                return la.solve_continuous_lyapunov(A,-BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.NAN
        elif version=='discrete':
            # If stable: solve using Lyapunov equation
            if(np.max(np.abs(u)) < 1):
                return la.solve_discrete_lyapunov(A,BB)
            else:
                print("cannot compute infinite-time Gramian for an unstable system!")
                return np.NAN

            
    # If time horizon is finite, perform numerical integration
    else:
        # check version
        if version=='continuous':
            ## Compute required number of steps for desired precision
            # Prefactors and matrix powers
            # M1 = mm((mm(expm(A*0),B)), tp(mm(expm(A*0),B)))
            # M2 = mm((mm(expm(A*T),B)), tp(mm(expm(A*T),B)))
            # A2 = mm(A,A)
            # A3 = mm(A2,A)
            # A4 = mm(A3,A)
            # # Fourth derivative at start and end
            # D1 = mm(A4,M1) + 4*mm(mm(A3,M1),tp(A)) + 6*mm(mm(A2,M1),tp(A2)) + 4*mm(mm(A,M1),tp(A3)) + mm(M1,tp(A4))
            # D2 = mm(A4,M2) + 4*mm(mm(A3,M2),tp(A)) + 6*mm(mm(A2,M2),tp(A2)) + 4*mm(mm(A,M2),tp(A3)) + mm(M2,tp(A4))
            # # Get maximum error
            # u1,s1,v1 = svd(D1)
            # u2,s2,v2 = svd(D2)
            # mmax = np.max([s1,s2])
            # n = pow(pow(T,5)*mmax/(180*tol),1/4)
            # n = np.int(np.ceil(n))
            # print(n)
            
            # Number of integration steps
            STEP = 0.001
            t = np.arange(0,(T+STEP/2),STEP)
            # Collect exponential difference
            dE = sp.linalg.expm(A * STEP)
            dEa = np.zeros((n_parcels,n_parcels,len(t)))
            dEa[:,:,0] = np.eye(n_parcels)
            # Collect Gramian difference
            dG = np.zeros((n_parcels,n_parcels,len(t)))
            dG[:,:,0] = mm(B,B.T)
            for i in np.arange(1, len(t)):
                dEa[:,:,i] = mm(dEa[:,:,i-1],dE)
                dEab = mm(dEa[:,:,i],B)
                dG[:,:,i] = mm(dEab,dEab.T)

            # Integrate
            if sp.__version__ < '1.6.0':
                G = sp.integrate.simps(dG,t,STEP,2)
            else:
                G = sp.integrate.simpson(dG,t,STEP,2)
            return G
        elif version=='discrete':
            Ap = np.eye(n)
            Wc = np.eye(n)
            for i in range(T):
                Ap = mm(Ap,A)
                Wc = Wc + mm(Ap,tp(Ap))
            return Wc
                