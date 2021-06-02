import numpy as np 
import scipy as sp

def sim_state_eq( A, B, xi, U):
    """This function caclulates the trajectory for the network given our model
     if there are no constraints, and the target state is unknown, using the
     control equation precess dx = Ax(t) + BU(t). x(t) is the state vector, A is
     the adjacency matrix, U(t) is the time varying input as specified by the
     user, and B selects the control set (stimulating electrodes)
    
    Original matlab code by Jason Kim
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
     
     xi            : NxM initial state (numpy array) of your system where N is the number of
                   nodes, and M is the number of independent states (ie 
                   frequency bands). xi MUST have N rows, however, the number
                   of state measurements can change (and can be equal to 1). 
    
     U             : NxMxT matrix of Energy (numpy array), where N is the number of nodes, M
                   is the number of state measurements, and T is the number of
                   time points. For example, if you want to simulate the
                   trajectory resulting from stimulation, U could have
                   log(StimFreq)*StimAmp*StimDur as every element. You can
                   also enter U's that vary with time, or are different
                   accross frequency bands.
    
      Returns:
     x             : x is the NxMxT trajectory (numpy array) that results from simulating
                   x(t+1) = Ax(t) + Bu(t) the equation with the parameters
                   above.
    
     @author JStiso 
     June 2017
    """

    # Simulate trajectory
    T = np.size(U,2)
    N = np.size(A,0)
    M = np.size(xi,1)

    # initialize x
    x = np.zeros((N, M, T))
    xt = xi
    for t in range(T):
        x[:,:,t] = xt
        dx = np.matmul(A,xt) + np.matmul(B,(U[:,:,t]) )# state equation
        xt = xt + dx
    return x

def optimal_energy(A, T, B, x0, xf, rho, S):
    """This is a python adaptation of matlab code originally written by Tomaso Menara and Jason Kim
     compute optimal inputs/trajectories for a system to transition between two states
     Fabio, Tommy September 2017

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
    t = np.arange(0,(T+STEP),STEP)

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

def minimum_energy(A, T, B, x0, xf):
    """ This function computes the minimum energy required to transition between two states
     
     This is a python adaptation of code originally written by Jason Kim
     
     Computes minimum control energy for state transition.
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

    # Compute Matrix Exponential
    AT = np.concatenate((np.concatenate((A, -.5*(B.dot(B.T))), axis=1), 
                         np.concatenate((np.zeros(np.shape(A)), -A.T), axis=1)), axis=0)

    E = sp.linalg.expm(AT*T)

    # Compute Costate Initial Condition
    E12 = E[0:n,n:]
    E11 = E[0:n,0:n]
    p0 = np.linalg.pinv(E12).dot(xf - E11.dot(x0))

    # Compute Costate Initial Condition Error Induced by Inverse
    n_err = np.linalg.norm(E12.dot(p0) - (xf - E11.dot(x0)))

    # Prepare Simulation
    nStep=1000
    t = np.linspace(0,T,nStep+1)

    v0 = np.concatenate((x0, p0), axis=0)          # Initial Condition
    v = np.zeros((2*n,len(t)))          # Trajectory
    Et = sp.linalg.expm(AT*T/(len(t)-1))
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
