import numpy as np 

def openLoopControl( A, B, xi, U):
    ''' 
    This function caclulates the trajectory for the network given our model
     if there are no constraints, and the target state is unknown, using the
     control equation precess dx = Ax(t) + BU(t). x(t) is the state vector, A is
     the adjacency matrix, U(t) is the time varying input as specified by the
     user, and B selects the control set (stimulating electrodes)
    
       Inputs
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
    
       Outputs
     x             : x is the NxMxT trajectory (numpy array) that results from simulating
                   x(t+1) = Ax(t) + Bu(t) the equation with the parameters
                   above.
    
     @author JStiso 
     June 2017
    '''

    # Simulate trajectory
    T = np.size(U,2)
    N = np.size(A,0)
    M = np.size(xi,0)

    # initialize x
    x = np.zeros((N, M, T))
    xt = xi
    for t in range(T):
        x[:,:,t] = xt
        dx = np.matmul(A,xt) + np.matmul(B,np.squeeze(U[:,:,t])) # state equation
        xt = xt + dx
    return x