from tkinter import E
import numpy as np
from scipy.linalg import eig

def T(R):
    """Generate the transition matrix from rate parameters
    
    Parameters
    ----------
    R : ndarray (5,)
        - Rate and Interaction Parameters
        - [0] index: p / stayOff
        - [1] index: q / stayOn
        - [3] index: r / off
        - [4] index: s / on

    Returns
    -------
    TM : NDArray (4,4)
        - Transition Matrix as defined in the model

    Raises
    ------
    N/A
    """
    # Generate matrix from rates
    TM = np.array(
        [[R[0]**2,R[0]*R[2],R[0]*R[2],R[2]**2],
        [R[0]*R[3],R[0]*R[1],R[4]*R[2]*R[3],R[4]*R[1]*R[2]],
        [R[0]*R[3],R[4]*R[2]*R[3],R[0]*R[1],R[4]*R[1]*R[2]],
        [R[3]*R[3],R[4]*R[1]*R[3],R[4]*R[1]*R[3],(R[4]*R[1])**2]])
    # Normalize each column to 1 (for probability)
    for i in range(4):
        TM[:,i] /= TM[:,i].sum()
    return TM

def randomer(n:int):
    """Generate random set of rate parameters R and interaction parameter \lambda

    Parameters
    ----------
    n : int
        - No. of sets to generate

    Returns
    -------
    R : NDArray (5,)
        - Rate and Interaction Parameters
        - [0] index: p / stayOff
        - [1] index: q / stayOn
        - [3] index: r / off
        - [4] index: s / on

    Raises
    ------
    N/A
    """
    # Generates uniform random sets
    R = np.random.rand(n,5)
    # Converts \lambda to uniform over log scale 
    # R[:,4] = np.power(10,4*(R[:,4]-0.5))
    R[:,4] = 0.01 + 99.99*R[:,4] 
    return R

def markeig(TM,eps:float=0.001):
    """Finds & manipulates the eigenvector that satisfies the Markov property
    
    Parameters
    ----------
    TM : NDArray (4, 4)
        - Transition matrix
    eps : float, optional
        - To account for precision error, by default 0.001

    Returns
    -------
    evec : NDArray (4,)
        - Eigenvector equivalent to steady state probabilities

    Raises
    ------
    Eigenvector not found
    """
    # Find all eigenvalues and eigenvectors of TM
    eval,evec = eig(TM) 
    # Magnitude of sum of sign is 4 if all are same. Except 0 prob??
    sign = abs(np.sum(np.sign(evec),axis=0))==4
    # Check if eigenvalue is 1
    eval1=(abs(eval[sign]-1)<eps) 
    if (np.sum(sign) != 1 and eval1):
        raise Exception("Proper Eigenvector not found")
    # Choose eigenvector with elements of same sign
    evec = evec[:,sign]
    # Normalize eigenvector
    evec = evec/np.sum(evec) 
    # Convert to 1D
    evec = np.concatenate(evec)         
    return evec

def statechange(TM,s:int):
    """Advances the Markov chain to state at next time step
    
    Parameters
    ----------
    TM : NDArray (4, 4)
        - Transition matrix
    s : int
        - Current state number [@ t_i)]

    Returns
    -------
    s : int
        - Next state number [@t_{i+1}]

    Raises
    ------
    N/A
    """
    # Define state space
    states=np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    # Get current state vector
    S=states[s]
    # Advance markov chain to next step
    S=TM@S 
    # Collapse to particular state number depending on probability
    s=np.random.choice(4,p=S)
    return s

def timeseries(tmax:int,TM):
    """Generates a time-series of states for the provided Markov chain 

    Parameters
    ----------
    tmax : int
        - No. of timestep iterations to make
    TM : NDArray (4, 4)
        - Transition matrix

    Returns
    -------
    TS : NDArray (2, tmax)
        - Time-series data

    Raises
    ------
    N/A
    """
    # Get a random initial state
    s=np.random.randint(4)
    TS=np.array([[0,s]])
    for t in range(tmax):
        s=statechange(TM,s)
        TS=np.append(TS,[[t+1,s]],axis=0)
    return TS