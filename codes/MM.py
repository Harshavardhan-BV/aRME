import numpy as np
from scipy.linalg import eig

def T(R):
    """Generate the transition matrix from rate parameters
    
    Parameters
    ----------
    R : ndarray (6,)
        - Rate and Interaction Parameters
        - [0] index: p / stayOff
        - [1] index: q / stayOn
        - [2] index: r / off
        - [3] index: s / on
        - [4] index: l / coordination parameter (lambda)
        - [5] index: d / competition parameter (delta)


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
        [[R[0]**2,R[0]*R[2],R[0]*R[2],0],
        [R[0]*R[3],R[0]*R[1],0,R[5]*R[1]*R[2]],
        [R[0]*R[3],0,R[0]*R[1],R[5]*R[1]*R[2]],
        [0,R[4]*R[1]*R[3],R[4]*R[1]*R[3],(R[4]*R[1])**2]])
    # Normalize each column to 1 (for probability)
    for i in range(4):
        TM[:,i] /= TM[:,i].sum()
    return TM

def randomer(n:int,lamdis:str):
    """Generate random set of rate parameters R and interaction parameter \lambda

    Parameters
    ----------
    n : int
        - No. of sets to generate
    lamdis : str
        - Sample distribution for lambda and delta
            - loguni
            - uni

    Returns
    -------
    R : NDArray (6,n)
        - Rate and Interaction Parameters
        - [0] index: p / stayOff
        - [1] index: q / stayOn
        - [2] index: r / off
        - [3] index: s / on
        - [4] index: l / coordination parameter (lambda)
        - [5] index: d / competition parameter (delta)

    Raises
    ------
    \lambda distribution undefined
    """
    # Generates uniform random sets
    R = np.random.rand(n,6)
    if lamdis=='loguni':
        # Converts \lambda and \delta to uniform in log scale over 0 to 100
        R[:,4:6] = np.power(10,4*(R[:,4:6]-0.5))
    elif lamdis=='uni':
        # Converts \lambda and \delta to uniform over 0 to 100
        R[:,4:6] = 0.01 + 99.99*R[:,4:6]
    else:
        raise Exception("\lambda distribution undefined")
    return R

def parm_sweeper(n:int,parm:str,lamdis:str='loguni',f_val:float=1.0,l_val=1):
    """Generate parameters set R such that it sweeps over parameter parm with others being constant

    Parameters
    ----------
    n : int
        - No. of sets to generate
    parm : str
        - parameter
            - p / stayOff
            - q / stayOn
            - r / off
            - s / on
            - l / coordination parameter (lambda)
            - d / competition parameter (delta)
    lamdis : str, optional
        - Sample distribution for lambda
            - loguni / Uniform over log scale (default)
            - uni / Uniform over linear scale
    f_val : int, optional
        - Default value for fixed parameters (other than lambda)
    l_val : int, optional
        - Default value for l if not varied

    Returns
    -------
    R : NDArray (6,n)
        - Rate and Interaction Parameters
        - [0] index: p / stayOff
        - [1] index: q / stayOn
        - [2] index: r / off
        - [3] index: s / on
        - [4] index: l / coordination parameter (lambda)
        - [5] index: d / competition parameter (delta)

    Raises
    ------
    \lambda distribution undefined
    """
    R = np.full((n,6),f_val)
    p = np.linspace(0,1,n)
    i = {'p':0 , 'q':1, 'r':2, 's':3, 'l':4, 'd':5}
    i = i[parm]
    if (parm == 'l' or parm == 'd'):
        if lamdis=='loguni':
            # Converts \lambda or \delta to uniform in log scale over 0 to 100
            R[:,i] = np.power(10,4*(p-0.5))
        elif lamdis=='uni':
            # Converts \lambda or \delta to uniform over 0 to 100
            R[:,i] = 0.01 + 99.99*p
        else:
            raise Exception("\lambda distribution undefined")
    else:
        R[:,i] = p
        R[:,4:6] = np.full((n,2),l_val) 
    return R

def parm_sweeper2D(n:int,parm,lamdis:str='loguni',f_val:float=1.0):
    """Generate parameters set R such that it sweeps over 2 parameters in parm with others being constant

    Parameters
    ----------
    n : int
        - No. of sets to generate per parameter 
    parm : tuple 
        - parameter
            - p / stayOff
            - q / stayOn
            - r / off
            - s / on
            - l / interaction parameter (lambda)
    lamdis : str, optional
        - Sample distribution for lambda
            - loguni / Uniform over log scale (default)
            - uni / Uniform over linear scale
    f_val : int, optional
        - Default value for fixed parameters

    Returns
    -------
    R : NDArray (6,n^2)
        - Rate and Interaction Parameters
        - [0] index: p / stayOff
        - [1] index: q / stayOn
        - [2] index: r / off
        - [3] index: s / on
        - [4] index: l / coordination parameter (lambda)
        - [5] index: d / competition parameter (delta)

    Raises
    ------
    \lambda distribution undefined
    """
    R = np.full((n**2,6),f_val)
    p = np.linspace(0,1,n+1)
    p=p[1:] # Drop 0
    if ('l' in parm) or ('d' in parm):
        if lamdis=='loguni':
            # Converts \lambda to uniform in log scale over 0 to 100
            q = np.power(10,4*(p-0.5))
        elif lamdis=='uni':
            # Converts \lambda to uniform over 0 to 100
            q = 0.01 + 99.99*p
        else:
            raise Exception("\lambda distribution undefined")
        if ('l' in parm) and ('d' in parm):
            p=q
    else:
        q=p
    i = {'p':0 , 'q':1, 'r':2, 's':3, 'l':4, 'd':5}
    parm=sorted(parm)
    parms=np.meshgrid(q,p)
    for j in range(2):
        R[:,i[parm[j]]] = parms[j].reshape(-1)  
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
    evec=evec.astype(float) #Temporary workaround 
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