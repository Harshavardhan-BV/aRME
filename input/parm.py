import numpy as np
import pandas as pd
from itertools import combinations
from MM import *

## Input Parameters
nP=100 # No. of parameter sets 
parms=('p','q','r','s','l')
for parm in parms:
    fname='sweep-'+parm+'.csv' # Filename to save as

    fname='../raw_output/'+fname
    RS=parm_sweeper(nP,parm,f_val=1.0)
    y=np.empty([0,4]) #Empty array to append eigenvector to
    RS=RS[1:]

    for i in range(len(RS)):
        tm=T(RS[i]) # Transition matrix for a set of rate
        ev=markeig(tm) # Get eigenvector
        y=np.append(y,[ev],axis=0) # Eigenvectors added row wise

    y=np.append(RS,y,axis=1) # Rates and respective eigenvectors column wise
    df=pd.DataFrame(y,columns=['p','q','r','s','l','p00','p01','p10','p11']) # Convert to dataframe
    df.to_csv(fname,index=False)



