import numpy as np
import pandas as pd
from MM import T, randomer, timeseries

## Input Parameters
t_max=1000 # Maximum time steps to simulate for 
fname='uni-sim.csv' # Filename to save as

fname='../raw_output/'+fname
RS=randomer(1)
tm=T(RS[0])

ts=timeseries(t_max,tm)
ts=pd.DataFrame(ts,columns=['Time','State'])
ts.to_csv(fname,index=False)
