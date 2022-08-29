#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['text.usetex'] = True
#%%
parm='uni'
# %%
datname='../raw_output/'+parm+'-sim.csv'
figname='../figures/timeseries-'+parm+'.svg'
df=pd.read_csv(datname)
df['G1']=2*(df.State%2)
df['G2']=df.State//2

# %%
plt.figure(1,[10,10])
plt.plot(df.Time[0:100],df.G1[0:100],label='G1')
plt.plot(df.Time[0:100],df.G2[0:100],label='G2')

plt.xlabel('Time')
plt.ylabel('On/Off')
plt.legend()
plt.savefig(figname)