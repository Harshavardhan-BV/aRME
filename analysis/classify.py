#%%
import numpy as np
import pandas as pd
#%%
parm='uni'
eps=0.1
cats=['p2ex','p0ex','comp','coor','semi','indep']
#%%
rdatname='../raw_output/'+parm+'.csv'
#%%
df=pd.read_csv(rdatname)
def ptsin(df,cat,eps):
    if cat=='p0ex':
        return df.p00>=1-2*eps
    elif cat=='p2ex':
        return df.p11>=1-2*eps
    elif cat=='comp':
        return df.p11<np.square(1-np.sqrt(df.p00)-eps)
    elif cat=='coor':
        return df.p11>=1-df.p00-eps
    elif cat=='semi':
        return df.p11>np.square(1-np.sqrt(df.p00)+eps)
    elif cat=='indep':
        return df.p11<=np.square(1-np.sqrt(df.p00)+eps)
    else:
        return 0
    
for cat in cats:
    adatname='../analysed_data/'+parm+'/'+cat+'.csv'
    yes=ptsin(df,cat,eps)
    ndf,df=df[yes],df[~yes]
    ndf.to_csv(adatname,index=False)

# %%
