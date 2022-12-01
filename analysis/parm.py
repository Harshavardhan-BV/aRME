#%%
import common_fn as cf
import matplotlib.pyplot as plt
from itertools import combinations
#%%
cats=['p2ex','p0ex','comp','coor','semi','indep']
parms=('p','q','r','s','l')
for parm in parms:
    parm_name='sweep-'+parm
    #%%
    cf.mkdirs(parm_name)
    # %%
    cf.classify(parm_name,cats,normalize='max',force=True)
    # %%
    cf.p2vp0(parm_name,force=True)
    # %%
    cf.p2vp0_cat(parm_name,cats,force=True)
    # %%
    cf.parm_box(parm_name,cats,force=True)
    # %%
    cf.p2vp0_parm(parm_name,parm,force=True)