#%%
import common_fn as cf
import matplotlib.pyplot as plt
from itertools import combinations
# plt.style.use('dark_background')
#%%
cats=['p2ex','p0ex','comp','coor','semi','indep']
parms=('p','q','r','s','l')
for parm in combinations(parms,2):
    parm_name='sweep-'+''.join(parm)
    #%%
    cf.mkdirs(parm_name)
    # %%
    cf.classify(parm_name,cats,normalize='max')
    # %%
    cf.p2vp0(parm_name)
    # %%
    cf.p2vp0_cat(parm_name,cats)
    # %%
    cf.parm_box(parm_name,cats)
    # %%
    cf.p2vp0_2Dparm(parm)
