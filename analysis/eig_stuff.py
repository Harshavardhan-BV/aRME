#%%
import common_fn as cf
import matplotlib.pyplot as plt
# plt.style.use('dark_background')
#%%
parm_name='sweep-l'
cats=['p2ex','p0ex','comp','coor','semi','indep']
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
