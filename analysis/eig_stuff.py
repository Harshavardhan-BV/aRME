#%%
import common_fn as cf
#%%
parm_name='uni'
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
cf.parm_box(parm_name,cats,force=True)
