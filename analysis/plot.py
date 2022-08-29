#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['text.usetex'] = True
#%%
parm='loguni'
#%%
datname='../raw_output/'+parm+'.csv'
figname='../figures/'+parm+'-p2vp0.svg'
df=pd.read_csv(datname)
#%%
x=np.linspace(0,1,100)
y_coord=1-x
y_indep=np.square(1-np.sqrt(x))
#%%
plt.figure(1,[10,10])
plt.scatter(df.p00,df.p11,color='tab:red',s=1,label='Sim')
plt.plot(x,y_coord,label='Coord')
plt.plot(x,y_indep,label='Indep')

plt.xlabel(r'$p_0$')
plt.ylabel(r'$p_2$')
plt.legend()
plt.savefig(figname)

# %%
cats=['p2ex','p0ex','comp','coor','semi','indep']
figname='../figures/'+parm+'-classified.svg'
plt.clf()
plt.figure(1,[10,10])
for cat in cats:
    datfname='../analysed_data/'+parm+'/'+cat+'.csv'
    df=pd.read_csv(datfname)
    plt.scatter(df.p00,df.p11,s=1,label=cat)
plt.xlabel(r'$p_0$')
plt.ylabel(r'$p_2$')
plt.legend()
plt.savefig(figname)
# %%
cats=['p2ex','p0ex','comp','coor','semi','indep']
figname='../figures/'+parm+'-parms.svg'
ndf=pd.DataFrame()
for cat in cats:
    datfname='../analysed_data/'+parm+'/'+cat+'.csv'
    df=pd.read_csv(datfname)
    tdf=pd.melt(df.loc[:,'p':'l'])
    tdf['Category']=cat
    ndf=pd.concat([ndf,tdf])
#%%
tdf=ndf[ndf.variable!='l']
plt.clf()
plt.figure(1,[10,10])
ax = sns.catplot(x="variable", y="value",hue="Category", data=tdf, kind='box')
plt.savefig(figname)

figname='../figures/'+parm+'-lambda.svg'

tdf=ndf[ndf.variable=='l']
tdf.loc[:,'value']=np.log10(tdf['value'])
plt.clf()
plt.figure(1,[10,10])
ax = sns.catplot(x="variable", y="value",hue="Category", data=tdf, kind='box')
plt.xlabel(r'$\lambda$')
plt.ylabel('log(value)')
plt.savefig(figname)
