# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams["svg.hashsalt"]=''
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"]=30
plt.style.use('dark_background')
# %%
def sweep2D(parms,i=1):
    parm_name='sweep-'+''.join(parms)
    rdatname='../raw_output/'+parm_name+'.csv'
    df=pd.read_csv(rdatname)
    df['l']=np.log10(df['l'])
    figname='../writing/Slides/figures/'+parm_name+'-'+parms[i]+'.pdf'
    fig=plt.figure(figsize=(10,10))
    x=np.linspace(0,1)
    y=np.square(1-np.sqrt(x))
    plt.plot(x,y,label='indp',color='w')
    sns.scatterplot(data=df,x='p00',y='p11',hue=parms[i],edgecolor='none',palette='coolwarm')
    plt.xlabel(r'$p_0$')
    plt.ylabel(r'$p_2$')
    plt.legend()
    plt.savefig(figname)
    fig.clf()
    plt.close(fig)
# %%
def sweep(parms):
    parm_name='sweep-'+parms
    rdatname='../raw_output/'+parm_name+'.csv'
    df=pd.read_csv(rdatname)
    df['l']=np.log10(df['l'])
    figname='../writing/Slides/figures/'+parm_name+'.pdf'
    fig=plt.figure(figsize=(10,10))
    x=np.linspace(0,1)
    y=np.square(1-np.sqrt(x))
    plt.plot(x,y,label='indp',color='w')
    sns.scatterplot(data=df,x='p00',y='p11',hue=parms,edgecolor='none',palette='coolwarm')
    plt.xlabel(r'$p_0$')
    plt.ylabel(r'$p_2$')
    plt.legend()
    plt.savefig(figname)
    fig.clf()
    plt.close(fig)
# %%
sweep2D(('p','l'))
sweep2D(('q','l'))
sweep2D(('r','l'))
sweep2D(('s','l'))
# %%
sweep2D(('p','l'),0)
sweep2D(('q','l'),0)
sweep2D(('r','l'),0)
sweep2D(('s','l'),0)
# %%
sweep('p')
sweep('q')
sweep('r')
sweep('s')
sweep('l')
# %%
figname='../writing/Slides/figures/loguni-p2vp0-classified.pdf'
fig=plt.figure(figsize=(10,10))
cats=['p2ex','p0ex','comp','coor','semi','indep']
for cat in cats:
    datfname='../analysed_data/loguni/'+cat+'.csv'
    df=pd.read_csv(datfname)
    plt.scatter(df.p00,df.p11,label=cat)
plt.xlabel(r'$p_0$')
plt.ylabel(r'$p_2$')
plt.legend()
fig.tight_layout()
fig.savefig(figname)
fig.clf()
plt.close(fig)

# %%
figname='../writing/Slides/figures/loguni-parms.pdf'
fig, ax = plt.subplots(1,2,figsize=(20,10),gridspec_kw={'width_ratios': [4, 1]})
df=pd.DataFrame()
for cat in cats:
    datfname='../analysed_data/loguni/'+cat+'.csv'
    ndf=pd.read_csv(datfname)
    tdf=pd.melt(ndf.loc[:,'p':'l'])
    tdf['Category']=cat
    df=pd.concat([df,tdf])
yes=df.variable!='l'
df,dfl=df[yes],df[~yes]
b1=sns.boxplot(x="variable", y="value",hue="Category", data=df,ax=ax[0])
ax[1].set_yscale('log')
b2=sns.boxplot(x="variable", y="value",hue="Category", data=dfl,ax=ax[1])
ax[0].legend([],[], frameon=False)
b2.set(xlabel=None,ylabel=None,xticklabels=[r'$\lambda$'])
ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
fig.tight_layout()
fig.savefig(figname)
fig.clf()
plt.close(fig)