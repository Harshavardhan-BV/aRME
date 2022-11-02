import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
plt.rcParams["svg.hashsalt"]=''
plt.rcParams['text.usetex'] = True
plt.rcParams["font.size"]=22

def mkdirs(parm_name,):
    try:
        os.makedirs("../figures/"+parm_name)
    except:
        pass
    try:
        os.makedirs("../analysed_data/"+parm_name)
    except:
        pass

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

def normdf(df,normalize):
    if normalize=='max':
        return df.loc[:,'p':'s'].max(axis=1)
    if normalize=='p':
        return df.loc[:,'p']
    if normalize=='q':
        return df.loc[:,'q']
    if normalize=='r':
        return df.loc[:,'r']
    if normalize=='s':
        return df.loc[:,'s']

def classify(parm_name,cats,eps=0.1,force=False,normalize=''):
    # Read the raw data
    rdatname='../raw_output/'+parm_name+'.csv'
    df=pd.read_csv(rdatname)
    # For every category
    for cat in cats:
        adatname='../analysed_data/'+parm_name+'/'+cat+'.csv'
        # Find points that satisfy category conditions
        yes=ptsin(df,cat,eps)
        # New dataframe with points that satisfy & remove from big dataframe 
        ndf,df=df[yes],df[~yes]
        if normalize:
            ndf.loc[:,'p':'s']=ndf.loc[:,'p':'s'].div(normdf(ndf,normalize),axis=0)
        if not os.path.isfile(adatname) or force: 
            ndf.to_csv(adatname,index=False)

def p2vp0(parm_name,plot_lines=True,force=False):
    # Read the raw data
    rdatname='../raw_output/'+parm_name+'.csv'
    df=pd.read_csv(rdatname)
    # Create a figure and plot points
    figname='../figures/'+parm_name+'/p2vp0.svg'
    fig=plt.figure(figsize=(10,10))
    plt.scatter(df.p00,df.p11,color='tab:red',s=1,label='Sim')
    # Plot reference lines for coordinated and independent
    if plot_lines:
        x=np.linspace(0,1,100)
        y_coord=1-x
        y_indep=np.square(1-np.sqrt(x))
        plt.plot(x,y_coord,label='Coord')
        plt.plot(x,y_indep,label='Indep')
    # Labels and legends
    plt.xlabel(r'$p_0$')
    plt.ylabel(r'$p_2$')
    plt.legend()
    fig.tight_layout()
    # Save figure
    if not os.path.isfile(figname) or force: 
        fig.savefig(figname)
    # Clear figure
    fig.clf()
    plt.close(fig)

def p2vp0_cat(parm_name,cats,force=False):
    # Create a figure
    figname='../figures/'+parm_name+'/p2vp0-classified.svg'
    fig=plt.figure(figsize=(10,10))
    # For every category
    for cat in cats:
        # Read categorized data in analysed_data
        datfname='../analysed_data/'+parm_name+'/'+cat+'.csv'
        df=pd.read_csv(datfname)
        # Plot points
        plt.scatter(df.p00,df.p11,s=1,label=cat)
    # Labels and legends
    plt.xlabel(r'$p_0$')
    plt.ylabel(r'$p_2$')
    plt.legend()
    fig.tight_layout()
    # Save figure
    if not os.path.isfile(figname) or force: 
        fig.savefig(figname)
    # Clear figure
    fig.clf()
    plt.close(fig)

def longdf(parm_name,cats):
    ndf=pd.DataFrame()
    for cat in cats:
        datfname='../analysed_data/'+parm_name+'/'+cat+'.csv'
        df=pd.read_csv(datfname)
        tdf=pd.melt(df.loc[:,'p':'d'])
        tdf['Category']=cat
        ndf=pd.concat([ndf,tdf])
    return ndf

def parm_box(parm_name,cats,force=False):
    # Create a figure
    figname='../figures/'+parm_name+'/parms.svg'
    fig, ax = plt.subplots(1,2,figsize=(25,10),gridspec_kw={'width_ratios': [4, 2]})
    df=longdf(parm_name,cats)
    yes=np.logical_and(df.variable!='l',df.variable!='d')
    df,dfl=df[yes],df[~yes]
    # Plot boxplots, lambda scale different so in a subplot
    b1=sns.boxplot(x="variable", y="value",hue="Category", data=df,ax=ax[0])
    ax[1].set_yscale('log')
    b2=sns.boxplot(x="variable", y="value",hue="Category", data=dfl,ax=ax[1])
    # Labels and legends
    ax[0].legend([],[], frameon=False)
    b2.set(xlabel=None,ylabel=None,xticklabels=[r'$\lambda$',r'$\delta$'])
    ax[1].legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    fig.tight_layout()
    # Save figure
    if not os.path.isfile(figname) or force: 
        fig.savefig(figname)
    # Clear figure
    fig.clf()
    plt.close(fig)

def timeseries(parm_name,force=False):
    rdatname='../raw_output/'+parm_name+'.csv'
    figname='../figures/'+parm_name+'/timeseries.svg'
    df=pd.read_csv(rdatname)
    df['G1']=2*(df.State%2)
    df['G2']=df.State//2
    fig,ax=plt.figure(1,[20,10])
    plt.plot(df.Time,df.G1,label='G1')
    plt.plot(df.Time,df.G2,label='G2')
    # Save figure
    if not os.path.isfile(figname) or force: 
        fig.savefig(figname)
    # Clear figure
    fig.clf()
    plt.close(fig)

def p2vp0_parm(parm_name,plt_parm,logl=True,force=False):
    # Read the raw data
    rdatname='../raw_output/'+parm_name+'.csv'
    df=pd.read_csv(rdatname)
    if logl:
        df['l']=np.log(df['l'])
    # Create a figure and plot points
    figname='../figures/'+parm_name+'/p2vp0-'+plt_parm+'.svg'
    fig=plt.figure(figsize=(10,10))
    sns.scatterplot(data=df,x='p00',y='p11',hue=plt_parm,edgecolor='none',palette='coolwarm')
    # Labels and legends
    plt.xlabel(r'$p_0$')
    plt.ylabel(r'$p_2$')
    plt.legend()
    fig.tight_layout()
    # Save figure
    if not os.path.isfile(figname) or force: 
        fig.savefig(figname)
    # Clear figure
    fig.clf()
    plt.close(fig)
