#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#written by Anna Liznar in jupyter


# In[ ]:


def make_dff():
    """
    make df from mapping events full length piRNA and seed
    """
    seed={'X1':{'all':22192, 
                'full': 12693,
                'full_one_mm':12721,
                'full_two_mm':12753,
                'full_three_mm':12775,
                'seed': 12693},
        'Xins':{'all':22192, 
                'full': 12693,
                'full_one_mm':12721,
                'full_two_mm':12753,
                'full_three_mm':12775,
                'seed': 12693}
         }
    df_seed=pd.DataFrame.from_dict(seed)
    #make df from dic
    df_s=pd.DataFrame()
    df_s['piRNAs']=df_seed['X1']+df_seed['Xins']/2 #average from X1 and Xins
    
    df_s=df_s.transpose() #transpose df
    #calc percentages
    df_s['full 0 mm[%]']=df_s['full']/df_s['all']*100
    df_s['seed [%]']=df_s['seed']/df_s['all']*100
    df_s['full 1 mm [%]']=df_s['full_one_mm']/df_s['all']*100
    df_s['full 2 mm [%]']=df_s['full_two_mm']/df_s['all']*100
    df_s['full 3 mm [%]']=df_s['full_three_mm']/df_s['all']*100

    pl=pd.DataFrame()
    pl=df_s[['full 0 mm[%]','full 1 mm [%]', 'full 2 mm [%]', 'full 3 mm [%]']]
    pl=pl.transpose()

    return pl
    
make_dff()


# In[ ]:


def plot_df():
    """
    plot from df from make df 
    """
    df=make_dff()

    df.plot(kind='line', rot=0)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.ylabel('Percentage of targets', fontsize=12)
    plt.legend(fontsize=12)
    axes = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(3, 2.5)
    axes.set_xticklabels(['0 MM ','0 MM ','1 MM', '2 MM', '3 MM'])
    plt.hlines(y=57.19628695,xmin=0, xmax=df['full 3 mm [%]':'piRNAs']+3, colors='r')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('/Seed/Line_new_0123.pdf', bbox_inches='tight')
    plt.show()
    plt.close()
    
plot_df()

