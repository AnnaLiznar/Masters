#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#written by Anna Liznar in jupyter


# In[ ]:


#make dic for Duplicates and Singles (PCR)
#obtained by using umitools
#then wc -l
Dup={'ctrl_1':{'Duplicates':52738812, 'Single': 6313196},
    'ctrl_2':{'Duplicates':64744568, 'Single': 1247232},
    'piRNA_1':{'Duplicates':64154364, 'Single': 1610884},
    'piRNA_2':{'Duplicates':62097884, 'Single': 1571516},
    'piRNA_3':{'Duplicates':51345432, 'Single': 1189556}}
ddf=pd.DataFrame.from_dict(Dup) #make df from dic
ddf=ddf.transpose() #transpose
ddf['total']=ddf['Duplicates']+ddf['Single'] #total 
ddf['Duplicates [%]']=ddf['Duplicates']/ddf['total']*100 #percentages
ddf['Single [%]']=ddf['Single']/ddf['total']*100


# In[ ]:


#plot kind bar
ddf[['Duplicates [%]','Single [%]']].plot.bar(stacked=True,rot=0)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/path/PCR_duplicates.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


#potential targets 
#all possible targets 
#targets with piRNA partners
#actual found targets
#dictionary:
ff={'ctrl_1':{'Potential':57529010, 'piRNA_partner':126852, 'Actual': 46},
    'ctrl_2':{'Potential':11124940, 'piRNA_partner':16263, 'Actual': 69},
    'piRNA_1':{'Potential':23750046, 'piRNA_partner':49842, 'Actual': 26143},
    'piRNA_2':{'Potential':21290694,  'piRNA_partner':61861, 'Actual': 18989},
    'piRNA_3':{'Potential':21149866, 'piRNA_partner':52530,  'Actual': 36086},}
podf=pd.DataFrame.from_dict(ff) #df from dic
podf=podf.transpose() #transpose
podf['total']=podf['Actual']+podf['Potential']+podf['piRNA_partner'] #total
podf['Actual [%]']=podf['Actual']/podf['total']*100 #percentage
podf['Potential [%]']=podf['Potential']/podf['total']*100
podf['piRNA partner [%]']=podf['piRNA_partner']/podf['total']*100


# In[ ]:


#plot
podf[['Potential [%]','Actual [%]','piRNA partner [%]']].plot.bar(stacked=True,rot=0)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
axes = plt.gca()
axes.set_ylim([95,100.1])
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
#plt.savefig('/path/Potential_targets.pdf', bbox_inches='tight')
plt.show()
plt.close()

