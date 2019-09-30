#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#written by Anna Liznar in jupyter


# In[ ]:


#make plots for plus or minus strand
Xins_all={0: {'plus': 3, 'minus': 25}, 
 1: {'plus': 0, 'minus': 20}, 
 2: {'plus': 0, 'minus': 5}, 
 3: {'plus': 0, 'minus': 2},
 4: {'plus': 0, 'minus': 1}, 
 5: {'plus': 0, 'minus': 1}}

Xins_sm3={0: {'plus': 1, 'minus': 22}, 
 1: {'plus': 0, 'minus': 18}, 
 2: {'plus': 0, 'minus': 4}, 
 3: {'plus': 0, 'minus': 1},
 4: {'plus': 0, 'minus': 1}, 
 5: {'plus': 0, 'minus': 1}}

#make dataframes
df_Xins_sm3=pd.DataFrame.from_dict(Xins_sm3)
df_Xins_sm3=df_Xins_sm3.transpose()

df_Xins_all=pd.DataFrame.from_dict(Xins_all)
df_Xins_all=df_Xins_all.transpose()
#concat dfs to one with new keys
df_Xins_strand=pd.concat([df_Xins_all, df_Xins_sm3], axis=1, sort=False, keys=['all', 'sm3'])


# In[ ]:


#make plots for plus or minus strand
X1_all={0: {'plus': 3, 'minus': 17},
 1: {'plus': 1, 'minus': 13},
 2: {'plus': 1, 'minus': 6}, 
 3: {'plus': 0, 'minus': 2},
 4: {'plus': 0, 'minus': 4}, 
 
 5: {'plus': 0, 'minus': 2},
 6: {'plus': 0, 'minus': 1}, 
 7: {'plus': 9, 'minus': 14},
 8: {'plus': 3, 'minus': 5},
 9: {'plus': 5, 'minus': 1}, 
 
 10: {'plus': 9, 'minus': 14}, 
 11: {'plus': 9, 'minus': 14},
 12: {'plus': 1, 'minus': 4}, 
 13: {'plus': 5, 'minus': 1}, 
 
 14: {'plus': 9, 'minus': 14}, 
 15: {'plus': 4, 'minus': 12}, 
 16: {'plus': 9, 'minus': 17}}

X1_SM3={0: {'plus': 3, 'minus': 15},
 1: {'plus': 1, 'minus': 11},
 2: {'plus': 1, 'minus': 3}, 
 3: {'plus': 0, 'minus': 2},
 4: {'plus': 0, 'minus': 4}, 
 
 5: {'plus': 0, 'minus': 2},
 6: {'plus': 0, 'minus': 0}, 
 7: {'plus': 7, 'minus': 14},
 8: {'plus': 2, 'minus': 5},
 9: {'plus': 4, 'minus': 1}, 
 
 10: {'plus': 7, 'minus': 14}, 
 11: {'plus': 7, 'minus': 14},
 12: {'plus': 1, 'minus': 4}, 
 13: {'plus': 1, 'minus': 4}, 
 
 14: {'plus': 7, 'minus': 14}, 
 15: {'plus': 3, 'minus': 12}, 
 16: {'plus': 7, 'minus': 17}}


# In[ ]:


#make dataframe from dic
df_X1_sm3=pd.DataFrame.from_dict(X1_SM3)
df_X1_sm3=df_X1_sm3.transpose()


df_X1_all=pd.DataFrame.from_dict(X1_all)
df_X1_all=df_X1_all.transpose()

df_strand=pd.concat([df_X1_all, df_X1_sm3], axis=1, sort=False, keys=['all', 'sm3'])
#concat dataframes with new keys -> all and SM3


# In[ ]:


#make plot  X1
df_strand.plot.bar(rot=0)


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(15, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/Stretches/Strandness_X1.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


#make plot Xins
df_Xins_strand.plot.bar(rot=0)


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(7, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/Stretches/Strandness_Xins.pdf', bbox_inches='tight')
plt.show()
plt.close()

