#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[ ]:


#dictionary with how many piRNAs were found in stretches
X1={0:{'allpiRNAs': 20, 'SM3': 18, 'length': 298}, 
 1: {'allpiRNAs': 14, 'SM3': 12, 'length': 367}, 
 2: {'allpiRNAs': 7, 'SM3': 4, 'length': 105}, 
 3: {'allpiRNAs': 2, 'SM3': 2, 'length': 41}, 
 4: {'allpiRNAs': 4, 'SM3': 4, 'length': 41},
    
 5: {'allpiRNAs': 2, 'SM3': 2, 'length': 33}, 
 6: {'allpiRNAs': 1, 'SM3': 0, 'length': 32}, 
 7: {'allpiRNAs': 23, 'SM3': 21, 'length': 251}, 
 8: {'allpiRNAs': 8, 'SM3': 7, 'length': 87}, 
 9:  {'allpiRNAs': 6, 'SM3': 5, 'length': 47}, 
    
 10: {'allpiRNAs': 23, 'SM3': 21, 'length': 251}, 
 11: {'allpiRNAs': 23, 'SM3': 21, 'length': 251}, 
 12: {'allpiRNAs': 5, 'SM3': 5, 'length': 71}, 
 13: {'allpiRNAs': 6, 'SM3': 5, 'length': 47}, 
    
 14: {'allpiRNAs': 23, 'SM3': 21, 'length': 251}, 
 15: {'allpiRNAs': 16, 'SM3': 15, 'length': 231}, 
 16: {'allpiRNAs': 26, 'SM3': 24, 'length': 444}} 

#dictionary with how many piRNAs were found in stretches
Xins={0: {'allpiRNAs': 28, 'SM3': 23, 'length': 299},
 1: {'allpiRNAs': 20, 'SM3': 18, 'length': 273}, 
 2:  {'allpiRNAs': 5, 'SM3': 4, 'length': 40},
 3:  {'allpiRNAs': 2, 'SM3': 1, 'length': 30}, 
 4:  {'allpiRNAs': 1, 'SM3': 1, 'length': 31},
 5:  {'allpiRNAs': 1, 'SM3': 1, 'length': 32}}


# In[ ]:


#X1
df_X1=pd.DataFrame.from_dict(X1) #make df from dic
df_X1=df_X1.transpose() #transpose
print(df_X1['SM3'].sum(), df_X1['allpiRNAs'].sum()) #print in jupyter 
#total amounts of SM3 piRNAs and piRNAs from unc RNAi kd

#Xins
df_Xins=pd.DataFrame.from_dict(Xins) #make df from dic
df_Xins=df_Xins.transpose()


# In[ ]:


#X1
from itertools import cycle, islice
my_colors = list(islice(cycle([ 'y', 'k']), None, len(df_X1[['allpiRNAs','SM3']])))   

df_X1[['allpiRNAs','SM3']].plot(kind='bar',  color=my_colors)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(9, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/Stretches/Amount_piRNAs_X1.pdf', bbox_inches='tight')
plt.show()
plt.close()

df_X1['length'].plot.bar(rot=0, color='g')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(9, 3)

plt.savefig('/Stretches/Length_X1.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


#Xins
from itertools import cycle, islice
my_colors = list(islice(cycle([ 'y', 'k']), None, len(df_Xins[['allpiRNAs','SM3']])))   

df_Xins[['allpiRNAs','SM3']].plot(kind='bar',  color=my_colors)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(5, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/Stretches/Amount_piRNAs_Xins.pdf', bbox_inches='tight')
plt.show()
plt.close()


df_Xins['length'].plot.bar(rot=0, color='g')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(5, 3)

plt.savefig('/Stretches/Length_Xins.pdf', bbox_inches='tight')
plt.show()
plt.close()

