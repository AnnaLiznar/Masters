#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
#written by Anna Liznar in jupyter


# In[102]:


path_1='/1_1u/'
df=pd.read_csv(path_1+'thresh_pairs.aligned.tab', header=None, sep='\t')

df=df.loc[:,[3,8,9]]
df=df.drop_duplicates()

path_1='/3_2u/'
ddf=pd.read_csv(path_1+'thresh_pairs.aligned.tab', header=None, sep='\t')

ddf=ddf.loc[:,[3,8,9]]
ddf=ddf.drop_duplicates()

path_1='/5_3u/'
dff=pd.read_csv(path_1+'thresh_pairs.aligned.tab', header=None, sep='\t')

dff=dff.loc[:,[3,8,9]] #9 sequence 8reverse complement matches mismacthes
ddf=ddf.drop_duplicates()

al=df.append(ddf, ignore_index=True)
al=al.append(dff, ignore_index=True) #column 8 reverse complement of binding


al['Number']=al.index #index as column
binding_df=al.loc[:, ['Number',8]]


#make list of reverse piRNA length
idx=list(range(0,33))
reverse=[]
for i in reversed(idx):
    reverse.append(i)
print(reverse) #reversed [[:-1]]

bd_df=binding_df[8].apply(lambda x: pd.Series(list(x))) #make columns for each position
tt=bd_df.transpose() #transpose
tt=tt.reindex(reverse) #use reverse index
tt=tt.transpose() #re-transpose

df=tt.replace({'1':3, 'm':2, 'g':0, 's':-1, '0':-1, '.':0}) #exchange to integers
df=df.fillna(value=0) #fill na with 0

df=df.apply(pd.to_numeric) #no strings!

##make heatmap
fig, ax = plt.subplots()
im = ax.pcolor(df, cmap='Greys', linewidths=1)
fig.colorbar(im)
ax.patch.set(hatch='', edgecolor='black')
fig.set_size_inches(8, 10)

#plt.savefig('/path/greys.heatmap_pattern.png',dpi=400)
plt.show()


# In[ ]:


#same for shuffled with file called:
df=pd.read_csv(path_1+'pairs.shuffled.aligned.tab', header=None, sep='\t')
df['Number']=df.index #index as column
binding_df=df.loc[:, ['Number',8]]


#make list of reverse piRNA length
idx=list(range(0,33))
reverse=[]
for i in reversed(idx):
    reverse.append(i)
print(reverse) #reversed [[:-1]]

bd_df=binding_df[8].apply(lambda x: pd.Series(list(x))) #make columns for each position
tt=bd_df.transpose() #transpose
tt=tt.reindex(reverse) #use reverse index
tt=tt.transpose() #re-transpose

df=tt.replace({'1':3, 'm':2, 'g':0, 's':-1, '0':-1, '.':0}) #exchange to integers
df=df.fillna(value=0) #fill na with 0

df=df.apply(pd.to_numeric) #no strings!

##make heatmap
fig, ax = plt.subplots()
im = ax.pcolor(df, cmap='Greys', linewidths=1)
fig.colorbar(im)
ax.patch.set(hatch='', edgecolor='black')
fig.set_size_inches(8, 10)

#plt.savefig('/path/greys.heatmap_pattern_shuffled.png',dpi=400)
plt.show()

