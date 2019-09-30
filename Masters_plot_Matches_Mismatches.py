#!/usr/bin/env python
# coding: utf-8

# In[1]:


import operator
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqIO.FastaIO import SimpleFastaParser
import matplotlib.pyplot as plt
import seaborn as sns


# In[ ]:


path_1='/path/piRNA_1/'

df=pd.read_csv(path_1+'pairs.aligned.tab', header=None, sep='\t')

df=df.loc[:,[3,8,9]] #take only sequences and events (matches/mismatches/gaps)
df=df.drop_duplicates()

path_1='/path/piRNA_2/'

ddf=pd.read_csv(path_1+'pairs.aligned.tab', header=None, sep='\t')

ddf=ddf.loc[:,[3,8,9]]
ddf=ddf.drop_duplicates()
#print(ddf)


path_1='/path/piRNA_3/'

dff=pd.read_csv(path_1+'pairs.aligned.tab', header=None, sep='\t')

dff=dff.loc[:,[3,8,9]]
ddf=ddf.drop_duplicates()


al=df.append(ddf, ignore_index=True)
al=al.append(dff, ignore_index=True) #column 8 reverse complement of binding

al['Number']=al.index
seq_df=al.loc[:,['Number',9]]

binding_df=al.loc[:, ['Number',8]]

#manipulate for plot
bd_df=binding_df[8].apply(lambda x: pd.Series(list(x))) #for each position make column

columns=list(bd_df.columns) #make list of all columns

final=pd.DataFrame() #new dataframe 
for c in columns:
    dd=pd.DataFrame(bd_df[c].value_counts()) #count each event for all positions
    dd=dd.transpose() #transpose df
    final=final.append(dd, sort=True) #append to final 
    

fif=final.iloc[:,1:4] #retain only matches, mismatches and gaps in df (fif)
fit=fif.fillna(0) #fill na with zero

fif=fif.rename(columns={'0':'Mismatch', '1': 'Match', 'g':'Gap'}) #rename
p=list(range(33)) #make list from 0 to 32
p=p[::-1] #reverse the list from 32 to 0
fif=fif.reindex(p) #make new index since this is from 3' and not 5'
fif=fif.reset_index(drop=True) #new: 0 is 5' end


# In[2]:


###make plot kind bar 
pal='rgbkymc'
p=fif.plot.bar(alpha=0.4, rot=90, title='Alignments',color=pal) #stacked=True, 
p.set_ylabel("Amount of Matches/Mismatches and Gaps", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig = plt.gcf()
fig.set_size_inches(15, 5)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
#plt.savefig('/path/MM_M_Gap.pdf', bbox_inches='tight')

plt.show()
plt.close()
    


# In[276]:


#plot with percentages
#matches and mismatches

p=list(range(33))
p=p[::-1]

fin=pd.DataFrame()
final=final.fillna(0)
fin['Matches [%]']=(final['1']/final['1'].sum())*100 #calc percentage
fin['Mismatches [%]']=final['0']/final['0'].sum()*100
fin['Matches outside binding area [%]']=final['m']/final['m'].sum()*100
fin['Mismatches outside binding area [%]']=final['s']/final['s'].sum()*100
fin=fin.reindex(p)
fin=fin.reset_index(drop=True)
#reindex with list reverse list 
#reset index 
#plot kind line

p=fin.plot( alpha=0.4, rot=90, title='Alignments', linewidth=3.0)
p.set_ylabel("Matches/Mismatches and Gaps [%]", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

leg = plt.legend()
# get the lines and texts inside legend box
leg_lines = leg.get_lines()
leg_texts = leg.get_texts()
# bulk-set the properties of all lines and texts
plt.setp(leg_lines, linewidth=4)
plt.setp(leg_texts, fontsize=12)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
fig = plt.gcf()
fig.set_size_inches(8.5, 5)
#plt.savefig('/path/MM_M.pdf', bbox_inches='tight')

plt.show()
plt.close()

