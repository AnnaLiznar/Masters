#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar in jupyter
# In[ ]:


import subprocess, os, sys
import re
from getopt import getopt
import operator
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqIO.FastaIO import SimpleFastaParser
from glob import glob
import matplotlib.pyplot as plt
from Bio.Alphabet import IUPAC
import string 
import seaborn as sns


# In[ ]:


path_1='/data1/Master/CLASH/results/'


df=pd.read_csv(path_1+'thresh_pairs.aligned.tab', header=None, sep='\t')
#print (df)

df=df.loc[:,[3,8,9]]
df=df.drop_duplicates()
#print(df)

path_1='/data1/Master/CLASH/Analysis/results/'

ddf=pd.read_csv(path_1+'thresh_pairs.aligned.tab', header=None, sep='\t')
#print (df)

ddf=ddf.loc[:,[3,8,9]]
ddf=ddf.drop_duplicates()
#print(ddf)


path_1='/data1/Master/CLASH/Analysis/results/'

dff=pd.read_csv(path_1+'thresh_pairs.aligned.tab', header=None, sep='\t')
#print (df)

dff=dff.loc[:,[3,8,9]]
ddf=ddf.drop_duplicates()
#print(dff)

al=df.append(ddf, ignore_index=True)
al=al.append(dff, ignore_index=True) #column 8 reverse complement of binding


al['Number']=al.index
seq_df=al.loc[:,['Number',9]]


binding_df=al.loc[:, ['Number',8]]


# In[ ]:


idx=list(range(0,33))
reverse=[]
for i in reversed(idx):
    reverse.append(i)


# In[ ]:


bd_df=binding_df[8].apply(lambda x: pd.Series(list(x)))
print(bd_df.head(5))
tt=bd_df.transpose()
tt=tt.reindex(reverse)
tt=tt.transpose()


df=tt.replace({'1':3, 'm':2, 'g':0, 's':-1, '0':-1, '.':0})
df=df.fillna(value=0)


# In[ ]:


df=df.apply(pd.to_numeric)


# In[ ]:


fig, ax = plt.subplots()
im = ax.pcolor(df, cmap='Greys', linewidths=1)
fig.colorbar(im)

ax.patch.set(hatch='', edgecolor='black')
fig.set_size_inches(8, 10)


plt.savefig('/data1/greys.heatmap_pattern.png',dpi=400)
plt.show()


# In[ ]:


columns=list(bd_df.columns)
#columns


# In[ ]:


final=pd.DataFrame()

for c in columns:
    dd=pd.DataFrame(bd_df[c].value_counts())
    dd=dd.transpose()
    final=final.append(dd, sort=True)
    


# In[ ]:


fif=final.iloc[:,1:4]
fif#['1']
fit=fif.fillna(0)


# In[ ]:


fif=fif.rename(columns={'0':'Mismatch', '1': 'Match', 'g':'Gap'})
p=list(range(33))
p=p[::-1]
fif=fif.reindex(p)
fif=fif.reset_index(drop=True)


# In[ ]:



pal='rgbkymc'
p=fif.plot.bar(stacked=True, alpha=0.4, rot=90, title='Alignments',color=pal)
p.set_ylabel("Amount of Matches/Mismatches and Gaps", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig = plt.gcf()
fig.set_size_inches(8.5, 5)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
#plt.savefig('/data1/MM_M_Gap.pdf', bbox_inches='tight')

plt.show()
plt.close()
    


# In[ ]:


finl=final.iloc[:,1:]
finl=finl.rename(columns={'0':'Mismatch', '1': 'Match', 'g':'Gap',
                              's': 'Match outside binding area',
                              'm': 'Mismatch outside binding area'})

pal='rgbkymc'
p=finl.plot.bar(stacked=True, alpha=0.4, rot=90, title='Alignments reverse complement',color=pal)
p.set_ylabel("Amount of Matches/Mismatches and Gaps", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig = plt.gcf()
fig.set_size_inches(8.5, 5)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
#plt.savefig('/data1/Figures_Anna_pattern/MM_M_Gap_S_M.pdf', bbox_inches='tight')

plt.show()
plt.close()


# In[ ]:



plt.plot(fit['1'])
plt.plot(fit['g'])
plt.plot(fit['0'])
plt.ylabel("Amount of Matches/Mismatches and Gaps", fontsize=12)
plt.xlabel("Position in piRNA", fontsize=12)
axes = plt.gca()
#axes.set_ylim([-5,55])
#axes.set_xlim([-5, 100])
plt.xticks(fontsize=12, rotation=90)
plt.yticks(fontsize=12)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
#plt.savefig('/data1/eCLASH//plot_piRNA2_range12_with_names.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


summe=final.sum(axis=0)
summe=pd.DataFrame(summe)
summe['Total']=summe[0]
summe=summe.transpose()
summe=summe.drop([0])
summe


# In[ ]:


p=list(range(33))
p=p[::-1]

fin=pd.DataFrame()
final=final.fillna(0)
fin['Matches [%]']=(final['1']/final['1'].sum())*100
fin['Mismatches [%]']=final['0']/final['0'].sum()*100
#fin['Gaps [%]']=final['g']/final['g'].sum()*100
fin['Matches outside binding area [%]']=final['m']/final['m'].sum()*100
fin['Mismatches outside binding area [%]']=final['s']/final['s'].sum()*100
fin=fin.reindex(p)
fin=fin.reset_index(drop=True)

print(fin.head())


# In[ ]:


p=list(range(33))
p=p[::-1]
fi=pd.DataFrame()
final=final.fillna(0)
fi['Matches [%]']=((final['1']+final['m'])/(final['1']+final['m']).sum())*100
#fi['Mismatches [%]']=(final['0']+final['s'])/(final['0']+final['s']).sum()*100

fi['Mismatches [%]']=final['0']/final['0'].sum()*100

#fin['Gaps [%]']=final['g']/final['g'].sum()*100
#fi['Matches outside binding area [%]']=final['m']/final['m'].sum()*100
#fi['Mismatches outside binding area [%]']=final['s']/final['s'].sum()*100
print(fi.head(5))
fi=fi.reindex(p)
fi=fi.reset_index(drop=True)
print(fi)


# In[ ]:


p=list(range(33))
p=p[::-1]


# In[ ]:


#pal='rgbkymc'
p=fi.plot( alpha=0.4, rot=90, title='Alignments', linewidth=3.0)#,color=pal)
p.set_ylabel("Matches/Mismatches and Gaps [%]", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()
#axes.set_ylim([0,110])
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
#plt.savefig('/data1/MM_M_new_outsidemadded.pdf', bbox_inches='tight')

plt.show()
plt.close()

#pal = sns.color_palette("tab20c")
#pal=["r","g","b"]
pal='rgbkymc'
p=fi.plot.bar(stacked=True, alpha=0.4, rot=90, title='Alignments',color=pal)
p.set_ylabel("Matches/Mismatches and Gaps [%]", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()
#axes.set_ylim([0,110])
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
fig = plt.gcf()
fig.set_size_inches(8.5, 5)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
#plt.savefig('/data1/MM_M_Gap.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


#pal='rgbkymc'
p=fin.plot( alpha=0.4, rot=90, title='Alignments', linewidth=3.0)#,color=pal)
p.set_ylabel("Matches/Mismatches and Gaps [%]", fontsize=14, labelpad=10)
p.set_xlabel("Position in piRNA", fontsize=14, labelpad=10)
axes = plt.gca()
axes.set_ylim([0,15])
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
#plt.savefig('/data1/MM_M_OutisideBDarea_PERC.pdf', bbox_inches='tight')

plt.show()
plt.close()


# In[ ]:


def mef():
    """
    reads rna cofold output 
    """
    #sed -i 's/( /(/g' /data1/eCLASH//Vienna/target_pirnas.ps
    path='/data1/eCLASH/7-Exon-chimeras/thresh/Vienna/tab_output_RNAcofold.txt'
    df=pd.read_csv(path, header=None, sep='\t')
    df=df.round(0)
    df=pd.DataFrame(df[0].value_counts())
    #df=df.transpose()
    df=df.sort_index()
    df['Frequencies Chimeras [%]']=(df[0]/df[0].sum())*100
    
    path_s='/data1/eCLASH/tab_output_shuffled_RNAcofold.txt'
    d=pd.read_csv(path_s, header=None, sep='\t')
    d=d.round(0)
    d=pd.DataFrame(d[0].value_counts())
    #df=df.transpose()
    d=d.sort_index()
    d['Frequencies Shuffled [%]']=(d[0]/d[0].sum())*100
    print(df.head(5))
    print(d.head(5))
    
    
    
    p=df['Frequencies Chimeras [%]'].plot(kind= 'line' ,alpha=0.4, rot=70, 
                                 title='Alignments reverse complement',color='red')
    p=d['Frequencies Shuffled [%]'].plot(kind= 'line' ,alpha=0.4, rot=70, 
                                 title='Alignments reverse complement',color='blue')
    p.set_ylabel("Frequency [%]", fontsize=14, labelpad=10)
    #p.set_xlabel("r'\delta$'G [kcal/mol]", fontsize=14, labelpad=10)
    p.set_xlabel("$\Delta$G [kcal/mol]", fontsize=14, labelpad=10)
    axes = plt.gca()
    #axes.set_ylim([0,110])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.title("Predicted binding energies")
    leg = plt.legend()
    # get the lines and texts inside legend box
    leg_lines = leg.get_lines()
    leg_texts = leg.get_texts()
    # bulk-set the properties of all lines and texts
    plt.setp(leg_lines, linewidth=4)
    plt.setp(leg_texts, fontsize=12)
    plt.legend( bbox_to_anchor=(0.73,0.79), fontsize=14) #loc='center left',
    fig = plt.gcf()
    fig.set_size_inches(6, 5)
    plt.savefig('/data1/Binding_energies_chimera_shuffled.pdf', bbox_inches='tight')

    plt.show()
    plt.close()
mef()


# In[ ]:




