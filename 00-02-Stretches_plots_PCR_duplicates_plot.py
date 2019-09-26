#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar
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


# In[ ]:


p={'piRNAs':{'all':0, 'SM3':0},'length':0,'Strand':{'plus':0, 'minus':0}}
le={'length':0}
pm={'plus':0, 'minus':0}
X1={}
Xins={}

amount_X1=list(range(0, 18))
amount_Xins=list(range(0, 7))


# In[ ]:


for i in amount_X1:
    X1[i]=pm
    
for x in amount_Xins:
    Xins[x]=pm


# In[ ]:


Xins_all={0: {'plus': 3, 'minus': 25}, 
 1: {'plus': 0, 'minus': 20}, 
 2: {'plus': 0, 'minus': 5}, 
 3: {'plus': 0, 'minus': 2},
 4: {'plus': 0, 'minus': 1}, 
 5: {'plus': 0, 'minus': 1}}


# In[ ]:


Xins_sm3={0: {'plus': 1, 'minus': 22}, 
 1: {'plus': 0, 'minus': 18}, 
 2: {'plus': 0, 'minus': 4}, 
 3: {'plus': 0, 'minus': 1},
 4: {'plus': 0, 'minus': 1}, 
 5: {'plus': 0, 'minus': 1}}


# In[ ]:


df_Xins_sm3=pd.DataFrame.from_dict(Xins_sm3)
df_Xins_sm3=df_Xins_sm3.transpose()
df_Xins_sm3
df_Xins_all=pd.DataFrame.from_dict(Xins_all)
df_Xins_all=df_Xins_all.transpose()
df_Xins_all
df_Xins_strand=pd.concat([df_Xins_all, df_Xins_sm3], axis=1, sort=False, keys=['all', 'sm3'])
df_Xins_strand


# In[ ]:


df_Xins=pd.DataFrame.from_dict(Xins)
df_Xins=df_Xins.transpose()
df_Xins


# In[ ]:


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


# In[ ]:


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


df_X1_sm3=pd.DataFrame.from_dict(X1_SM3)
df_X1_sm3=df_X1_sm3.transpose()
df_X1_sm3


# In[ ]:


df_X1_all=pd.DataFrame.from_dict(X1_all)
df_X1_all=df_X1_all.transpose()
df_X1_all
df_strand=pd.concat([df_X1_all, df_X1_sm3], axis=1, sort=False, keys=['all', 'sm3'])
df_strand


# In[ ]:


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


# In[ ]:


Xins={0: {'allpiRNAs': 28, 'SM3': 23, 'length': 299},
 1: {'allpiRNAs': 20, 'SM3': 18, 'length': 273}, 
 2:  {'allpiRNAs': 5, 'SM3': 4, 'length': 40},
 3:  {'allpiRNAs': 2, 'SM3': 1, 'length': 30}, 
 4:  {'allpiRNAs': 1, 'SM3': 1, 'length': 31},
 5:  {'allpiRNAs': 1, 'SM3': 1, 'length': 32}}


# In[ ]:


df_X1=pd.DataFrame.from_dict(X1)
df_X1=df_X1.transpose()
print(df_X1['SM3'].sum(), df_X1['allpiRNAs'].sum())
df_X1


# In[ ]:


df_Xins=pd.DataFrame.from_dict(Xins)
df_Xins=df_Xins.transpose()
df_Xins


# In[ ]:

df_X1['length'].plot.bar(rot=0, color='b')

# assign your bars to a variable so their attributes can be accessed
bars = plt.bar(df_X1.index, height=df_X1['length'], width=.4)

# access the bar attributes to place the text in the appropriate location
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x()-0.1, yval+1, yval)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(2, 1)


plt.show()
plt.close()



for i in df_X1.columns:
    #df_X1[['SM3', 'allpiRNAs']].plot.bar(subplots=True,rot=0)


   
    df_X1[i].plot.bar(rot=0, color='r')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=12)
    fig = plt.gcf()
    fig.set_size_inches(3, 2)


    plt.show()
    plt.close()
    


# In[ ]:


df_strand.plot.bar(rot=0)


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(15, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/data1/Stretches/Strandness_X1.pdf', bbox_inches='tight')
plt.show()
plt.close()

from itertools import cycle, islice
my_colors = list(islice(cycle([ 'y', 'k']), None, len(df_X1[['allpiRNAs','SM3']])))   
#df_X1[['allpiRNAs','SM3']].plot.bar(rot=0)
df_X1[['allpiRNAs','SM3']].plot(kind='bar',  color=my_colors)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(9, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/data1/Stretches/Amount_piRNAs_X1.pdf', bbox_inches='tight')
plt.show()
plt.close()


df_X1['length'].plot.bar(rot=0, color='g')


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(9, 3)

plt.savefig('/data1/Stretches/Length_X1.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


df_Xins_strand.plot.bar(rot=0)


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(7, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/data1/Stretches/Strandness_Xins.pdf', bbox_inches='tight')
plt.show()
plt.close()

from itertools import cycle, islice
my_colors = list(islice(cycle([ 'y', 'k']), None, len(df_Xins[['allpiRNAs','SM3']])))   
#df_X1[['allpiRNAs','SM3']].plot.bar(rot=0)
df_Xins[['allpiRNAs','SM3']].plot(kind='bar',  color=my_colors)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(5, 3)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/data1/Stretches/Amount_piRNAs_Xins.pdf', bbox_inches='tight')
plt.show()
plt.close()


df_Xins['length'].plot.bar(rot=0, color='g')


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(5, 3)

plt.savefig('/data1/Stretches/Length_Xins.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


df_Xins.plot.bar(subplots=True, rot=0)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(5, 9)


plt.show()
plt.close()


# In[ ]:


p=df.plot.bar(alpha=0.4, rot=90, title='Minus',color=pal)
    p.set_ylabel("Length of stretch", fontsize=14, labelpad=10)
    p.set_xlabel("Stretches", fontsize=14, labelpad=10)
    axes = plt.gca()
    #axes.set_ylim([0,110])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(8.5, 3)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    plt.savefig('/data1/Stretches/6-bedtools/Exon/Merged/350/ForPlot/Length_stretches.pdf', bbox_inches='tight')

    plt.show()
    plt.close()


# In[ ]:


Dup={'ctrl_1':{'Duplicates':52738812, 'Single': 6313196},
    'ctrl_2':{'Duplicates':64744568, 'Single': 1247232},
    'piRNA_1':{'Duplicates':64154364, 'Single': 1610884},
    'piRNA_2':{'Duplicates':62097884, 'Single': 1571516},
    'piRNA_3':{'Duplicates':51345432, 'Single': 1189556}}
ddf=pd.DataFrame.from_dict(Dup)
ddf=ddf.transpose()
ddf['total']=ddf['Duplicates']+ddf['Single']
ddf['Duplicates [%]']=ddf['Duplicates']/ddf['total']*100
ddf['Single [%]']=ddf['Single']/ddf['total']*100
ddf


# In[ ]:


ddf[['Duplicates [%]','Single [%]']].plot.bar(stacked=True,rot=0)


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
plt.savefig('/data1/Figures_preliminary/PCR_duplicates.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


#potential targets 
ff={'ctrl_1':{'Potential':57529010, 'piRNA_partner':126852, 'Actual': 46},
    'ctrl_2':{'Potential':11124940, 'piRNA_partner':16263, 'Actual': 69},
    'piRNA_1':{'Potential':23750046, 'piRNA_partner':49842, 'Actual': 26143},
    'piRNA_2':{'Potential':21290694,  'piRNA_partner':61861, 'Actual': 18989},
    'piRNA_3':{'Potential':21149866, 'piRNA_partner':52530,  'Actual': 36086},}
podf=pd.DataFrame.from_dict(ff)
podf=podf.transpose()
podf['total']=podf['Actual']+podf['Potential']+podf['piRNA_partner']
podf['Actual [%]']=podf['Actual']/podf['total']*100
podf['Potential [%]']=podf['Potential']/podf['total']*100
podf['piRNA partner [%]']=podf['piRNA_partner']/podf['total']*100
podf


# In[ ]:


podf[['Potential [%]','Actual [%]','piRNA partner [%]']].plot.bar(stacked=True,rot=0)


plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=12)
axes = plt.gca()
axes.set_ylim([95,100.1])
#axes.set_xlim([-1, 31])
fig = plt.gcf()
fig.set_size_inches(6, 4)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
#plt.savefig('/data1/Figures_preliminary/Potential_targets.pdf', bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


def calc_GC(filepath):
    """
    calculate GC content of sequences 
    and AT
    """
    liste=['small.exon.piRNA_2.fa', 'small.exon.piRNA_1.fa', 'small.exon.piRNA_3.fa']
    
    length=list(range(0,34))
    d={}
    for i in length:
        d[i]={'A':0, 'G':0, 'T':0, 'C':0}
    for i in liste:
        with open(filepath+'/'+i, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    pass
                else:
                    line_l=list(line)
                    for el in range(len(line_l)):
                        if line_l[el]=='A':
                            d[el]['A']+=1
                        elif line_l[el]=='T':
                            d[el]['T']+=1
                        elif line_l[el]== 'G':
                            d[el]['G']+=1
                        elif line_l[el]== 'C':
                            d[el]['C']+=1

    df=pd.DataFrame.from_dict(d)
    df=df.transpose()
    df.index = np.arange(1, len(df) + 1)
    

    df['A [%]']=df['A']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    df['G [%]']=df['G']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    df['T [%]']=df['T']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    df['C [%]']=df['C']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    print (df['A'].sum())
    print (df)
    plot(df)
def plot(df):
    """
    plot
    """
    df[['A [%]','T [%]', 'G [%]',  'C [%]']].plot(kind='bar', stacked=True, rot=0)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=12)
    axes = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(18, 4)
    plt.savefig('/data1/Figures_preliminary/piRNAs_stacked_seq_logo.pdf', bbox_inches='tight')
    plt.show()
    plt.close()


    
calc_GC('/data1/eCLASH/7-Exon-chimeras/thresh/')


# In[ ]:


def find_U_A(filepath):
    """
    calculate GC content of sequences 
    and AT
    """
    liste=['small.exon.piRNA_2.fa', 'small.exon.piRNA_1.fa', 'small.exon.piRNA_3.fa']
    

    d={'10A':{' '}, '5T':{' '}, '10ASeq':{' '}, '5TSeq':{' '}}
        
    print(d)
    for i in liste:
        with open(filepath+'/'+i, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    temp=line.split('\n')[0]
                    #print(temp)
                    
                elif line.startswith('T'):
                    #print(temp,'line: ',line)
                    d['5T'].update({temp.split('>')[1]})
                    d['5TSeq'].update({line.split('\n')[0]})
                elif line.startswith('A') or line.startswith('G') or line.startswith('C'):
                    d['10A'].update({temp.split('>')[1]})
                    d['10ASeq'].update({line.split('\n')[0]})
                    
                    #line_l=list(line)
                    #for el in range(len(line_l)):
                        #if el==9 and line_l[el]=='A':
                            #d[el]['10A'].update(temp)
                        #elif el==0 and line_l[el]=='T':
                            #d[el]['5T'].update(temp)
        #break

    #print(d)
    
    df=pd.DataFrame.from_dict(d, orient='index')
    df=df.transpose()
    df.index = np.arange(1, len(df) + 1)
    df_save=pd.DataFrame()
    df_save=df['10A']
    print(df_save)
    
    df_save.to_csv('/data1/eCLASH/7-Exon-chimeras/thresh/Bias/10A.csv', sep='\t', header=False, index=False)
    df_5=df['5T']
    df_5.to_csv('/data1/eCLASH/7-Exon-chimeras/thresh/Bias/5T.csv', sep='\t', header=False, index=False)
    #df=df.drop[['10A', '5T']]

    #df['A [%]']=df['A']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    #df['G [%]']=df['G']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    #df['T [%]']=df['T']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    #df['C [%]']=df['C']/(df['A'].sum()+df['G'].sum()+df['C'].sum()+df['T'].sum())*100
    #print (df['A'].sum())
    #print (df)
    #plot(df)
def plot(df):
    """
    plot
    """
    df[['A [%]','T [%]', 'G [%]',  'C [%]']].plot(kind='bar', stacked=True, rot=0)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=12)
    axes = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(18, 4)
    #plt.savefig('/data1/Figures_preliminary/piRNAs_stacked_seq_logo.pdf', bbox_inches='tight')
    plt.show()
    plt.close()


    
find_U_A('/data1/eCLASH/7-Exon-chimeras/thresh/')


# In[ ]:


def make_df():
    """
    make df from mapping events full length piRNA and seed
    """
    seed={'X1':{'full':154741, 
                'mapped_full_once': 90096, 
                'mapped_full_once_3mism': 91843,
                'seed': 74900, 
                'mapped_seed_once': 90096},
        'Xins':{'full':114898, 
                'mapped_full_once': 66940, 
                'mapped_full_once_3mism': 68101,
                'seed': 65595, 
                'mapped_seed_once': 66940}
         }
    df_seed=pd.DataFrame.from_dict(seed)
    df_seed=df_seed.transpose()
    df_seed['mapped_full_once [%]']=df_seed['mapped_full_once']/df_seed['full']*100
    df_seed['mapped_seed_once [%]']=df_seed['mapped_seed_once']/df_seed['seed']*100
    df_seed['mapped_full_once_3mism [%]']=df_seed['mapped_full_once_3mism']/df_seed['full']*100
    return df_seed
    
make_df()


# In[ ]:


def plot_make_df():
    """
    plot stacked from df from make df 
    """
    df=make_df()
    df[['mapped_full_once [%]','mapped_full_once_3mism [%]', 'mapped_seed_once [%]']].plot(kind='bar', rot=0)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=12)
    axes = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(2, 5)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=14)
    plt.savefig('/data1/Seed/Barplot_seed_full_3mm.pdf', bbox_inches='tight')
    plt.show()
    plt.close()
    
plot_make_df()


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
    #print(df_seed)
    df_s=pd.DataFrame()
    df_s['piRNAs']=df_seed['X1']+df_seed['Xins']/2

    df_s=df_s.transpose()
    df_s['full 0 mm[%]']=df_s['full']/df_s['all']*100
    df_s['seed [%]']=df_s['seed']/df_s['all']*100
    df_s['full 1 mm [%]']=df_s['full_one_mm']/df_s['all']*100
    df_s['full 2 mm [%]']=df_s['full_two_mm']/df_s['all']*100
    df_s['full 3 mm [%]']=df_s['full_three_mm']/df_s['all']*100
    #print(df_s)
    print(df_s.loc['piRNAs','seed [%]'])
    #df_seed['mapped_full_once_3mism [%]']=df_seed['mapped_full_once_3mism']/df_seed['full']*100
    #seed == 0 #horizontal line 
    pl=pd.DataFrame()
    pl=df_s[['full 0 mm[%]','full 1 mm [%]', 'full 2 mm [%]', 'full 3 mm [%]']]
    pl=pl.transpose()
    #pl=pl-df_s.loc['piRNAs','seed [%]']
    #print(pl)
    return pl
    
make_dff()


# In[ ]:


def plot_df():
    """
    plot from df from make df 
    """
    df=make_dff()
    #df[['full [%]', 'seed [%]']].plot(kind='bar', rot=0)
    df.plot(kind='line', rot=0)
    plt.xticks(fontsize=12)
    
    plt.yticks(fontsize=12)
    plt.ylabel('Percentage of targets', fontsize=12)
    plt.legend(fontsize=12)
    axes = plt.gca()
    fig = plt.gcf()
    fig.set_size_inches(3, 2.5)
    axes.set_xticklabels(['0 MM ','0 MM ','1 MM', '2 MM', '3 MM'])
    #plt.hlines(y=57.19628695,xmin=0, xmax=df['full 3 mm [%]':'piRNAs'], colors='r')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('/data1/Seed/Line_new_0123.pdf', bbox_inches='tight')
    plt.show()
    plt.close()
    
plot_df()

