#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#written by Anna Liznar in jupyter


# In[ ]:


def calc_GC(filepath):
    """
    calculate GCAT content of sequences 
 
    """
    liste=['small.exon.piRNA_2.fa', 'small.exon.piRNA_1.fa', 'small.exon.piRNA_3.fa']
    
    length=list(range(0,34))
    d={}
    for i in length:
        d[i]={'A':0, 'G':0, 'T':0, 'C':0}
    for i in liste:
        with open(filepath+'/'+i, 'r') as f:
            for line in f:
                #fasta header starts with >
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
    plt.savefig('/path/piRNAs_stacked_seq_logo.pdf', bbox_inches='tight')
    plt.show()
    plt.close()


    
calc_GC('/path/')


# In[ ]:


def find_U_A(filepath):
    """
    find names of sequences ahvin U bias and 10A bias
    """
    liste=['small.exon.piRNA_2.fa', 'small.exon.piRNA_1.fa', 'small.exon.piRNA_3.fa']
    

    d={'10A':{' '}, '5T':{' '}, '10ASeq':{' '}, '5TSeq':{' '}}
        
    print(d)
    for i in liste:
        with open(filepath+'/'+i, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    temp=line.split('\n')[0]

                elif line.startswith('T'):
                    d['5T'].update({temp.split('>')[1]})
                    d['5TSeq'].update({line.split('\n')[0]})
                elif line.startswith('A') or line.startswith('G') or line.startswith('C'):
                    d['10A'].update({temp.split('>')[1]})
                    d['10ASeq'].update({line.split('\n')[0]})

    
    df=pd.DataFrame.from_dict(d, orient='index')
    df=df.transpose()
    df.index = np.arange(1, len(df) + 1)
    df_save=pd.DataFrame()
    df_save=df['10A']
    print(df_save)
    
    df_save.to_csv('/path/Bias/10A.csv', sep='\t', header=False, index=False)
    df_5=df['5T']
    df_5.to_csv('/path/Bias/5T.csv', sep='\t', header=False, index=False)


    
find_U_A('/data1/eCLASH/7-Exon-chimeras/thresh/')

