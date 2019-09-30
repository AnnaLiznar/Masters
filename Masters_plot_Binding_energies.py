#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import motifs
from Bio.SeqIO.FastaIO import SimpleFastaParser
from glob import glob
import matplotlib.pyplot as plt


# In[2]:


def mef():
    """
    read ps file 
    """
    #sed -i 's/( /(/g' /data1/eCLASH/7-Exon-chimeras/thresh/Vienna/drei_target_pirnas.ps
    #before reading into python 
    #sed exchange '( ' with '(' ti enable reading into python as dataframe
    
    #actual chimeras
    path='/path/tab_output_RNAcofold.txt'
    df=pd.read_csv(path, header=None, sep='\t') #read into pandas df
    df=df.round(0) #round 
    df=pd.DataFrame(df[0].value_counts()) #count values (kind of make pivot)

    df=df.sort_index() #sort index
    df['Frequencies Chimeras [%]']=(df[0]/df[0].sum())*100 #calc percentages
    
    #shuffled
    path_s='/path/tab_output_shuffled_RNAcofold.txt'
    d=pd.read_csv(path_s, header=None, sep='\t')
    d=d.round(0)
    d=pd.DataFrame(d[0].value_counts())
    d=d.sort_index()
    d['Frequencies Shuffled [%]']=(d[0]/d[0].sum())*100

    
    p=df['Frequencies Chimeras [%]'].plot(kind= 'line' ,alpha=0.4, rot=70, 
                                 title='Alignments reverse complement',color='red')
    p=d['Frequencies Shuffled [%]'].plot(kind= 'line' ,alpha=0.4, rot=70, 
                                 title='Alignments reverse complement',color='blue')
    p.set_ylabel("Frequency [%]", fontsize=14, labelpad=10)

    p.set_xlabel("$\Delta$G [kcal/mol]", fontsize=14, labelpad=10) #make delta appear in xlabel
    axes = plt.gca()

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
    #plt.savefig('/path/Binding_energies_chimera_shuffled.pdf', bbox_inches='tight')

    plt.show()
    plt.close()
mef()


# In[ ]:




