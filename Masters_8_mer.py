#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


def find(dic_all):
    """
    take dic from read_fasta()
    -looks for stretches of identies 
    -looks for max stretch 
    """
    dic=dic_all
    dic_instances={}
    for key, value in dic.items():
        for t, p in value.items():
            if t=='piRNA':
                temp=p
            if t=='target':
                instances=[]
                for i in range(len(p)):
                    if i<len(p)-11:
                        if p[i]=='N':
                            pass
                        else:
                            instances.append(Seq(p[i]+p[i+1]+p[i+2]+p[i+3])+p[i+4]+p[i+5]+p[i+6]+p[i+7]+p[i+8])

                m = motifs.create(instances)
                r = m.reverse_complement()
                for pos, seq in r.instances.search(temp):
                    if key not in dic_instances.keys():
                        dic_instances[key]=[pos]
                    else:
                        dic_instances[key].append(pos)

    df=pd.DataFrame.from_dict(dic_instances, orient='index')#.T

    columns=list(df.columns)

    maximum=max(columns)
    length=len(columns)
    for i in columns:
        if maximum <=length*4:
            df[maximum+i+1]=df[i]+1
            df[maximum+i+2]=df[i]+2
            df[maximum+i+3]=df[i]+3
            df[maximum+i+4]=df[i]+4
            df[maximum+i+5]=df[i]+5
            df[maximum+i+6]=df[i]+6
            df[maximum+i+7]=df[i]+7
            df[maximum+i+8]=df[i]+8
            columns=list(df.columns)
            maximum=max(columns)

    plt.plot(df, 'ro')
    plt.ylabel("Position in target", fontsize=12)
    plt.xlabel("piRNAs", fontsize=12)
    axes = plt.gca()

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #plt.savefig('/thresh/8-mer.pdf', bbox_inches='tight')
    plt.show()
    plt.close()
#find()


# In[6]:


def read_fasta():
    """
    opens fasta from --target --pirna 
    links them by ID in dictionary
    """
    
    path="/thresh/"
    os.chdir(path)
    large=glob('large*') #targets
    small=glob('small*') #piRNAs

    dic_all={}
    for i in range(len(large)):
        with open (path+'/small'+large[i].split('large')[1]) as fasta: #read fasta
            for title, sequence in SimpleFastaParser(fasta): #fill dictionary
                dic_all[title.split(None, 1)[0]]={'piRNA':sequence}
        with open(path+'/'+large[i]) as fasta_file:  # Will close handle cleanly
            for title, sequence in SimpleFastaParser(fasta_file):
                dic_all[title.split(None, 1)[0]].update({'target':sequence}) #update dic

    find(dic_all)

            
read_fasta()    


# In[ ]:




