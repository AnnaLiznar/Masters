#!/usr/bin/env python
# coding: utf-8

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


def find(dic_all, ik):
    """
    take dic from read_fasta()
    -looks for stretches of identies 
    -looks for max stretch 
    """
    dic=dic_all
    dic_instances={}
    it=[2,8,14]
    ttemp=pd.DataFrame()
    count=0
    for key, value in dic.items():
        for t, p in value.items():
            if t=='piRNA':
                temp=p
            if t=='target':
                instances=[]
                for i in range(len(p)):
                    if i<len(p)-(ik+1):
                        if p[i]=='N': #error
                            pass
                        else:
                            kk=i+ik
                            count+=1
                            instances.append(Seq(p[i:kk])) #make instances
                    
                    m = motifs.create(instances)
                    r = m.reverse_complement() #reverse 
                    ll=0
                    for pos, seq in r.instances.search(temp): #search fro perfect complements
                        if key not in dic_instances.keys():
                            dic_instances[key]=pos+1
                        else:
                            dic_instances[key+str(ll)]=pos+1
                            ll+=1
    tt=pd.DataFrame.from_dict(dic_instances, orient='index')#
    tt['{}-mer'.format(ik)]=tt[0]
    tt=tt.drop(axis=1, labels=0)
    ttemp=pd.concat([ttemp, tt], axis=1, sort=False)

    print('instances:', count)
    print(len(dic_instances))
    return ttemp


# In[ ]:


def plot(ttemp):
    """
    
    """
    df=ttemp.reset_index(drop=True)
    columns=list(df.columns)
    #pivot 
    final=pd.DataFrame()
    columns=list(df.columns)
    for c in columns:
        dd=pd.DataFrame(df[c].value_counts())
        final=final.append(dd, sort=True)
    
    ###merge by index
    final=final.groupby(level=0).sum()
    idx=list(range(0,33))
    
    final=final.reindex(idx, fill_value=np.nan)
    final=final.fillna(value=0)
    
    k=final.plot.line(alpha=0.4, rot=90, title='Perfect complementarity')
    plt.ylabel("Frequency of k-mer sliding window", fontsize=12)
    plt.xlabel("Position in piRNA", fontsize=12)
    axes = plt.gca()
    axes.set_xlim([-1, 31])
    plt.xticks(fontsize=12, rotation=80)
    plt.yticks(fontsize=12)

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
    
    #plt.savefig('/path/k_mer_sliding_window.pdf', bbox_inches='tight')
   
    plt.show()
    plt.close()

#find()


# In[ ]:


def read_fasta():
    """
    opens fasta from --target --pirna 
    links them by ID in dictionary
    """
    
    path="/path/thresh/"
    os.chdir(path)
    large=glob('large*') #targets
    small=glob('small*') #piRNAs

    dic_all={} #new dic
    count=0
    c=0
    ikk=[2,4,6,10,12] #kmer sizes
    df=pd.DataFrame()
    for i in range(len(large)):
        i=i+1
        with open (path+'/small'+large[i].split('large')[1]) as fasta:
            for title, sequence in SimpleFastaParser(fasta):
                dic_all[title.split(None, 1)[0]]={'piRNA':sequence}
                count+=1
        with open(path+'/'+large[i]) as fasta_file:  # Will close handle cleanly
            for title, sequence in SimpleFastaParser(fasta_file):
                dic_all[title.split(None, 1)[0]].update({'target':sequence})
                c+=1
        for ik in ikk:
            d=find(dic_all, ik)
            df=df.join(d, how='outer')
        plot(df)
        
        break
    print ('count: ',count, 'c: ', c)

read_fasta()    

