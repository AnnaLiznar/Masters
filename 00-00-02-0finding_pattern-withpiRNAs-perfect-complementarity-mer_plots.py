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


def find(dic_all, ik):
    """
    take dic from read_fasta()
    -looks for stretches of identies 
    -looks for max stretch 
    """
    dic=dic_all
    dic_instances={}
    #it=list(range(2,14,4))
    it=[2,8,14]
    #idxe=list(range(0,33))
    ttemp=pd.DataFrame()
    #ttemp=ttemp.reindex(idxe)
    count=0
    for key, value in dic.items():
        #print (key, value)
        for t, p in value.items():
            #print(t, p)
            if t=='piRNA':
                #dic_instances[key]={t:p}
                temp=p
            if t=='target':
                instances=[]
                for i in range(len(p)):
                    if i<len(p)-(ik+1):
                        if p[i]=='N':
                            pass
                        else:
                            kk=i+ik
                            #print(Seq(p[i:kk]))
                            count+=1
                            #print('i:kk', i, kk)
                            instances.append(Seq(p[i:kk]))
#Seq(p[i]+p[i+1]+p[i+2]+p[i+3])+p[i+4])#+p[i+5]+p[i+6]+p[i+7]+p[i+8]+p[i+9]+p[i+10]+p[i+11]+p[i+12]+p[i+13])
                    
                    m = motifs.create(instances)
                    r = m.reverse_complement() #reverse 
                    #r.degenerate_consensus
                    ll=0
                    for pos, seq in r.instances.search(temp):
                        if key not in dic_instances.keys():
                            dic_instances[key]=pos+1
                            #print(dic_instances)
                        else:
                            dic_instances[key+str(ll)]=pos+1
                            ll+=1
    tt=pd.DataFrame.from_dict(dic_instances, orient='index')#.T
    tt['{}-mer'.format(ik)]=tt[0]
    tt=tt.drop(axis=1, labels=0)
    ttemp=pd.concat([ttemp, tt], axis=1, sort=False)
    #ttemp=ttemp.join(tt, how='outer')
    print('instances:', count)
    print(len(dic_instances))
    return ttemp


# In[ ]:


def plot(ttemp):
    """
    
    """
    #df=pd.DataFrame.from_dict(dic_instances, orient='index')#.T
    #df=df.astype(int)
    #print(df)
    df=ttemp.reset_index(drop=True)
    #print(df)
    columns=list(df.columns)
    #print(columns)
    #pivot 
    final=pd.DataFrame()
    columns=list(df.columns)
    for c in columns:
        dd=pd.DataFrame(df[c].value_counts())
        #dd=dd.transpose()
        final=final.append(dd, sort=True)
    #print(final)
    
    ###merge by index
    final=final.groupby(level=0).sum()
    idx=list(range(0,33))
    
    final=final.reindex(idx, fill_value=np.nan)
    #print(final)
    final=final.fillna(value=0)
    print(final)
    
    k=final.plot.line(alpha=0.4, rot=90, title='Perfect complementarity')
    plt.ylabel("Frequency of k-mer sliding window", fontsize=12)
    plt.xlabel("Position in piRNA", fontsize=12)
    axes = plt.gca()
    #axes.set_ylim([-5,55])
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
    
    #plt.savefig('/data1/eCLASH/k_mer_sliding_window.pdf', bbox_inches='tight')
   
    plt.show()
    plt.close()

#find()


# In[ ]:


def read_fasta():
    """
    opens fasta from --target --pirna 
    links them by ID in dictionary
    """
    
    path="/data1/eCLASH/7-Exon-chimeras/thresh/"
    os.chdir(path)
    large=glob('large*')
    small=glob('small*')
    print(large)
    print(small)
    dic_all={}
    count=0
    c=0
    ikk=[2,4,6,10,12]
    df=pd.DataFrame()
    for i in range(len(large)):
        i=i+1
        print(large[i])
        print('/small'+large[i].split('large')[1])
        with open (path+'/small'+large[i].split('large')[1]) as fasta:
            for title, sequence in SimpleFastaParser(fasta):
                dic_all[title.split(None, 1)[0]]={'piRNA':sequence}
                count+=1
        with open(path+'/'+large[i]) as fasta_file:  # Will close handle cleanly
            for title, sequence in SimpleFastaParser(fasta_file):
                dic_all[title.split(None, 1)[0]].update({'target':sequence})
                #dic_elemets[name.split('.txt')[0]].update({mo:line}
                c+=1
        #print(dic_all)
        for ik in ikk:
            
            d=find(dic_all, ik)
            df=df.join(d, how='outer')
        
        plot(df)
        
        break
    print ('count: ',count, 'c: ', c)
    #return dic_all

            
read_fasta()    


# In[ ]:




