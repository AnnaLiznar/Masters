#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar in jupyter
# In[ ]:


import pandas as pd 
import subprocess,sys, os 
import re
from glob import glob 


def make_2_col():
    """
    take planarian annotated smes file (new annotation of genome)
        has 10 columns and grep 0, 1, 2, 3 and 9 (zero based)
        column 9 has to be reduced to modified: 
        the whole column:
        transcript_id "SMEST026673001.1"; gene_id "SMESG000026673.1";
        this must stay:
        SMEST026673001.1
    """
    os.chdir('/data1/')
    path=os.getcwd()
    file=glob('*.bed')
    print(file)
    df=pd.read_csv(path+'/'+file[0], sep='\t', header=None, engine='python', usecols=[0,1,2,3,9])
    #print(df.head(5))
    df[9]=df[9].replace('transcript_id "', '', regex=True)
    #print(df.head(5))
    df[9]=df[9].replace('"; gene_id "SMESG[0-9]{8,9}.[0-9]{1}";', '', regex=True)
    print(df.head(20))
    output_path='/data1/Master/'
    
    #df.to_csv(output_path +'/smes_g.bed', sep='\t', header=False, index=False)


make_2_col()


# In[ ]:


def overlap_with_pandas():
    """
    get tsv file from interpro and smes_g.bed 
    find overlaps leave 2nd column of interpro in bed file 
    """
    
    print(os.getcwd())
    os.chdir('/data1/Master/')
    file=glob('*.tsv')
    df_ip=pd.read_csv('/data1/Master/'+file[0], sep='\t', 
                      header=None, engine='python')
    df_ip.rename(columns={0:'SMEST',1:'Interpro'}, inplace=True)
    
    df_col=pd.read_csv('/data1/Master/smes_g.bed', 
                      sep='\t', header=None, engine='python')
    df_col.rename(columns={0:'SMESG', 1:'SMEST'}, inplace=True)
    

    
    df_m=df_col.append(df_ip)

    
    col={'SMESG': 'first', 'Interpro':'first'}#,'sequences':'first'}

    df_g=df_m.groupby(df_m['SMEST'], sort=True).aggregate(col)
    print(df_g.head(20))
    
    
    #df_g.to_csv('/data1/Master/CLASH/Interpro.bed', 
     #             sep='\t', header=False, index=True)
        
    
    drop=df_g.dropna()
    print(drop.head(5))
    
    drop.to_csv('/data1/Master/for_heatmap.csv', 
                header=True, sep='\t', index=True)
    
overlap_with_pandas()

