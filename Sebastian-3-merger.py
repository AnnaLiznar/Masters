#!/usr/bin/env python
# coding: utf-8

# In[2]:


import subprocess, os, sys
import re
import argparse
import operator
import numpy as np
import pandas as pd
import multiprocessing
from glob import glob

import matplotlib.pyplot as plt


# In[ ]:


def main():
    args = get_args()
    num_cores = 7
    pool = multiprocessing.Pool(processes=num_cores)

    df = merger(args.infile_path, args.motif, args.outfile)

########## END MAIN ###########


# In[8]:


def merger(path, motif, op):
    """
    keep those which occur at least in 3 replicates 
    make new df with stretches
        -> the longest stretch from the 4 reps is saved rest is discarded and considered duplicates
    """
    #-------------------------
    os.chdir(path)
    os.getcwd()
    no=glob('*'+motif)
    print('files: ',no)
    #-------------------------
    di ={}
    #iterate over dfs and load gene ids in dictionary 
    for i in no:
        df=pd.read_csv(filepath_or_buffer=path+i,
                      sep= "\t")
        for ele in df['Gene']:
            if ele not in di:
                di[ele]={i:True}
            elif ele in di:
                di[ele].update({i:True})

    #make list which gene ids should be kept in df 
    #thresh is at least in 3 dfs
    keep=[]
    for key, value in di.items():
        if len(value)>=3:
            keep.append(key)
    new_df = pd.DataFrame({'Contig': [0],
                       'Start': [0],
                       'Stop':[0] , 
                       'Strand':[0], 
                       'Gene': [0], 
                       'Length': [0], 
                       'Amount':[0]})
    #print(keep)
    for gene in keep:
        #print(gene)
        for p in no:
            df=pd.read_csv(filepath_or_buffer=path+p,
                      sep= "\t")
            size=df.shape[0]
            kk=list(range(0, size))
            for i in kk:
                #maxi=df.loc[i,'Length']
                if gene in df.loc[i,'Gene']:
                        #if maxi <= new_df.loc[i,'Length']:
                    new_df=new_df.append({'Contig': df.loc[i,'Contig'],
                                          'Start': df.loc[i,'Start'], 
                                          'Stop':df.loc[i,'Stop'] ,
                                          'Strand':df.loc[i,'Strand'],
                                          'Gene': df.loc[i,'Gene'],
                                          'Length': df.loc[i,'Length'], 
                                          'Amount':df.loc[i,'Amount']},
                                          ignore_index=True)
                    new_df=new_df.sort_values('Length', 
                                              ascending=False).drop_duplicates('Gene').sort_index() 
                    #print(new_df)
    new_df=new_df.iloc[1:]
    new_df=new_df[new_df['Length'] > 1000]
    #return new_df
    new_df.to_csv(path_or_buf=op, sep = '\t', header=True, index=False)


# Get command line arguments
def get_args():
    parser = argparse.ArgumentParser(description="Find unique")
    parser.add_argument('--infile_path', type = str, 
                        help = 'path to infiles\ninfiles: custom bedfiles\nhere only path ')
    parser.add_argument('--motif', type = str, 
                       help = 'Motif for infile\nMotif: files ends with same string')
    parser.add_argument('--outfile', type = argparse.FileType('w'), 
                        help = 'Output custom bedfile filename')

    args = parser.parse_args()

    return args


# In[ ]:


# run main
if __name__ == "__main__":
    main()

