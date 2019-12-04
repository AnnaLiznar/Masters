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

import matplotlib.pyplot as plt


# In[ ]:


def main():
    args = get_args()
    num_cores = 7
    pool = multiprocessing.Pool(processes=num_cores)

    df = keep(args.infile, args.outfile)

########## END MAIN ###########


# In[ ]:


def keep(path, op):
    """
    keep those genes, which occur only once in df
        -> potential kd templates 
    """
    df=pd.read_csv(filepath_or_buffer=path, 
                  sep='\t') #read in file

    d=df['Gene'].value_counts().to_dict() #make dictionary with counts of genes
    keep=[] #inititate list
    for key, value in d.items(): #iterate over dic with counts 
        if value == 1: #keep only those which occur only onces 
              keep.append(key) #append to list 

    df=df[df['Gene'].isin(keep)]#.reset_index(drop=True)
    
    df.to_csv(path_or_buf=op, sep = '\t', header=True, index=False)


# In[19]:


# Get command line arguments
def get_args():
    parser = argparse.ArgumentParser(description="Find unique")
    parser.add_argument('--infile', type = argparse.FileType('r'), 
                        help = 'Input: custom bed file after clustering ')
    
    parser.add_argument('--outfile', type = argparse.FileType('w'), 
                        help = 'Output custom bedfile filename')
    args = parser.parse_args()

    return args


# In[ ]:


# run main
if __name__ == "__main__":
    main()

