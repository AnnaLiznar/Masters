#!/usr/bin/env python
#written by Anna Liznar


import subprocess, os, sys
import re
import argparse
import operator
import numpy as np
import pandas as pd
import multiprocessing

import matplotlib.pyplot as plt


def main():
    args = get_args()
    num_cores = 7
    pool = multiprocessing.Pool(processes=num_cores)

    df4 = calc(args.infile, args.min_piRNAs, args.outfile)

########## END MAIN ###########


def calc(small, mind, out):
    """
    read in bed file as df
    count piRNAs in clustering groups
    
    make doctionary 
        use threshold for amount of piRNAs in groups 
        keep in df only those from above step 
        
    groupby of df by groups in dictionary 
    calc the length of stretch in df 
        include the length, gene, position, strand, contigue and start stop of stretch 
    """
    df4=pd.read_csv(sep='\t', filepath_or_buffer=small, header = None)
    elements=df4[19].tolist() #make list of all elements of column 19 -> clustered groups
    elements=list(set(elements)) #order list then make list 
    
    dd={} #initiate dic
    for el in elements: #iterate over list of all elements 
        dd[el]={'Amount':df4[19].value_counts()[el]} #count the occurence of groups and update the dic
        
    new={} #initiate dic 
    for key, value in dd.items(): #iterate over dic
        for k, v in value.items(): #ietarte over dic in dic (nested dic)
            if v > mind and key not in new: #threshhold for minimum of piRNAs
                new[key]={k:v} #update new dic
                
    keys=list(new.keys()) #new keys after threshold
    rslt = df4[df4[19].isin(keys)] #take only those rows, which are in new 
    
    kk=rslt.groupby(by=19, as_index=False).nth([0,-1]) #groupby 
    #keep firat and last of group 
    groups_list=list(range(0,kk.shape[0])) #shape : length of df make list
    for i in groups_list: #iterate
        if i%2==0: #if even 
            if kk.iloc[i,19] in new: #if group in new 

                new[kk.iloc[i,19]].update({'Length':kk.iloc[i+1,2]-kk.iloc[i,1]})
                #substract stop of last piRNA from start of first piRNA -> length
                new[kk.iloc[i,19]].update({'Gene':kk.iloc[i,11]})
                #gene annotation
                new[kk.iloc[i,19]].update({'Contig':kk.iloc[i,0]})
                #contigue
                new[kk.iloc[i,19]].update({'Start':kk.iloc[i,1]})
                #start of first piRNA
                new[kk.iloc[i,19]].update({'Stop':kk.iloc[i+1,2]})
                #stop of last piRNA
                new[kk.iloc[i,19]].update({'Strand':kk.iloc[i+1,5]})
                #strand of piRNA
                
    df = pd.DataFrame.from_dict(new, orient='index')
    #make df from dic
    columnsTitles = ['Contig', 'Start', 'Stop', 'Strand', 'Gene', 'Length', 'Amount']

    df = df.reindex(columns=columnsTitles)
    df[['Start', 'Stop', 'Length', 'Amount']] = df[['Start', 'Stop', 'Length', 'Amount']].astype(int)
    df.to_csv(path_or_buf=out, sep = '\t', header=True, index=False)

# Get command line arguments
def get_args():
    parser = argparse.ArgumentParser(description="Find high piRNA densities")
    parser.add_argument('--infile', type = argparse.FileType('r'), help = 'Input bedtools cluster bed file')
    parser.add_argument('--min_piRNAs', type = int, help =  'Minimum amount of piRNAs in stretch')
    parser.add_argument('--outfile', type = argparse.FileType('w'), help = 'Output bedfile filename')
    args = parser.parse_args()

    return args

# run main
if __name__ == "__main__":
    main()

