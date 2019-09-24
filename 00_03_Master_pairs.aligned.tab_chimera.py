#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd 
import subprocess, sys, os 
import re
from glob import glob 


# In[ ]:


def get_chimera():
    """
    chimera: sequence made of two distinct origins 
    in annotated bedfiles all found chimeras were used 
    in pairs.aligned.tab file are chimera names are stored, which actually are chimeras 
    -> make list of all chimeras from pairs.aligned.tab file 
    """
    pairs_aligned_path='/data1/Master/CLASH/Analysis/6-chimeraSearch_clipseqtool/results/'
    dir_list=['1_1u', '3_2u', '5_3u', '7_4u', '9_5u'] #iterate 
    pairs_aligned_file='pairs.aligned.tab'
    
    
    for ordner in dir_list:
        os.chdir(pairs_aligned_path+'/'+ordner)
        print(os.getcwd())
        path=os.getcwd()
        file=glob(pairs_aligned_file)
        print(file)
        command_cat='cat {} | cut -f1 > {}/{}_input_names.txt'.format(
                    path+'/'+file[0], path, ordner)
        print(command_cat)
        subprocess.check_call(command_cat, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
get_chimera()


# In[ ]:


def make_bedfile():
    """
    -> grep only those chimeric sequences which are stored in list from get_chimera
    """
    bedtools_path='/data1/Master/CLASH/Analysis/7-bedtools/'
    elements=['LINE', 'DNA', 'LTR', 'Unk', 'Intron', 'No_annot', 'Cluster']
    #path and list for input names list
    pairs_aligned_path='/data1/Master/CLASH/Analysis/6-chimeraSearch_clipseqtool/results/'
    dir_list=['1_1u', '3_2u', '5_3u', '7_4u', '9_5u']
    input_names_list=[]
    for ordner in dir_list:
        input_names_list.append(pairs_aligned_path+ordner+'/'+ordner+'_input_names.txt')
    #print(input_names_list)
    
    for ele in elements:
        os.chdir(bedtools_path+'/'+ele)
        files=glob(ele+'*.bed')
        for file in files:
            name_one=file.split('fragments.large.filtered.Aligned.bed')[0]
            if 'align' in name_one:
                name=name_one.split(ele+'_align_')[1]
            else:
                name=name_one.split(ele+'_')[1]
            #print(name)


            for ink in input_names_list:
                if name in ink:
                    comd='grep -Fwf {} {} > {}'.format(
                            ink, 
                            bedtools_path+ele+'/'+file, 
                            bedtools_path+ele+'/grep_chimera_'+file)
                    try:
                        subprocess.check_call(comd, 
                                              stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE, 
                                              shell=True)
                    except subprocess.CalledProcessError as e:
                        if e.returncode > 1:
                            raise
                    print(comd)
                else:
                    continue
    
make_bedfile()  

