#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar
# In[ ]:


import subprocess, os, sys
import re
input_path= input('Input Path :')
name_of_outputfolder= input('Name of output folder :')


# In[ ]:


def get_info(input_path: str)-> list: 
    """
    files with .log extension are checked for motifs
        motif 1: percentage of rRNA (and with 100-percentage of rRNA -> percentage of RNA)
        motif 2: name of file -> Nr_123
    should be run in path of output files of trimming 
    """
    os.chdir('/data1/eCLASH/1-1-sortmerna/rRNA/')
    os.getcwd()
    #----------------------------------------------------------
    list_files=os.popen('ls').readlines()
    keep_phrases = [ "Total reads passing E-value threshold",
                                       "Reads file = " ]
    global keep_motifs
    keep_motifs=[]
    log_files=[]
    #----------------------------------------------------------
    for i in range(len(list_files)): #goes through the files in the listed files
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.log') and file not in log_files:
            log_files.append(file)
    for d in range(len(log_files)): #goes into all files
        with open(log_files[d]) as f: #opens the files
                    f = f.readlines() #reads each lines
                    motifs_list=[]
                    for line in f:
                        #print(line)
                        for phrase in keep_phrases:
                            if phrase in line:
                                #print(line, phrase)
                                motif_name=re.findall(r'([N][r]\d{2}_\w{4,5}_\d+)', line)
                                motif_percentage=re.findall(r'\d{0,2}\.\d{0,2}', line)
                                #print(motif_name)
                                #print(motif_percentage)
                                if len(motif_name)==1 and motif_name not in motifs_list:
                                    motifs_list.append(motif_name[0])
                                else:
                                    pass
                                if len(motif_percentage)==1:
                                    motif_p=motif_percentage[0]
                                    #print(motif_p)
                                    motifs_list.append(motif_p)
                                    this_one=motifs_list[0]
                                    for sublist in motifs_list:
                                        if sublist == this_one:
                                            pass
                                        else:
                                            keep_motifs.append(motifs_list)
                                else:
                                    pass    
        #break
    return keep_motifs

    
def export_info(input_path: str, name_of_outputfolder: str, keep_motifs: list):
    """
    1. change path -> output path, new dir is created
        name needs to be inputed from user 
    2. creates file (txt)
        input is the content of the list generated in prior function
    """
    #---------------------------------------------------------
    if not os.path.exists(input_path+'/'+name_of_outputfolder):
        os.makedirs(input_path+'/'+name_of_outputfolder)
    os.chdir(input_path+'/'+name_of_outputfolder)
    os.getcwd()
    #---------------------------------------------------------
    with open("summary_rRNA_RNA.txt", "w") as f: #opens//generates txt file
        for listl in keep_motifs: #iteration needs to be under the with open function
            f.write('rRNA\t{}\t{}\nRNA\t{}\t{}\n\n'.format(float(listl[1]), 
                                listl[0],
                                100-float(listl[1]), 
                                listl[0]
                                ))
    f.close() #close AFTER iteration 
    
    
def run_them(input_path, name_of_outputfolder):
    """
    run the functions 
    """
    get_info(input_path)
    export_info(input_path, name_of_outputfolder, keep_motifs)
    
run_them(input_path, name_of_outputfolder)


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import subprocess, os, sys
import pandas as pd 
import re
from glob import glob 
from zipfile import ZipFile


# In[ ]:


def make_stacked_bp():
    """
    from df make stacked bar plot 
    """
    
    df=pd.DataFrame.from_dict(ppp, orient='index')
    df=df.apply(pd.to_numeric)
    df['rRNA']=df[0]*60
    #df[]
    df['Other']=100-df['rRNA']
    
    df=df.drop([0], axis=1)
    #df=df.transpose()
    print(df)
    
    pal = sns.color_palette("tab20")
    p=df.plot.bar(stacked=True , color=pal, alpha=0.4, rot=35, title='rRNA content')
    p.set_ylabel("[%]", fontsize=12)
    p.set_xlabel("Replicate", fontsize=12)
    axes = plt.gca()
    axes.set_ylim([0,110])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('/data1/eCLASH/1-1-sortmerna/rRNA/output/rRNA_plot.pdf', 
                #bbox_inches='tight')
    plt.show()
    plt.close
make_stacked_bp()


# In[ ]:




