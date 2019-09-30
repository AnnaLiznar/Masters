#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar in jupyter
# In[ ]:


import subprocess, os, sys
import re
import optparse, shutil


# In[ ]:


path_to_input_fasta=input('Path to input bed2 files : ') 

output_path_RNA=input('Path to output w/ new dir name : ')

def exchange_U_w_T(path_to_input_fasta, output_path_RNA):
    """
    exchanges Uridin with Thymin --> DNA to RNA
    """
    fasta_files=[]
    #----------------------------------------------------------
    if not os.path.exists(output_path_RNA): #make dir for output
        os.mkdir(output_path_RNA)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_fasta)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #------------------------------------
    for i in range(len(list_files)): #goes through the files in the listed files
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.fa') and file not in fasta_files:
            fasta_files.append(file)

   
    for file in fasta_files:
            op_file=file.split('25_cutoff_') #looses the \n= newline
            op_file=op_file[1]
            #print(op_file)            
            command="sed 's/T/U/g' {} > {}/U_RNA_{}".format(
                    path_to_input_fasta+file, output_path_RNA, op_file
                                            )

            subprocess.check_call(command, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, shell=True)
            print(command)   
        
            
exchange_U_w_T(path_to_input_fasta, output_path_RNA)


# In[ ]:




