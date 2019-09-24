#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar in jupyter
# In[ ]:


import subprocess, os, sys

path_to_piPipes=input('Path to piPipes :') 
input_file_path=input('Input File Path :') 

output_dir_name=input('Output Dir Name :') 


def get_fq_files_from_output_dir(path_to_piPIpes, input_file_path, output_dir_name):
    """
    get files from dir
    run piPipes_fastq_to_insert on shell 
    """
    re_re, re_re_namesss=[], []
    #------------------------------------
    if not os.path.exists(input_file_path+'/'+output_dir_name):
        os.mkdir(input_file_path+'/'+output_dir_name)
    os.chdir(input_file_path)
    os.getcwd()
    #------------------------------------
    list_files=os.popen('ls').readlines()
    #------------------------------------
    for i in range(len(list_files)): #goes through the files in the listed files
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.startswith('rejected_reads') and file not in re_re:
            re_re.append(file)
    for jo in re_re:
        if jo.endswith('.fq') or jo.endswith('.fq.fq'):
            try:
                file=jo.split('.fq')#looses the .fq
            except:
                print('extension .fq not splittable')
            try:
                file=jo.split('.fq.fq') #looses the .fq.fq
            except:
                print('extension .fq.fq not splittable')
            file=file[0]
            if jo.startswith('rejected_reads') and file not in re_re_namesss:
                re_re_namesss.append(file)
    os.chdir(path_to_piPipes+'/')
    os.getcwd()            
    try:
        for r in range(len(re_re_namesss)):
            command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq',
                input_file_path+'/'+output_dir_name+'/'+re_re_namesss[r]+'.insert'
                                                            )
            subprocess.check_call(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            print(command)
    except:
        command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq',
                input_file_path+'/'+output_dir_name+'/'+re_re_namesss[r]+'.insert'
                                                            )
        print('{} command couldnt be run'.format(command))
    try:
        for r in range(len(re_re_namesss)):
            command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq.fq',
                input_file_path+'/'+output_dir_name+'/'+re_re_namesss[r]+'.insert'
                                                            )
            subprocess.check_call(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            print(command)
    except:
        command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq.fq',
                input_file_path+'/'+output_dir_name+'/'+re_re_namesss[r]+'.insert'
                                                            )
        print('{} command couldnt be run'.format(command))

    #print(re_re, '\n',re_re_namesss)
get_fq_files_from_output_dir(path_to_piPipes, input_file_path, output_dir_name)


# In[ ]:




