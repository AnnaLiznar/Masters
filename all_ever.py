#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess, os, sys
import ipywidgets as widgets
from IPython.display import display
import optparse, shutil

import cutadapt

from ipywidgets import Layout, HBox, VBox
from ipywidgets import interact, interactive, fixed, interact_manual
import gzip


def getDirectoryContent():
    """
    lists directory content
    """
    global test   
    test = False
    directories = ['../', './'] #go back to upper dir
    files = []
    for item in os.listdir(os.getcwd()): #
        if os.path.isdir(item):
            if item == '.ipynb_checkpoints':
                pass
            elif item[0] == '.':
                pass
            else:
                directories.append(item)
        else:
            files.append(item)     
    main.options =  directories
    main.value = './'
    sub.options = sorted(files)
    curDir.value = os.getcwd()

    test = True
    
    
def mainObserver(sender):
    '''Activates with observe from getDirectoryContent.
     If 'new' is chosen, it returns. 
    '''
    global test
    if test:
        if sender['type'] == 'change' and sender['name'] == 'value':
            try:
                os.chdir('./' + sender['new'])
                getDirectoryContent()
            except:
                raise
    else:
        pass
#-----------------------------------------------------------------------------
curDir = widgets.Label(value=os.getcwd())

print('current Directory:')

main = widgets.Select(description='Directories:')
sub = widgets.Select(options=['MMM'], value='MMM', 
                     description='Files:')
output_Dir= widgets.Text(
    placeholder='i.e. :/data1/AnLi_FM_MA/data/name',
        description='Output Path+name:', disabled=False)
#----outputdir
box = widgets.HBox([main, sub, output_Dir])             
#--->output_Dir.value

display(curDir, box)

test = True

getDirectoryContent()

main.observe(mainObserver)


# In[ ]:


print('input: ', curDir.value,'\n' 'output: ', output_Dir.value)


# In[ ]:


#vincent widgets -> for the path (input and output)
#go to directory where planarian rRNA full sequences are stored
#choose a name for dir for indexed rRNA 
def get_DirectoryContent():
    """
    lists directory content
    """
    global test   
    test = False
    directories = ['../', './'] #go back to upper dir
    files = []
    for item in os.listdir(os.getcwd()): #
        if os.path.isdir(item):
            if item == '.ipynb_checkpoints':
                pass
            elif item[0] == '.':
                pass
            else:
                directories.append(item)
        else:
            files.append(item)     
    main.options =  directories
    main.value = './'
    sub_rRNA.options = sorted(files)
    curDir_rRNA.value = os.getcwd()

    test = True
    
    
def main_Observer(sender):
    '''Activates with observe from get_DirectoryContent.
     If 'new' is chosen, it returns. 
    '''
    global test
    if test:
        if sender['type'] == 'change' and sender['name'] == 'value':
            try:
                os.chdir('./' + sender['new'])
                get_DirectoryContent()
            except:
                raise
    else:
        pass
#-----------------------------------------------------------------------------
curDir_rRNA = widgets.Label(value=os.getcwd())

print('current Directory:')

main = widgets.Select(description='Directories:')
sub_rRNA = widgets.Select(options=['MMM'], value='MMM',
                          description='Files:')
output_Dir_rRNA= widgets.Text(
    placeholder='name of dir for indexing rRNA',
        description='dir name:', disabled=False) 
#----outputdir
box = widgets.HBox([main, sub_rRNA, output_Dir_rRNA])  
#--->output_Dir_rRNA.value

display(curDir_rRNA, box)

test = True

get_DirectoryContent()

main.observe(main_Observer)


# In[ ]:


print('input: ', curDir_rRNA.value,'\n' 'output: ', 
      output_Dir_rRNA.value, '\n' 'file: ', sub_rRNA.value)


# In[ ]:


###widgets for adapter trimming 
pre_prepared_trimming_option = widgets.Select(options=[None, 
        "regular 3' adapter", "regular 5'adapter", 
        "non-internal 3' adapter", "non-internal 5' adapter",
        "anchored 3' adapter", "anchored 5' adapter",
        "5' or 3' (botch possible)", 
        "linked adapter", "non-anchored linked adapter"
        ],
        value=None, rows=8,
        description='Trim-Type:',disabled=False,
        )
       
txt_ind1= widgets.Text(
        value='GTAC', placeholder='Your Sequence',
        description='Input here:', disabled=False
        )

accordion1 = widgets.Accordion(children=[VBox([txt_ind1])])

accordion2 = widgets.Accordion(children=[VBox(
    [pre_prepared_trimming_option])])

accordion1.set_title(0, 'Paste Your Sequence below')
accordion2.set_title(0, 'Choose from options below')

tab_nest = widgets.Tab()
tab_nest.children = [accordion1, accordion2]
tab_nest.set_title(0, 'Unique')
tab_nest.set_title(1, 'definite')
       
display (tab_nest)
    


# In[ ]:


print (pre_prepared_trimming_option.value, txt_ind1.value)


# In[ ]:


#create output dir
def get_fasta_names():
    """
    function nees path
    goes through dir - 
        reads file names with certain extension 
        (first cuts the \n) - loads them in list 
    parameter 
        filenames = basenames 
        list of strings
    return 
        output_dir 
    """
    os.chdir(curDir.value)
    path_list = os.popen('ls').readlines()

    fastq_dir_list, Nr_list=[], [] #creates a new list 
    possible_ext = [".fastq", ".fq.gz",
                    ".fastq.gz", ".fasta", ".fa", 
                    ".fa.gz",".fasta.gz"]
    
    for i in range(len(path_list)):
        if path_list[i].endswith('\n'):
            path_list[i] = path_list[i][:-1]

        for e in possible_ext:
            if path_list[i].endswith(e):
                Nr_list.append('_'+path_list[i])             
                #list for selection widget
                fastq_dir = path_list[i][:-len(e)]+ "_fastqc"
                fastq_dir_list.append(fastq_dir) 
                #list of all files with _fastqc as endings 
    return fastq_dir_list, Nr_list  


def create_each_dir(fastq_dir_list:list):
    """
    an overall dir was created in prior subfunction
    here:
    each fast.gz file will have its own dir for 
    the output//creates dir w/ basenames of files 
    """
    #change directory
    output_path=output_Dir.value
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    #os.mkdir (output_path) #creates a directory 
    os.chdir(output_path)
    cwd_n=os.getcwd()
    os.chdir( cwd_n+'/' ) 
    #os.chdir (os.path.join( output_path,'/'))
    output_subpath = os.getcwd()+'/'
    #read files in current directory (here it should be empty)
    subpath_list = os.popen('ls').readlines()
    new_fast_dir_list=[]
    #iterate over the list of files created in get_fasta_names 
    for g in range(len(fastq_dir_list)):
        new_fast_dir_list.append(output_subpath+'_'+fastq_dir_list[g])
    for dire in new_fast_dir_list:
        #actually make the directory 
        if not os.path.exists(dire):
            os.mkdir(dire)
 
    return output_subpath

fastq_dir_list, Nr_list=get_fasta_names()
output_subpath=create_each_dir(fastq_dir_list)




# In[ ]:


print(output_subpath)


# In[ ]:


def run_all_fastqc(Nr_list:list, output_subpath: str):
    """
    it runs all fastqc quality analysis and stores output at 
    """
    print (Nr_list, '\n', output_subpath)
    for no in range(len(Nr_list)):
        run (curDir.value+'/'+Nr_list[no][1:],
             output_subpath+'_'+fastq_dir_list[no])


# In[ ]:


#run fastgc
def run(fastq_path, output_dir):
    """
    Run the fastqc program on a specified fastq file and return the 
    output directory path.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file to analyze.##############
        ##############-from widget

    options : list of str (optional)
        Command-line flags and options to fastqc. 

    Returns
    -------
    output stored in pre-saved dir 
        fastqc format
    """
    command = "fastqc -o {} {}".format( output_dir, fastq_path) 
    #output then shere its stored
    subprocess.check_call(command, stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE, shell=True) 
    print (command)
    subprocess.check_output(command, shell=True)



# In[ ]:


run_all_fastqc(Nr_list, output_subpath)


# In[ ]:


def run_cutadapt(fastq_path, output_name, kind1, kind2):
    """
    Run cutadapt

    Parameters
    ----------
    fastq_path : str
        path to raw data -> doesnt matter if comressed or not

    ourput_dir : str
        same directory as for fastqc analysis 
        
    kind : str
        individual or pre-prepared
        if kind == ind -> input has to be the sequence 
        if kind == pre-prepared 
            -> input has to be what kind of trimming exactly 

    Returns
    -------
    output stored in pre-saved dir 
        txt format
    """
    command_line_options=[[],] #for pre-prepared
    command_line_options=[["regular 3' adapter", '-a ADAPTER'],
                          ["regular 5'adapter", '-g ADAPTER'], 
                          ["non-internal 3' adapter", '-a ADAPTERX'], 
                          ["non-internal 5' adapter", '-g XADAPTER'], 
                          [ "anchored 3' adapter", '-a ADAPTER$'],
                          ["anchored 5' adapter", '-g ^ADAPTER'], 
                          ["5' or 3' (botch possible)", '-b ADAPTER'], 
                          ["linked adapter", '-a ADAPTER1...ADAPTER2'],
           ["non-anchored linked adapter", '-g ADAPTER1...ADAPTER2']
                         ]
    for elem in range(len(command_line_options)):
        if kind1 == command_line_options[elem][0]: 
            try:
                command = "cutadapt -a {} -o {} {}".format(
                    kind1, 
                    output_name, 
                    fastq_path)
                subprocess.check_call(command, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, shell=True) 
            except subprocess.CalledProcessError as e: 
                print("Error executing cutadapt\npreprep")
                pass 
        
        else:

            try:
                command = "cutadapt -a {} -m 18 -M 40 -o {} {}".format(
                    kind2,
                    output_name, 
                    fastq_path)
                subprocess.check_call(command, stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE, shell=True) 
            except subprocess.CalledProcessError as e: 
                print("Error executing cutadapt\nindividual")
                pass
    print(command)


# In[ ]:


def move_files_cutadapt(destination_path: str):
    """
    literally just moves files
    """
    os.chdir(curDir.value+'/')
    files_in_cwd=os.popen('ls').readlines() #reads dir
    print (files_in_cwd)
    files_in_cwd_wo=[] #new list to append files from dir
    for wh in files_in_cwd:
        files_in_cwd_wo.append(wh.rstrip('\n')) #strip the \n 
    print (files_in_cwd_wo)
    for f in range(len(files_in_cwd_wo)):
        
        if files_in_cwd_wo[f].startswith('trimmed_'):
            shutil.move(files_in_cwd_wo[f], destination_path)
            print ('{} was moved into following directory {}'.format (
                files_in_cwd_wo[f], destination_path)
                  )
        else: 
            print (files_in_cwd_wo[f]+' couldnt be moved or found')
            pass


# In[ ]:


def run_all_cutadapt(input_dir, output, 
                     pre_prepared_trimming_option, 
                     txt_ind1):
    """
    it runs all cutadapt quality analysis and stores output at 
    def run_cutadapt(fastq_path, output_name, kind1, kind2)
    fastq_path----       
    #curDir.value wehere all the raw data is stored 
    #Nr_list[no][1:] Nr_1.fq.gz without the '_' infront
    output---- only the name
    kind1--- trimming option -> pre-prepared
    kind2--- individual adapter sequence
    """
    for no in range(len(Nr_list)): 
        run_cutadapt (curDir.value+'/'+Nr_list[no][1:], 
                      'trimmed_'+Nr_list[no][1:], 
                      pre_prepared_trimming_option.value, txt_ind1.value)
    
    

run_all_cutadapt(curDir.value, 
                 Nr_list, 
                 pre_prepared_trimming_option, 
                 txt_ind1)


# In[ ]:


def run_all_move_files_cutadapt(Nr_list:list,
                                output_subpath: str, 
                                fastq_dir_list:list):
    """
    it runs all move_files_cutadapt stores output at output_subpath
    """
    for no in range(len(Nr_list)):
        move_files_cutadapt (output_subpath+'/')


# In[ ]:


run_all_move_files_cutadapt(Nr_list, output_subpath, fastq_dir_list)


# In[ ]:


main.observe(main_Observer)


# In[ ]:


print('input: ', 
      curDir_rRNA.value,'\n' 'output: ',
      output_Dir_rRNA.value, '\n' 'file: ', 
      sub_rRNA.value)


# In[ ]:


def get_trimmed_names(output_subpath):
    """
    gets names of trimmed outputfiles
    """
    trimmed_List=[] #with the names of trimmed files
    os.chdir(output_subpath+'/')
    path_list_rRNA_filter = os.popen('ls').readlines()

        
    for i in range(len(path_list_rRNA_filter)):
        if path_list_rRNA_filter[i].startswith('trimmed_'):
            trimmed_List.append(path_list_rRNA_filter[i][:-1]) 
    return trimmed_List

def run_rRNA_index(rRNA_full_sequences_path, 
                   output_Dir_rRNA, sub_rRNA, output_subpath):
    """
    Build an index for rRNA.

    Parameters
    ----------
    rRNA_full_sequences_path : str
        Path to planarian rRNA full sequences.
        ############################-from widget

    output_dir_rRNA : str
        command line to outputpath. 
        
    sub_rRNA.value : str
        file chosen from widget --> full planarian rRNA sequences
        
    Returns
    -------
    index for rRNA 
    """
    if not os.path.exists(output_subpath+'/'+output_Dir_rRNA.value):
        os.makedirs(output_subpath+'/'+output_Dir_rRNA.value)
    
    try:
        command = "indexdb_rna --ref {},{} -v".format(
            curDir_rRNA.value+'/'+sub_rRNA.value,
            output_subpath+output_Dir_rRNA.value+'/'+sub_rRNA.value) 
            #input path, output Dir path
        subprocess.check_call(command, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, shell=True) 
        subprocess.check_output(command, shell=True)
        print (command)
    except:
        print ('Error')

#inititate the function
trimmed_List=get_trimmed_names(output_subpath)
run_rRNA_index(curDir_rRNA.value, 
               output_Dir_rRNA, 
               sub_rRNA, 
               output_subpath)


# In[ ]:


def unzip_all_gz(output_subpath: str):
    """
    unzips all fq.gz files in output_subpath...> 
    /data1/AnLi_FM_MA/data/output_dir_xxyyxx/
    """
    os.chdir(output_subpath)
    os.getcwd()
    os.popen('ls').readlines()
    unzip_me='gunzip *fq.gz'
    subprocess.check_call(unzip_me, 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE, 
                          shell=True)   
unzip_all_gz(output_subpath)


# In[ ]:


def run_rRNA_filter(rRNA_full_sequences_path, 
                    output_Dir_rRNA, 
                    sub_rRNA, 
                    trimmed_Location, 
                    name_aligned, 
                    name_other_rejected_reads,
                    output_subpath):
    """
    Build an index for rRNA.

    Parameters
    ----------
    rRNA_full_sequences_path : str
        Path to planarian rRNA full sequences.
        ############################-from widget

    output_dir_rRNA : list of str (optional)
        Command-line to output path. 
        
    sub_rRNA.value : str
        file chosen from widget --> 
        full planarian rRNA sequences
        
    Returns
    -------
    index for rRNA 
    """   
    os.chdir(output_subpath+'/')
    print (os.getcwd())

    
    c="sortmerna --ref {},{}--reads'{}'--fastx "+
    "-a 7--aligned{}--other{}--log -v".format(
        curDir_rRNA.value+'/'+sub_rRNA.value,
                output_subpath+output_Dir_rRNA.value+'/'+sub_rRNA.value, 
                trimmed_Location,
                name_aligned,
                name_other_rejected_reads
                ) 
    try:
        subprocess.check_call(c, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True) 
    except:
        print('something is wrong')
        pass
    subprocess.check_output(c, 
                            shell=True)
    


# In[ ]:


def run_all_rRNA_filter(output_subpath, 
                        trimmed_List, 
                        output_Dir_rRNA, 
                        curDir_rRNA, 
                        sub_rRNA):
    """
    run above function (filter rRNA)
    """
    #destination_path --> output_subpath+'/'
    for no in range(len(trimmed_List)):
        run_rRNA_filter(curDir_rRNA.value+sub_rRNA.value,
                    output_Dir_rRNA, sub_rRNA, 
                    output_subpath+trimmed_List[no][:-3], 
                    'rRNA_hits_aligned_reads_'+trimmed_List[no][:-3],
                    'rejected_reads_'+trimmed_List[no][:-3],
                    output_subpath
                      )
        


# In[ ]:


run_all_rRNA_filter(output_subpath, trimmed_List, 
                    output_Dir_rRNA, curDir_rRNA, sub_rRNA)


# In[ ]:


#run fastgc
def run_fastqc_5(fastq_report_path, output_subpath):
    """
    Run the fastqc program on a specified fastq file and return the 
    output directory path.

    Parameters
    ----------
    fastq_path : str
        Path to the fastq file to analyze.
        ############################-from widget

    options : list of str (optional)
        Command-line flags and options to fastqc. 

    Returns
    -------
    output stored in pre-saved dir 
        fastqc format
        --> fastq_report_path
    """
    command = "fastqc -o {} {}".format(fastq_report_path, 
                                       output_subpath) 
    subprocess.check_call(command, stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE, shell=True) 
    print (command)
    subprocess.check_output(command, shell=True)


# In[ ]:


def run_all_fastqc_5(output_subpath):
    """
    initiates fastqc analysis of *fq in all dirs in output dir
    """
    os.chdir(output_subpath)

    if not os.path.exists(output_subpath+'/'+'fastq_report_fq_all'):
        os.mkdir(output_subpath+'/'+'fastq_report_fq_all')
        fastq_report_path=os.getcwd()
    fastq_report_path=output_subpath+'fastq_report_fq_all'
    os.chdir(output_subpath)
    files_in_output_dir = os.popen('ls').readlines()

    
    for file in files_in_output_dir:
        file=file.split('\n')

        print(file)
        try:
            if file[0].endswith('fq'):
                run_fastqc_5(fastq_report_path, output_subpath+file[0])
        except:
            print("something's wrong")
            
    return fastq_report_path 
    


# In[ ]:


fastq_report_path=run_all_fastqc_5(output_subpath)


# In[ ]:


path_to_piPipes=input('Path to piPipes :') #//piPipes-master/bin
input_file_path=input('Input File Path :') #//output_dir_sm2_180918 or 
#//data/output_sm3_KD
output_dir_name=input('Output Dir Name :') #insert_file_piPIpes


def get_fq_files_from_output_dir(path_to_piPIpes, input_file_path,
                                 output_dir_name):
    """
    get files from dir
    run piPipes_fastq_to_insert on shell by iterating over fq files
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
                input_file_path+'/'+output_dir_name+'/'+
                re_re_namesss[r]+'.insert'
                                                            )
            subprocess.check_call(command, stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE, shell=True)
            print(command)
    except:
        command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq',
                input_file_path+'/'+output_dir_name+'/'+
                re_re_namesss[r]+'.insert'
                                          )
        print('{} command couldnt be run'.format(command))
    try:
        for r in range(len(re_re_namesss)):
            command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq.fq',
                input_file_path+'/'+output_dir_name+'/'+
                re_re_namesss[r]+'.insert'
                                          )
            subprocess.check_call(command, 
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE, 
                                  shell=True)
            print(command)
    except:
        command="./piPipes_fastq_to_insert {} {}".format(
                input_file_path+'/'+re_re_namesss[r]+'.fq.fq',
                input_file_path+'/'+output_dir_name+'/'+
                re_re_namesss[r]+'.insert'
                                         )
        print('{} command couldnt be run'.format(command))

get_fq_files_from_output_dir(path_to_piPipes, 
                             input_file_path, 
                             output_dir_name)


# In[ ]:


path_to_genome=input('Path to genome/genome_name : ') 
#path/final_dd_Smed_g4.fa
path_to_bowtie=input('Path to Bowtie-build : ')
#path/piPipes-master/bin
input_file_path=input('Input File Path : ') 
#path/insert_file_piPIpes
output_dir_path=input('Output File Path : ') 
#path/bowtie_index
name_of_output_index_genome=input('Name of indexed genome : ')
#Smed_g4 



def get_index_genome_create_sam_bam_files_resp(
    path_to_genome, path_to_bowtie, 
    input_file_path, output_dir_path, 
    name_of_output_index_genome):
    """
    generate index for planarian genome
    create bam files--->mapping to genome
    """
    insert=[]
    #----------------------------------------------------------
    if not os.path.exists(output_dir_path): #make dir for output
        os.mkdir(output_dir_path)
    
    #generate index for planarian genome-----------------------
    os.chdir(path_to_bowtie)
    os.getcwd()
    command="bowtie-build {} {}".format(
                path_to_genome, 
        output_dir_path+'/'+name_of_output_index_genome
                                        )
    if name_of_output_index_genome not in output_dir_path:
        subprocess.check_call(command, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(command)
    else:
        print('Index for genome was already generated')
    #----------------------------------------------------------
    os.chdir(input_file_path)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #------------------------------------
    for i in range(len(list_files)): 
        #goes through the files in the listed files
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.insert') and file not in insert:
            insert.append(file)   
    for file in insert:
        bowtie_sam='(bowtie -r -p7 -v0 -a --best --strata'+
        ' {} {} -S {}) 2>> {}'.format(
            output_dir_path+'/'+name_of_output_index_genome,
            input_file_path+'/'+file, 
            output_dir_path+'/'+file+'.sam',
            output_dir_path+'/'+'bowtie_'+file+'.log'
                                                    )
        
        samtools_bam='samtools view -bS -@ 7 {} |'+
        ' samtools sort -@ 7 -> {}'.format(
        output_dir_path+'/'+file+'.sam',
        output_dir_path+'/'+file+'.bam')
        remove_sam='rm {}'.format(output_dir_path+'/'+file+'.sam')
                  
        subprocess.check_call(bowtie_sam, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, shell=True)
            
        subprocess.check_call(samtools_bam, stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, shell=True)

        subprocess.check_call(remove_sam, stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, shell=True)
          
get_index_genome_create_sam_bam_files_resp(path_to_genome, 
                            path_to_bowtie,
                            input_file_path,
                            output_dir_path,
                            name_of_output_index_genome)


# In[ ]:


path_to_piPipes=input('Path to bedtools_piPipes : ')
#//piPipes-master/bin
path_to_input_bam=input('Path to bam files (input) : ') 
#//planarian_genome/final_dd_Smed_g4.fa
path_to_output_bed=input('Path to output w/o name of dir : ') 
#//planarian_genome/final_dd_Smed_g4.fa
output_dir_name_bed=input('Output File Name (bed files) : ') 
#//bowtie_index
import re

def get_bam_files_run_bedtools_piPipes(path_to_piPipes,
                                        path_to_input_bam, 
                                        path_to_output_bed, 
                                        output_dir_name_bed):
    """
    generate bedfiles out of bam files with piPipes
    """
    bed_files=[]
    #----------------------------------------------------------
    if not os.path.exists(path_to_output_bed+'/'+
                          output_dir_name_bed): #make dir for output
        os.mkdir(path_to_output_bed+'/'+
                 output_dir_name_bed)
    #----------------------------------------------------------
    os.chdir(path_to_input_bam)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #----------------------------------------------------------

    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.bam') and file not in bed_files:
            bed_files.append(file)
    print((bed_files))

    os.chdir(path_to_piPipes)
    ozzy=os.getcwd()
    print(ozzy)
    

    #----------------------------------------------------------
    for file in bed_files:

        motif_name=re.findall(r'(\w{1}\d{1}-\w{1,2}\d{1,2})', file)
        #_x1-sm31 32 33
        motif_name=re.findall(r'(\w{1}\d{1}-\w{3}\d{1})', file) 
        #x1-unc1 2 3
        motif_name=re.findall(r'(\w{4}-\w{1,2}\d{1,2})', file) 
        #xins-sm31 32 33
        motif_name=re.findall(r'(\w{4}-\w{3}\d{1})', file) 
        #xins-unc1 2 3
        motif=motif_name[0]
        command="./bedtools_piPipes bamtobed -i {} > {}.bed".format(
                    path_to_input_bam+file, 
            path_to_output_bed+output_dir_name_bed+'/'+motif
                                            )

        subprocess.check_call(command, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE,
                              shell=True)
        print(command)            
            
            
            
get_bam_files_run_bedtools_piPipes(path_to_piPipes,
                                   path_to_input_bam, 
                                   path_to_output_bed,
                                   output_dir_name_bed)


# In[ ]:


path_to_piPipes=input('Path to bedtools_piPipes : ')
#//piPipes-master/bin
sortmeRNA_insert_path=input('Path to sortmeRNA.insert files : ')
#//insert_file_piPIpes/
path_to_input_bed=input('Path to input bed files : ') 
#//data/...
path_to_output_bed2=input('Output path bed2 : ') 
#/data/...

def get_bed_files_runbed2(path_to_piPipes, sortmeRNA_insert_path,
                                path_to_input_bed, 
                                               path_to_output_bed2):
    """
    generate bed2files out of bed files with piPipes
    """
    bed_files, insert_files=[], []
    #----------------------------------------------------------
    os.chdir(path_to_input_bed)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #----------------------------------------------------------

    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.bed') and file not in bed_files:
            bed_files.append(file)
    print((bed_files))
    #----------------------------------------------------------
    os.chdir(sortmeRNA_insert_path)
    os.getcwd()
    list_insert_files=os.popen('ls').readlines()
    #----------------------------------------------------------

    for i in range(len(list_insert_files)): 
        files=list_insert_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.insert') and file not in insert_files:
            insert_files.append(file)
    print((insert_files))
    
    #change dir to piPipes-------------------------------------
    os.chdir(path_to_piPipes)
    ozzy=os.getcwd()
    print(ozzy)
    #----------------------------------------------------------
    for i in range(len(bed_files)):
        for fi in range(len(insert_files)):
            if i==fi:
                cd="./piPipes_insertBed_to_bed2 {} {} > {}2".format(
                    sortmeRNA_insert_path+insert_files[i],
                            path_to_input_bed+bed_files[i], 
                    path_to_output_bed2+bed_files[i]
                                                    )

                subprocess.check_call(cd, 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE, 
                                      shell=True)  
            else:
                pass

            
            
get_bed_files_runbed2(path_to_piPipes, 
                      sortmeRNA_insert_path,
                      path_to_input_bed, 
                      path_to_output_bed2)


# In[ ]:


path_to_input_bed2=input('Path to input bed2 files : ') 
#//planarian_genome/final_dd_Smed_g4.fa
output_path_length=input('Path to output w/ new dir name : ')

import re

def get_bed2_print_length(path_to_input_bed2, output_path_length):
    """
    input
        bed2 files -> contain sequences 
        outputpath_length: dir where output is stored
        
    get sequences form seventh column and only leave those, 
    which are unique, then print the lengths of 
    first column of these, sort and get only the unique ones 
    return
    txt file with the length of sequences and theit occurence
    """
    bed2_files=[]
    #----------------------------------------------------------
    if not os.path.exists(output_path_length): #make dir for output
        os.mkdir(output_path_length)
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #------------------------------------
    for i in range(len(list_files)):
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.bed2') and file not in bed2_files:
            bed2_files.append(file)   
    for file in bed2_files:
        motif_name=re.findall(r'(\w{1}\d{1}-\w{1,2}\d{1,2})', file)
        #_x1-sm31 32 33
        motif_name=re.findall(r'(\w{1}\d{1}-\w{3}\d{1})', file) 
        #x1-unc1 2 3
        motif_name=re.findall(r'(\w{4}-\w{1,2}\d{1,2})', file) 
        #xins-sm31 32 33
        motif_name=re.findall(r'(\w{4}-\w{3}\d{1})', file)
        #xins-unc1 2 3
        motif=motif_name[0]
        if not os.path.exists(output_path_length+'/output_length_'+motif):
            command="cat '{}' | cut -f7 | sort | uniq |"+
            " awk '{}' - | sort | uniq -c > output_length_{}".format(
                path_to_input_bed2+'/'+file,
                '{print length($1)}',
                motif
                     )

            subprocess.check_call(command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, 
                                  shell=True)
            print(command)   
        else:
            
            command="cat '{}' | cut -f7 | sort | "+
            "uniq | awk '{}' - | sort | uniq -c >> output_length_{}".format(
                path_to_input_bed2+'/'+file,
                '{print length($1)}',
                motif
                     )

            subprocess.check_call(command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, 
                                  shell=True)
            print(command)

            
get_bed2_print_length(path_to_input_bed2, output_path_length)


def move_files(path_to_input_bed2, output_path_length):
    """
    literally moves files
    input: 
        path to the dir, where files are 
        stored and the path, where output is supposed to be stored 
    creates new dir and stores files in new dir
    return:
        new dir with output files stored there
    """
    destination_path=output_path_length
    os.chdir(path_to_input_bed2)
    files_in_cwd=os.popen('ls').readlines() #reads dir
    print (files_in_cwd)
    files_in_cwd_wo=[] #new list to append files from dir
    for wh in files_in_cwd:
        files_in_cwd_wo.append(wh.rstrip('\n')) #strip the \n 
    print (files_in_cwd_wo)
    for f in range(len(files_in_cwd_wo)):
        
        if files_in_cwd_wo[f].startswith('output_length'):
            shutil.move(files_in_cwd_wo[f], destination_path)
            print ('{} was moved into following directory {}'.format (
                files_in_cwd_wo[f],
                destination_path))
        else: 
            print (files_in_cwd_wo[f]+' couldnt be moved or found')
            pass
move_files(path_to_input_bed2, output_path_length)


# In[ ]:


path_to_input_bed2=input('Path to input bed2 files : ') 
#/bowtie_index/bed_files/
output_path_cutoff25=input('Path to output w/ new dir name : ')
#//bowtie_index/bed_files/greater_than_25/

def get_bed2_cutoff_25(path_to_input_bed2, output_path_cutoff25):
    """
    stores bed2 files with only the lengths over 25 to
    reassure only piRNAs are in bed2 files 
        -rest are miRNA, Poly A Tails or other sequences
    """
    bed2_files=[]
    #----------------------------------------------------------
    if not os.path.exists(output_path_cutoff25): #make dir for output
        os.mkdir(output_path_cutoff25)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #------------------------------------
    for i in range(len(list_files)):
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.bed2') and file not in bed2_files:
            bed2_files.append(file)
  
    for file in bed2_files:

            command="awk 'length($7)>25' {} > {}/25_cutoff_{}".format(
                path_to_input_bed2+file,
                output_path_cutoff25,
                file
                    )

            subprocess.check_call(command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, 
                                  shell=True)
            print(command)   
        
            
get_bed2_cutoff_25(path_to_input_bed2, output_path_cutoff25)

def create_fasta_files(output_path_cutoff25):
    """
    creates fasta files out of bed2 files
    """
    cutoff_bed2=[]
    #----------------------------------------------------------
    if not os.path.exists(output_path_cutoff25+'/fasta_cutoff'): 
        os.mkdir(output_path_cutoff25+'/fasta_cutoff')
    #----------------------------------------------------------
    os.chdir(output_path_cutoff25)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #------------------------------------ 
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') #looses the \n= newline
        file=file[0]
        if file.endswith('.bed2') and file not in cutoff_bed2:
            cutoff_bed2.append(file)
    #----------------------------------------------------------
    for file in cutoff_bed2:
            op_file=file.split('.bed2')
            op_file=op_file[0]
            com="cat {} | cut -f7 |sort|uniq| awk '{}'->{}.fa".format(
                output_path_cutoff25+'/'+file,
                '{print ">"$1"\\n"$1}',
                output_path_cutoff25+'/fasta_cutoff/'+op_file
                                            )
            print(com) 

            subprocess.check_call(com, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, shell=True)
    
create_fasta_files(output_path_cutoff25) 


# In[ ]:


path_to_input_fasta=input('Path to input bed2 files : ') 
#//fasta_cutoff/
output_path_RNA=input('Path to output w/ new dir name : ')
#//T_U_exchg

def exchange_U_w_T(path_to_input_fasta, output_path_RNA):
    """
    exchanges Uridin with Thymin --> DNA to RNA
    """
    fasta_files=[]
    #----------------------------------------------------------
    if not os.path.exists(output_path_RNA): 
        os.mkdir(output_path_RNA)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_fasta)
    os.getcwd()
    list_files=os.popen('ls').readlines()
    #------------------------------------
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.endswith('.fa') and file not in fasta_files:
            fasta_files.append(file)
    
    for file in fasta_files:
            op_file=file.split('25_cutoff_') 
            op_file=op_file[1]
            #print(op_file)            
            command="sed 's/T/U/g' {} > {}/U_RNA_{}".format(
                path_to_input_fasta+file, 
                output_path_RNA, 
                op_file
                        )
            subprocess.check_call(command, 
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, 
                                  shell=True)
            print(command)   
        
exchange_U_w_T(path_to_input_fasta, output_path_RNA)


# In[ ]:


path_to_input_txt=input('Path to input txt files (unique KD/WT) : ') 
#//T_U_exchg/'
path_to_input_bed2=input('Path to input bed2 files : ') 
#//bed_files/'
output_path_mapped_bed2=input('Path to output w/ new dir name : ')
#path_to_input_txt+ /mapped_unique

def grep(path_to_input_txt, path_to_input_bed2, 
         output_path_mapped_bed2):
    """
    input parameter
        txt file: only unique sequences 
        from either the KD or WT, which dont overlap in both (from vendiagram)
        bed2 file: mapped sequences to genome
        output path: dir for the output
    return
        bed2 file, includes only the sequences, 
        that are included in either the KD or WT
    """
    txt,bed2=[], []
    
    #----------------------------------------------------------
    if not os.path.exists(output_path_mapped_bed2):
        os.mkdir(output_path_mapped_bed2)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_txt)
    os.getcwd()
    list_files=os.popen('ls').readlines()       
    #------------------------------------
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.startswith('sm3_unique_from_'):
            file=file.split('sm3_unique_from_')
            file=file[1]
        if file.endswith('.txt') and file not in txt:
            txt.append(file)   
    #print(txt)
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    list_files=os.popen('ls').readlines()   
    #------------------------------------
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.startswith('d_'):
            file=file.split('d_')
            file=file[1]
        if file.endswith('.bed2') and file not in bed2:
            bed2.append(file)
    dictup={}
   #----------------------------------------------------------
    for bedfile in bed2:
        for txtfile in txt:
            if bedfile.startswith('xins-') and txtfile.endswith('_Xins.txt'):
                if bedfile.endswith('sm31.bed2') and txtfile.startswith('sm31'):
                        if bedfile not in dictup:
                                dictup[bedfile]=txtfile
                elif bedfile.endswith('sm32.bed2') and txtfile.startswith('sm32'):
                        if bedfile not in dictup:
                                dictup[bedfile]=txtfile
                elif bedfile.endswith('sm33.bed2') and txtfile.startswith('sm33'):
                        if bedfile not in dictup:
                                dictup[bedfile]=txtfile                        
            elif bedfile.startswith('x1') and txtfile.endswith('X1.txt'):
                if bedfile.endswith('sm31.bed2') and txtfile.startswith('sm31'):
                        if bedfile not in dictup:
                                dictup['d_'+bedfile]=txtfile                        
                elif bedfile.endswith('sm32.bed2') and txtfile.startswith('sm32'):
                        if bedfile not in dictup:
                                dictup['d_'+bedfile]=txtfile                        
                elif bedfile.endswith('sm33.bed2') and txtfile.startswith('sm33'):
                        if bedfile not in dictup:
                                dictup['d_'+bedfile]=txtfile                        

    for key, val in dictup.items():
        
        op_file=val.split('.txt')[0] 


        command = "grep -Fwf {} {} > {}/mapped_{}.bed2".format(
            path_to_input_txt+'/'+'sm3_unique_from_'+val, 
            path_to_input_bed2+'/'+key,
            output_path_mapped_bed2, op_file)
        sort="sort -k1,1 -k2,2n {}/mapped_{}.bed2 > "+
        "{}/sorted_mapped_{}.bed2".format(
            output_path_mapped_bed2,
            op_file,output_path_mapped_bed2, 
            op_file) 

                    
                                                                  
        subprocess.check_call(command, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(command)
        subprocess.check_call(sort, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(sort)
    
grep(path_to_input_txt, path_to_input_bed2, output_path_mapped_bed2)  





def grep_unc(path_to_input_txt, path_to_input_bed2, output_path_mapped_bed2):
    """
    input parameter
        txt file: only unique sequences from 
        either the KD or WT, which dont overlap in both (from vendiagram)
        bed2 file: mapped sequences to genome
        output path: dir for the output
    return
        bed2 file, includes only the sequences, 
        that are included in either the KD or WT
    """
    txt,bed2=[], []
    #----------------------------------------------------------
    if not os.path.exists(output_path_mapped_bed2): 
        os.mkdir(output_path_mapped_bed2)
    #----------------------------------------------------------
    os.chdir(path_to_input_txt)
    os.getcwd()
    list_files=os.popen('ls').readlines()      
    #------------------------------------
    for i in range(len(list_files)): 
        if list_files[i].startswith('unc_'):
            files=list_files[i]
            file=files.split('\n')
            file=file[0]
            if file.startswith('unc_unique_from_'):
                file=file.split('unc_unique_from_')
                file=file[1]
            if file.endswith('.txt') and file not in txt:
                txt.append(file)   
        else:
            pass
    #print('txt :', txt)
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    list_files=os.popen('ls').readlines()  
    #------------------------------------
    for i in range(len(list_files)): 
        if list_files[i].startswith('d_x1-unc') or 
        list_files[i].startswith('xins-unc'):
            files=list_files[i]
            file=files.split('\n') 
            file=file[0]
            if file.startswith('d_'):
                file=file.split('d_')
                file=file[1]
            if file.endswith('.bed2') and file not in bed2:
                bed2.append(file)
    
    #print('bed :', bed2)
    dictup={}
   #----------------------------------------------------------
    for bedfile in bed2:
        for txtfile in txt:
            if bedfile.startswith('xins-') and txtfile.endswith('_Xins.txt'):
                if bedfile.endswith('1.bed2') and txtfile.startswith('sm31'):
                        if bedfile not in dictup:
                                dictup[bedfile]=txtfile
                elif bedfile.endswith('2.bed2') and txtfile.startswith('sm32'):
                        if bedfile not in dictup:
                                dictup[bedfile]=txtfile
                elif bedfile.endswith('3.bed2') and txtfile.startswith('sm33'):
                        if bedfile not in dictup:
                                dictup[bedfile]=txtfile                        
            elif bedfile.startswith('x1') and txtfile.endswith('X1.txt'):
                if bedfile.endswith('1.bed2') and txtfile.startswith('sm31'):
                        if bedfile not in dictup:
                                dictup['d_'+bedfile]=txtfile                        
                elif bedfile.endswith('2.bed2') and txtfile.startswith('sm32'):
                        if bedfile not in dictup:
                                dictup['d_'+bedfile]=txtfile                        
                elif bedfile.endswith('3.bed2') and txtfile.startswith('sm33'):
                        if bedfile not in dictup:
                                dictup['d_'+bedfile]=txtfile                        

    print('dictionary :', dictup)
    for key, val in dictup.items():
        
        op_file='unc_'+val.split('.txt')[0] 
        #op_file=op_file[0]
        print(op_file)
        command = "grep -Fwf {} {} > {}/mapped_{}.bed2".format(
            path_to_input_txt+'/'+'unc_unique_from_'+val, 
            path_to_input_bed2+'/'+key,
            output_path_mapped_bed2, 
            op_file)
        sort="sort -k1,1 -k2,2n {}/mapped_{}.bed2 > {}/"+
        "sorted_mapped_{}.bed2".format(
            output_path_mapped_bed2,
            op_file,
            output_path_mapped_bed2, 
            op_file) 

                    
                                                                  
        subprocess.check_call(command, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(command)
        subprocess.check_call(sort, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(sort)
        
grep_unc(path_to_input_txt, path_to_input_bed2, output_path_mapped_bed2)  


# In[ ]:


import sys,os
import csv
import dask.dataframe as dd
import pandas as pd
import numpy as np
from glob import glob
import subprocess

path='/vendiagram/mapped_unique/intersected/'
os.chdir(path)
path=os.getcwd()

list_data=os.popen('ls').readlines() #list of filenames

intersect=[]
for i in list_data:
    if i.startswith('intersected_sorted_mapped'):
        file=i.split('\n')
        file=file[0]
        intersect.append(file)
        
files=[]
for i in intersect:
    try:
        df=dd.read_csv(path+'/'+i, sep= '\t', engine= 'python',
                       header=None, dtype={11:'object'})
        df.header=i #name 
        files.append(df)
    except pd.io.common.EmptyDataError:
        print (i, " is empty")
        
for col in files:
    col['17']=col[3]/col[4]
    col.to_csv(path+'/'+'weighted_counts_{}.*.bed'.format(
        col.header.split('.bed')[0]), 
               index=False, header=False, sep= '\t')

#write all chunks into one dataframe
for file in files:
    with open('all_{}'.format(file.header), 'a+') as out:
        for fn in filenames:
                #print(fn[42:-8])
                #print(file.header[26:-4])
                if fn[42:-8] == file.header[26:-4]:
                    print('file.header', file.header[26:-4])
                    with open(fn) as f:
                        out.write(f.read())   
#remove all chunks
for fi in filenames:
    
    command_remove='rm {}'.format(path+'/'+fi)
    subprocess.check_call(command_remove, stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, shell=True)
    print(command_remove)
    


# In[ ]:


path_to_input_all_annot=input('Path to input TE file (bed file) : ') 
#/data1/AnLi_FM_MA/data/        ---> combined_all_annot_te_sorted.bed
path_to_input_bed2=input('Path to input bed2 files (sorted and mapped) : ') 
#//mapped_unique
output_path_inters_bed=input('Path to output w/ new dir name : ')
#//intersected

import subprocess, os, sys
import optparse, shutil

def inters(path_to_input_all_annot, path_to_input_bed2, 
           output_path_inters_bed):
    """
    input parameter
        bed file: 
            here--> all annotated TE files were combined in one -> 
            and sorted! 
        bed2 file: 
            mapped (and sorted) sequences to genome and  unique 
            to KD/control 
        output path: dir for the output
    return
        bed file:
            has info about where the sequences 
            originated from(-> mapped to all annotated TEs)
    """
    bed2=[]
    
    #----------------------------------------------------------
    if not os.path.exists(output_path_inters_bed):
        os.mkdir(output_path_inters_bed)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_all_annot)
    os.getcwd()
    list_files=os.popen('ls').readlines()      
    #------------------------------------
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.startswith('combined_all_annot_te_sorted'):
            all_TE_file=file     #all annot TE file
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    list_files=os.popen('ls').readlines()      
    #------------------------------------
    for i in range(len(list_files)):
        files=list_files[i]
        file=files.split('\n')
        file=file[0]
        if file.endswith('.bed2') and file.startswith('sorted') 
        and file not in bed2:
            bed2.append(file)
    print(bed2)
    #----------------------------------------------------------
    for bedfile in bed2:
        op_bed=bedfile.split('.bed2')
        op_bed=op_bed[0]
        command = "bedtools intersect -wa -wb"+
        "-split -sorted -a {} -b {} > {}/intersected_{}.bed".format(
                    path_to_input_bed2+'/'+bedfile, 
                    path_to_input_all_annot+all_TE_file,
                    output_path_inters_bed, 
                    op_bed
                          )


        subprocess.check_call(command, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(command)   
inters(path_to_input_all_annot, path_to_input_bed2, output_path_inters_bed)


# In[ ]:


path_to_input_bed=input('Path to input bed files (weighed counts) : ') 
#//DNA_mapped_unique/intersected/
output_path_sense_bed=input('Path to output w/ new dir name (sense) : ')
#//DNA_mapped_unique/intersected/sense_anti

def inters(path_to_input_all_annot, path_to_input_bed2, 
           output_path_inters_bed):
    """
    input parameter
        bed file: 
            intersected (mapped to TE) 
            files with weighed counts as column #17  
        output path: dir for the output
            sense
            antisense
    return
        bed file:
            seperated by if they are sense or antisense
            plus-> sum of column #17 to count percentage later on
    """
    bed=[]
    
    #----------------------------------------------------------
    if not os.path.exists(output_path_sense_bed): 
        os.mkdir(output_path_sense_bed)  
    #----------------------------------------------------------
    os.chdir(path_to_input_bed)
    os.getcwd()
    list_files=os.popen('ls').readlines()      
    #------------------------------------
    for i in range(len(list_files)):
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.endswith('.bed') and file.startswith('weighed') and
        file not in bed:
            bed.append(file)
    #print(bed)
    #----------------------------------------------------------
    for bedfile in bed:
        op_bed=bedfile.split('.bed')
        op_bed=op_bed[0]
        #df['one'] == df['two']
        command_sense = "sense_{} =`cat {} |"+
        " mawk '$6==$13 {print $0}' | mawk '{if ($14=="maker") "+
        "{print $0}}' | mawk -F "\\t" '{ sum+=$18} END {print sum}'`".format(
                        bedfile,   
                        path_to_input_bed+'/'+bedfile
                                                       )
        command_antisense = "antisense_{} =`cat {} |"+
        " mawk '$6!=$13 {print $0}' | mawk '{if ($14=="maker")"+
        " {print $0}}' | mawk -F "\\t" '{ sum+=$8} END {print sum}'`".format(
                        bedfile,   
                        path_to_input_bed+'/'+bedfile


inters(path_to_input_all_annot, path_to_input_bed2, 
       output_path_inters_bed)


# In[ ]:





# In[ ]:


path_to_input_all_annot=input('Path to input TE file (bed file) : ') 
#/data1/AnLi_FM_MA/data/        ---> combined_all_annot_te_sorted.bed
path_to_input_bed2=input('Path to input bed2 files (sorted and mapped) : ') 
#/vendiagram/mapped_unique/
output_path_inters_bed=input('Path to output w/ new dir name : ')
#/fasta_cutoff/vendiagram/mapped_unique/intersected/non_overlap

def inters(path_to_input_all_annot, path_to_input_bed2, output_path_inters_bed):
    """
    
    non overlapping sequences are intersected to the TEs 
    
    input parameter
        bed file: 
            here--> all annotated TE files were combined in one -> and sorted! 
        bed2 file: 
            mapped (and sorted) sequences to genome and  unique to KD/control 
        output path: dir for the output
    return
        bed file:
            has info about where the sequences originated from
            (-> mapped to all annotated TEs)
    """
    bed2=[]
    
    #----------------------------------------------------------
    if not os.path.exists(output_path_inters_bed): 
        os.mkdir(output_path_inters_bed)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_all_annot)
    os.getcwd()
    list_files=os.popen('ls').readlines()      
    #------------------------------------
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.startswith('combined_all_annot_te_sorted'):
            all_TE_file=file                            
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    list_files=os.popen('ls').readlines()   
    #------------------------------------
    for i in range(len(list_files)): 
        files=list_files[i]
        file=files.split('\n') 
        file=file[0]
        if file.endswith('.bed2') and file.startswith('sorted') and 
        file not in bed2:
            bed2.append(file)
    print(bed2)
    #----------------------------------------------------------
    for bedfile in bed2:
        op_bed=bedfile.split('.bed2')
        op_bed=op_bed[0]
        command = "bedtools intersect -wa -sorted"+
        " -a {} -b {} -v > {}/non_overlap_intersected_{}.bed".format(
                             path_to_input_bed2+'/'+bedfile, 
                             path_to_input_all_annot+all_TE_file,
                             output_path_inters_bed, 
                             op_bed
                                    )


        subprocess.check_call(command, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE,
                              shell=True)
        print(command)

inters(path_to_input_all_annot, path_to_input_bed2, output_path_inters_bed)


# In[ ]:


#calculate weighted count
import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob


path='/vendiagram/mapped_unique/intersected/non_overlap'
os.chdir(path)
path=os.getcwd()

list_data=os.popen('ls').readlines() #list of filenames
print(list_data)

intersect=[]
for i in list_data:
    if i.startswith('non_over'):
        file=i.split('\n')
        file=file[0]
        intersect.append(file)
        
#read in all files with pandas        
files=[]
for i in intersect:
    try:
        df=pd.read_csv(path+'/'+i, sep= '\t',
                       engine= 'python', header=None)
        df.header=i
        files.append(df)
    except pd.io.common.EmptyDataError:
        print (i, " is empty")
        
#create new columns with weighted count 
for col in files:
    col[7]=col[3]/col[4]
    col[7] = pd.Series([round(val, 8) for val in col[7]], index = col.index)
    col.to_csv('weighted_counts_{}'.format(col.header), 
               index=False, header=False, sep= '\t')
    
#non overlapping -> from venn diagram 
sm2_d={}
with open('sm2_non_overl_summary.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(['summary'])
for file in files_rounded:
    if file.header not in sm2_d:
        sm2_d[file.header] = round(file[7].sum(),4)
        
with open('sm2_non_overl_summary.csv', 'a') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow([file.header])
    for key, value in sm2_d.items():
        writer.writerow([key, value])
        

#replace string with ''       
path='//vendiagram/mapped_unique/intersected/'
os.chdir(path)
path=os.getcwd()
    
intersect=[]
for i in list_data:
    if i.startswith('all_') or i.startswith('rounded_weighted_counts'):
        file=i.split('\n')
        file=file[0]
        intersect.append(file)   
    
for i in intersect:
    number=0
    for chunk in pd.read_csv(path+'/'+i,chunksize=1000000, 
                             sep= '\t', engine= 'python', 
                             header=None):
       
        number+=1
        print(number)
        sense=pd.DataFrame(chunk.loc[chunk[5]==chunk[12]])
        maker=pd.DataFrame(sense.loc[sense[13]=='maker'])

        replace=maker[10].replace(
            to_replace=(r'[-][A-Z]{2}[:][a-z]{4}[:]\d{1,4}'),
            value= '', regex=True)

        makers=maker
        makers[10]=replace #new column in maker 

        maker.to_csv('replaced_{}_{}'.format(i, number), 
                     sep="\t", index=False)    
    
#concat all files with same name 
for file in intersect:
    with open('all_repla_{}'.format(file), 'a+') as out:
        for fn in filenames:
                #print(file.header[26:-4])
                if file[28:] in fn:
                    #print('file.header', file[30:])
                    with open(fn) as f:
                        out.write(f.read())
                        
list_data=os.popen('ls').readlines() #list of filenames
print(list_data)     
    
replace=[]
for i in list_data:
    if i.startswith('all_repla_'):
        file=i.split('\n')
        file=file[0]
        replace.append(file)
#sort and keep only uniquely occuring sequences         
for i in replace:
    op=i.split('all_repla_all_intersected_sorted_mapped_')
    op=op[1]
    command="cat '{}' | cut -f11,18 | "+
    "sort -u > '{}/for_ven_maker_sort_uniq_{}'".format(path+'/'+i, 
                                               path, op
                           )
    subprocess.check_call(command, 
                          stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE, 
                          shell=True)
    print(command)
    
    


# In[ ]:


#write summary 
#how many of TE and other annotated elemnts are present 
#for stacked barplot 
import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob


os.chdir('/data1/AnLi_FM_MA/data/output_sm3_KD/normalized/')
path=os.getcwd()
intersect=glob('norm*bed.csv')

col_13=['schMed4_repeatMaskerDNA_', 'schMed4_repeatMaskerLINE_',
        'schMed4_repeatMaskerLTR_', 
       'schMed4_repeatMaskerUnknown_', 'maker']
col_14=['exon', 'five_prime_UTR', 'three_prime_UTR']
output_sense={} #dict
output_antisense={} #dict

for Te in col_13:
    if Te not in output_sense:
        output_sense[Te] = 0
    if Te not in output_antisense:
        output_antisense[Te] = 0
        
for kind in col_14:
    if kind not in output_sense:
        output_sense[kind] = 0
    if kind not in output_antisense:
        output_antisense[kind] = 0
        
summary_string=['summary_sm3']

with open('sm3_summary.csv', 'w') as csv_file:
    writer = csv.writer(csv_file)
    writer.writerow(summary_string)

for i in intersect: 

    chunk=pd.read_csv(path+'/'+i, 
                      sep= '\t', 
                      engine= 'python', 
                      header=None)

    sense=pd.DataFrame(chunk.loc[chunk[5]==chunk[12]]) #sense

    antisense=pd.DataFrame(chunk.loc[chunk[5]!=(chunk[12])]) #antisense
    #--------------------------------------------------------
    maker_sense=pd.DataFrame(sense.loc[sense[13]=='maker'])
    maker_antisense=pd.DataFrame(antisense.loc[antisense[13]=='maker'])
    #---------------------------------------------------------
    for Te in col_13:
        for key, val in output_sense.items():
                if Te == key:#output_sense.keys():
                    print(Te)
                    te=pd.DataFrame(sense.loc[sense[13]==Te]) 
                    #df nut mit bestimmnten TE
                    sum_17=te[17].sum()
                    print(sum_17)
                    output_sense[Te]=(round(sum_17, 4)) 
    #---------------------------------------------------------------
        for key, val in output_antisense.items():
                #print(key)
                if Te == key:#output_sense.keys():
                    print(Te)
                    te=pd.DataFrame(antisense.loc[antisense[13]==Te]) 
                    #df nut mit bestimmnten TE
                    sum_17=te[17].sum()
                    print(sum_17)
                    output_antisense[Te]=(round(sum_17, 4)) 
    #----------------------------------------------------------------

    for kind in col_14:
        for key, value in output_sense.items():
            if kind == key:
                Kind=pd.DataFrame(
                    maker_sense.loc[sense[14]==kind])
                sum_17_14=Kind[17].sum()
                output_sense[kind]=(round(sum_17_14, 4)) 
        #--------------------------------------------------------------
        for key, value in output_antisense.items():
            if kind == key:
                Kind=pd.DataFrame(
                    maker_antisense.loc[antisense[14]==kind])
                sum_17_14=Kind[17].sum()
                output_antisense[kind]=(round(sum_17_14, 4))  

    with open('sm3_summary.csv', 'a') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([i])
        writer.writerow(['sense'])
        for key, value in output_sense.items():

            writer.writerow([key, value])
        writer.writerow(['antisense'])
        for key, value in output_antisense.items():

            writer.writerow([key, value])


# In[ ]:


#cutoff 18 till 26 nts 
from Bio import SeqIO
import subprocess, os, sys
from glob import glob

fasta_file_dir=str(input('fastq files dir: '))
fasta_cutoff=str(input('path to new_dir: '))

def keep_RNA_smaller_26(fasta_file_dir, fasta_cutoff):
    """
    parse through fasta files
    """
    if not os.path.exists(fasta_cutoff): 
        os.mkdir(fasta_cutoff)    
    #-------------------------
    os.chdir(fasta_file_dir)
    os.getcwd()
    fasta_files=glob('rejected_*')
    #-----------------------
    for fasta in fasta_files:
        os.chdir(fasta_cutoff)
        file_txt=open('pot_miRNA_{}'.format(fasta), 'w')
        os.chdir(fasta_file_dir)
        file= SeqIO.parse(fasta, "fastq")
        for seq_record in file:
            os.chdir(fasta_cutoff)
            with open('pot_miRNA_{}'.format(fasta), 'a') as file_txt:
                if len(seq_record.seq)>=18 and 
                len(seq_record.seq)<=26:
                    SeqIO.write(seq_record, file_txt, 'fastq')
    print(fasta)
    
keep_RNA_smaller_26(fasta_file_dir, fasta_cutoff)


# In[ ]:


#run mappel.pl from mirdeep2 

from Bio import SeqIO
import subprocess, os, sys
from glob import glob
import optparse, shutil

fasta_file_dir=str(input('file to input fastq: ')) #//miRNA_output/
fasta_cutoff=str(input('path to output dir : ')) #/mapper_output
genome=str(input('path to genome : ')) #//data/output_sm3_KD/bowtie_index/

def run_mapper_pl(fasta_file_dir, fasta_cutoff):
    """
    starts mapper.pl from miRDeep2
        -e input file is fastq format
        -h parse to fasta format
        -j remove all entries that have a sequence that contains letters
        -m collapse reads (needed for miRDeep2.pl)
        -p genome (which genome to map to)
        -s print processed reads to this file
        -t print read mappings to this file
        -v outputs progress report
        
    input: fastq file 
           genome of the organsim 
    output:fasta file, mapped to genome 
           arf file -> needed for mirdeep2.pl (finding miRNAs)
    """
    if not os.path.exists(fasta_cutoff): #make dir for output
        os.mkdir(fasta_cutoff)    
    #-------------------------
    os.chdir(fasta_file_dir)
    os.getcwd()
    files=glob('pot_miRNA_*')
    fasta_files=[]
    for file in files:
        if file.endswith('.fq.fq'):
            fasta_files.append(file)
    #-----------------------
    for fasta in fasta_files:
        new_fasta=fasta.split('.fq.fq')[0]
        #print(fasta)
        command = "mapper.pl {} -e -h -j -m -p {} -s"+
            " {}.fa -t {}{}_reads_vs_mapper.arf -v ".format(
            fasta_file_dir+fasta,
            genome+'btin_Smed_g4',
            fasta_cutoff+'/'+new_fasta,
            fasta_cutoff+'/',
            new_fasta.split('pot_miRNA_rejected_reads_trimmed_')[1]
            
            ) 


        subprocess.check_call(command, 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(command)
    
run_mapper_pl(fasta_file_dir, fasta_cutoff)


# In[ ]:


#run mirdeepl.pl 
from Bio import SeqIO
import subprocess, os, sys
from glob import glob
import optparse, shutil

fasta_file_dir=str(input('path to input fasta and arf: ')) 
#//mapper_output/
genome=str(input('path to genome : '))
#//planarian_genome/
mature_miRNA_path=str(input('path to mature planarian miRNA file : '))
#//mature_miRNA_planaria  -> from mirbase
 
def move_files_n(new_fasta):
    """
    literally just moves files and folders
    into order to keep output files under input file names
    """
    os.chdir(fasta_file_dir)
    os.getcwd()
    input_files=[] #they should stay in the dir

    #---------------------------------------------------- 
    if not os.path.exists(fasta_file_dir+'/{}'.format(
        new_fasta.split('.fa')[0])): 
        os.mkdir(fasta_file_dir+'/{}'.format(
            new_fasta.split('.fa')[0])) 
        input_files.append(new_fasta.split('.fa')[0])
        #make dir for outputwith the name of the input file 
    os.chdir(fasta_file_dir)
    output_files=glob('*') #list of all files
    #print(output_files)
    for file in output_files:
        fasta=file #.split('pot_miRNA_rejected_reads_trimmed_')[1]
        if file.endswith('.fa') or file.endswith('.arf'): 
            
            input_files.append(file)

    files_to_move=[]
    for fil in output_files:
        if fil not in input_files:
            files_to_move.append(fil)
        elif fil.startswith('mirna_reslults'):
            files_to_move.append(fil)
        elif fil.startswith('pdfs_'):
            files_to_move.append(fil)
        #elif fil.startswith('Nr'):
            #files_to_move.append(fil)

    #----------------------------------------------------
    for item in files_to_move:
        if item not in input_files:# or os.path.isdir(item):
            shutil.move(item, 
                        fasta_file_dir+'/{}'.format(
                            new_fasta.split('.fa')[0]))
            print ("'{}' moved to: '{}'".format (
                item, './{}'.format(
                new_fasta.split('.fa')[0])))
            #print(item)
        
        elif item in input_files: 
            print (item +' not moved')
            pass


        
def run_mirdeep(fasta_file_dir, genome, mature_miRNA_path):
    """
    run mirdeep.pl 
    input 
        fasta file
        reference genome 
        arf file from mapper.pl
        mature miRNAs from mirbase
        none (could be mature miRNAs from other organsims)
    """ 
    #-------------------------
    os.chdir(fasta_file_dir)
    os.getcwd()
    fasta_files=glob('*.fa')
    arf_files=glob('*.arf')
    mature_miRNA='DNA_mature_mirna_mirbase.fa'
    #print(fasta_files)
    #print(arf_files)
    #-----------------------

    for arf in arf_files:        
        for fasta  in fasta_files: 
            if arf.split(
                '_reads_vs_mapper.arf')[0]+'.fa' == fasta.split(
                'trimmed_')[1]:
                #print(arf.split('_reads_vs_mapper.arf')[0])
                #print(fasta.split('trimmed_')[1])
                new_fasta=fasta.split(
                    'pot_miRNA_rejected_reads_trimmed_')[1]
                #----------------------------------------------------
                #print(new_fasta)
                c = "miRDeep2.pl {} {} {} {} "+
                "none none -t Planaria 2> {}.log ".format(
                    fasta_file_dir+fasta,
                    genome+'final_dd_Smed_g4.fa',
                    fasta_file_dir+arf,
                    mature_miRNA_path+'/'+mature_miRNA,
                    new_fasta.split('.fa')[0]

                    ) #name of genome is saved as string, 
                #change here if needed

                print(c)
                subprocess.check_call(c, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, 
                                      shell=True)

                move_files_n(new_fasta)
        
    
run_mirdeep(fasta_file_dir, genome, mature_miRNA_path)


# In[ ]:


#grep those bedfiles, which include the mature and novel miRNA sequence 


import subprocess, os, sys
import re
import optparse, shutil
from glob import glob

path_to_input_txt=input('Path to input txt files (unique KD/WT) : ') 
#//mapper_output/mirdeep-pl_output/insert/
path_to_input_bed2=input('Path to input bed2 files : ') 
#//bed_files/less_than_26/
output_path_mapped_bed2=input('Path to output w/ new dir name : ')
#/data/output_sm3_KD/grep_miRNA

def grep(path_to_input_txt, path_to_input_bed2, output_path_mapped_bed2):
    """
    input parameter
        txt file: only unique sequences from 
        either the KD or WT, which dont overlap in both (from vendiagram)
        bed2 file: mapped sequences to genome
        output path: dir for the output
    return
        bed2 file, includes only the sequences, 
        that are included in either the KD or WT
    """
    
    #----------------------------------------------------------
    if not os.path.exists(output_path_mapped_bed2): 
        os.mkdir(output_path_mapped_bed2)
    
    #----------------------------------------------------------
    os.chdir(path_to_input_txt)
    os.getcwd()
    txt=glob('*mature.txt')
    #------------------------------------
 
    print(txt)
    #----------------------------------------------------------
    os.chdir(path_to_input_bed2)
    os.getcwd()
    bed2=glob('lesser*.bed2')     
    #------------------------------------

    
    print(bed2)
    for file in txt:
        
        motif=file.split('_')[0]
        if 'x1' in file:#
        
            command = "grep -Fwf {} {} > {}/m_{}.bed2".format(
                path_to_input_txt+'/'+file, 
                path_to_input_bed2+'/lesser_26_d_'+motif+'.bed2',
                output_path_mapped_bed2, file.split('.txt')[0])
            sort="sort -k1,1 -k2,2n {}/m_{}.bed2 > {}/sort_m_{}.bed2".format(
                output_path_mapped_bed2,
                file.split('.txt')[0],output_path_mapped_bed2, 
                file.split('.txt')[0]) 
        else:
            command = "grep -Fwf {} {} > {}/m_{}.bed2".format(
                path_to_input_txt+'/'+file, 
                path_to_input_bed2+'/lesser_26_'+motif+'.bed2',
                output_path_mapped_bed2, file.split('.txt')[0])
            sort="sort -k1,1 -k2,2n {}/m_{}.bed2 > {}/sort_m_{}.bed2".format(
                output_path_mapped_bed2,
                file.split('.txt')[0],output_path_mapped_bed2, 
                file.split('.txt')[0]) 


                                                                  
        subprocess.check_call(command, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(command)
        subprocess.check_call(sort, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        print(sort)

grep(path_to_input_txt, path_to_input_bed2, output_path_mapped_bed2)  


# In[ ]:


import subprocess,sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

path_to_bed=input('Path to bed file: ') 
#//vendiagram/mapped_unique/intersected/
path_protar=input('Path to targeted and protected piRNAs: ') 
#/data1/AnLi_FM_MA/data/output_sm3_KD/protected_targeted_piRNAs/

def output_folder(path_t):
    """
    creates new output folder 
    """
    if not os.path.exists(path_t):
        os.makedirs(path_t)
    os.chdir(path_t)
    path_norm=os.getcwd()
    return path_norm

def exon_wc():
    """
    get -> files (bed files) after intersecting with 
    annoated elements start with all_repl
    read in files and keep only those mapped to exons 
    output file contains only name *and* weighted count (wc)
    """
    path_protar='/data1/AnLi_FM_MA/data/output_sm3_KD/grep_miRNA/'
    path_op=output_folder(path_protar+'/exon_bed_files_for_venn')
    #get files
    os.chdir(path_to_bed) #change current working dir
    os.getcwd()
    all_repl=glob('all_repl*.bed')
    for file in all_repl:
        for chunk in pd.read_csv(path_to_bed+file, 
                                 sep='\t', 
                                 engine='python', 
                                 header=None, 
                                 usecols=[10,14,17], 
                                 chunksize=100000): #10:JG_000 14:exon 17:wc
            exon=pd.DataFrame(chunk.loc[chunk[14]=='exon']) #only exons
            exon=exon.drop([14], axis=1)

            if not os.path.exists('{}/name_wc_{}'.format(path_op, 
                file.split(
                    'all_repla_all_intersected_sorted_mapped_')[1])):
                exon.to_csv('{}/name_wc_{}'.format(
                path_op, 
                file.split(
                    'all_repla_all_intersected_sorted_mapped_')[1]), 
                sep='\t', header=None, index=False)
            if os.path.exists('{}/name_wc_{}'.format(path_op, 
                file.split(
                    'all_repla_all_intersected_sorted_mapped_')[1])):
                exon.to_csv('{}/name_wc_{}'.format(
                path_op, 
                file.split(
                    'all_repla_all_intersected_sorted_mapped_')[1]), 
                sep='\t', header=None, index=False, mode='a')
exon_wc()     

def exon():
    """
    get -> files (bed files) after intersecting with 
    annoated elements start with all_repl
    save files containing *only* the name 
    """
    path_op=output_folder(path_protar+'/exon_bed')
    #get files
    os.chdir(path_to_bed) #change current working dir
    os.getcwd()
    all_repl=glob('all_repl*.bed')
    for file in all_repl:
        for chunk in pd.read_csv(path_to_bed+file, sep='\t', engine='python', 
                                 header=None, 
                                 usecols=[10,14], 
                                 chunksize=100000): #10:JG_000 14:exon
            exon=pd.DataFrame(chunk.loc[chunk[14]=='exon']) #only exons
            exon=exon[10]


            if not os.path.exists('{}/name_{}'.format(path_op, 
                file.split('all_repla_all_intersected_sorted_mapped_')[1])):
                exon.to_csv('{}/name_{}'.format(
                path_op, file.split(
                    'all_repla_all_intersected_sorted_mapped_')[1]), 
                sep='\t',
                header=None, index=False)
            if os.path.exists('{}/name_{}'.format(path_op, 
                file.split('all_repla_all_intersected_sorted_mapped_')[1])):
                exon.to_csv('{}/name_{}'.format(
                path_op, 
                file.split('all_repla_all_intersected_sorted_mapped_')[1]), 
                sep='\t', header=None, index=False, mode='a')
exon()     

def open_excel():
    """
    leave only names in csv 
    """
    path_op=output_folder(path_protar+'/protar')
    os.chdir(path_protar+'/protar') #change current working dir
    os.getcwd()
    excel=glob('*.xlsx')
    for file in excel:
        df=pd.read_excel(file, sheet_name='Sheet1',
                         header=0, index_col=None)
        #print (df.head(5))
        save=df['target_id']
        #print(save.head(5))
        save.to_csv('{}/names_{}.bed'.format(path_op, 
                    file.split('.xlsx')[0]), 
                    header=None, sep='\t',
                    index=False)
        
open_excel()



def in_reps():
    """
    fill dataframe with all three reps 
    """
    path_op=output_folder('/data/output_sm3_KD/grep_miRNA/names_norm'+
                          '/thresh_allreps')
    os.chdir('/data/output_sm3_KD/grep_miRNA/names_norm/') 
    path_protar='/data/output_sm3_KD/grep_miRNA/names_norm/'
    
    #path_op=output_folder(path_protar+'/thresh_allreps')
    #os.chdir(path_protar+'/normalized') #change current working dir
    os.getcwd()
    exon=glob('norm*.csv')
    
    sm3_X1=pd.DataFrame(columns=['names','1', '2', '3',
                                 'wc_1', 'wc_2', 'wc_3'])
    sm3_Xins=pd.DataFrame(columns=['names','1', '2', '3',
                                   'wc_1', 'wc_2', 'wc_3'])
    sm3_unc_X1=pd.DataFrame(columns=['names','1', '2', '3',
                                     'wc_1', 'wc_2', 'wc_3'])
    sm3_unc_Xins=pd.DataFrame(columns=['names','1', '2', '3',
                                       'wc_1', 'wc_2', 'wc_3'])
    
    
    for file in exon:
        #columnnames
        p=file.split('_')
        if len(p)==4:
            columnname=p[2].split('sm3')[1]
            #norm unc sm33 xins.bed #['', 3]
        else:
            columnname=p[1].split('sm3')[1]
            #norm sm33 xins.bed
        print (columnname)
        df=pd.read_csv(path_protar+file, sep='\t', 
                       engine='python', header=None, usecols=[0,2])
        #df=pd.read_csv(path_protar+'/normalized/'+file, sep='\t', 
                       #engine='python', header=None, usecols=[0,3])
        df=df.drop_duplicates(subset=[0], 
                              keep='first', ).reset_index(drop=True)
        #print(df.head(2))
        df['names']=df[0] #have two columns -> zero and 'names'
        df.rename(columns={0:columnname}, inplace=True) 
        #columnname -> number of replicate
        df.rename(columns={2:'wc_'+columnname}, inplace=True) 
        #for wc -> wc plus numb of rep

        if 'norm_sm3' in file and 'X1' in file:
            #print(df.head(2))
            sm3_X1=sm3_X1.append(df, ignore_index=True)

        elif 'norm_sm3' in file and 'Xins' in file:
            sm3_Xins=sm3_Xins.append(df, ignore_index=True)

        elif 'norm_unc_sm3' in file and 'X1' in file:
            sm3_unc_X1=sm3_unc_X1.append(df, ignore_index=True)
            
        elif 'norm_unc_sm3' in file and 'Xins' in file:
            sm3_unc_Xins=sm3_unc_Xins.append(df, ignore_index=True)
            
    print(sm3_X1.head(2))
    sm3_X1.to_csv('{}/123_sm3_X1.csv'.format(path_op), 
                    header=True, sep='\t', index=False)
    sm3_Xins.to_csv('{}/123_sm3_Xins.csv'.format(path_op), 
                    header=True, sep='\t', index=False)
    sm3_unc_X1.to_csv('{}/123_unc_X1.csv'.format(path_op), 
                    header=True, sep='\t', index=False)
    sm3_unc_Xins.to_csv('{}/123_unc_Xins.csv'.format(path_op), 
                    header=True, sep='\t', index=False)
    
in_reps()

def thresh():
    """
    drop all names, if they occur less than twice in all three reps
    count the sum of wc and drop all below upper quartile number
    save descriptions of df (all) and after threshold
    """
    path_protar='//data/output_sm3_KD/grep_miRNA/names_norm/' #temp
    
    path_op=output_folder(path_protar+'/thresh_allreps')
    
    
    
    os.chdir(path_protar+'/thresh_allreps') #change current working dir
    os.getcwd()
    files=glob('123*.csv')
    #later for aggregate
    col={'1':'first', '2': 'first', '3':'first', 
         'wc_1':'first', 'wc_2':'first', 'wc_3':'first'}
    
    for file in files:
        print(file)
        df=pd.read_csv(path_protar+'thresh_allreps/'+file, sep='\t', 
                        engine='python')
        df=df.groupby(df['names']).aggregate(col)
        df=df.dropna(subset=['1', '2', '3'], thresh=2)
        df=df.reset_index()
        df=df.drop(['1', '2','3'], axis=1)

        ###keep only those appearing in at least two reps
        col_list=list(df)
        col_list.remove('names')
        #print(col_list)
        df['sum_wc']=df[col_list].sum(axis=1)
        
        df_upper_quartile=df.loc[:,'sum_wc'].quantile(q=0.75)#, axis=1)
        print('up qu:', df_upper_quartile) #for Iana
        df_thresh5=df.drop(df[df['sum_wc']<df_upper_quartile].index)
        df_13=df.drop(df[df['sum_wc']<1.3].index)#, inplace=True)

        #print(df_thresh5.head(2))


        
        df_13=df_13.drop(['wc_1', 'wc_2', 'wc_3'], axis=1)
        df_13.to_csv('{}/thresh_1.3_{}'.format(path_op, file.split('123_')[1]), 
                      header=True, sep='\t', index=False)

        
thresh()

def protar_py():
    """
    check for overlaps 
    """
    path_op=output_folder(path_protar+'/intersected')
    os.chdir(path_protar+'/all_reps_sm3') #change current working dir
    os.getcwd()
    exon=glob('*')
    #print(exon)
    os.chdir(path_protar+'/protar') #change current working dir
    os.getcwd()
    excel=glob('*.bed')
    one=pd.read_csv(excel[0], sep='\t', engine='python', 
                                 header=None)
    two=pd.read_csv(excel[1], sep='\t', engine='python', 
                                 header=None)
    one_v=excel[0].split('_')[1]
    two_v=excel[1].split('_')[1]
    
    
    for file in exon:
        ff=file.split('.csv')[0]
        print(ff)
        df=pd.read_csv(path_protar+'/all_reps_sm3/'+file, 
                       sep='\t',
                       engine='python', 
                       usecols=['names'])
        #print(df[0].head())
        df=df.drop_duplicates( keep='first', ).reset_index(drop=True)
        print(df.head())
        df[one_v]=one
        df[two_v]=two
        df['unique_{}'.format(one_v)]=df[one_v][df[one_v].isin(df['names'])]
        df['unique_{}'.format(two_v)]=df[two_v][df[two_v].isin(df['names'])]
        print(df)
        #df['unique_{}'.format(two_v)]=df[two_v][~df[two_v].isin(df[0])]mit schlange
        #---> nur die aufnehmen, 
        #die nicht in df[two_v] vorkommen!genau das gegenteil
       
        unique_protected=df['unique_{}'.format(one_v)].dropna(how='all')
        unique_pingpong=df['unique_{}'.format(two_v)].dropna(how='all')
        #print(unique_protected)
        #print(unique_pingpong)
        df.to_csv('{}/whole_df_{}.csv'.format(path_op, 
                    ff), 
                    header=True, sep='\t', index=False)
        unique_protected.to_csv('{}/protected_{}.csv'.format(
            path_op, ff))
        unique_pingpong.to_csv('{}/pingpong_{}.csv'.format(
            path_op, ff))
protar_py()      


# In[ ]:


#make dataframe fro edge R 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob


import matplotlib.pyplot as plt 

def get_files_2(path_2, motif)-> list:
    """
    load all files names into a list
    """
    os.chdir(path_2)
    path_2=os.getcwd()
    files=glob(motif)
    return(files)

def output_folder(path):
    """
    creates new output folder 
    """
    if not os.path.exists(path):
        os.makedirs(path)
    os.chdir(path)
    path=os.getcwd()
    unique_output=path
    return unique_output
def plot_miRNA(data, name,path):
    """
    barplot
    """
    summary=data
    p=summary.plot.bar(y='sum', 
                       x='name' ,
                       title=name, legend=False)
    p.set_ylabel('Total reads', fontsize=14)
    p.set_xlabel('miRNAs', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    fig_name=name.replace(' ', '_')
    plt.savefig(path+'/{}_sm3.pdf'.format(fig_name),
                bbox_inches='tight')
    plt.show()
    plt.close()
    
def get_known_unique(): #here
    """
    create dataframes from unique *known* miRNAs
    
    1-> sequence unc1 unc2 unc3 sm31 sm32 sm33 
    2-> leave only sequences, which are present 
    at least in two reps, discard rest 
    """
    unique_path=output_folder(path+'/edgeR_dataframe')
    files= get_files_2(path+'/unique_miRNA', '*known_mature.bed2')
    print(path)
    #unique_path=output_folder('/Users/annaliznar/Downloads/DELETEMEE')
    #files= get_files_2('/Users/annaliznar/Downloads/', 'unique_*')
    print(files)
    
    x1=pd.DataFrame(columns=['sequences', 'unc1', 'unc2', 'unc3',
                             'sm31', 'sm32', 'sm33'])
    xins=pd.DataFrame(columns=['sequences', 'unc1', 'unc2', 'unc3',
                               'sm31', 'sm32', 'sm33'])
    merged_xins=pd.DataFrame(columns=['sequences', 'unc1', 'unc2', 
                                      'unc3', 'sm31', 'sm32', 'sm33'])
    
    for file in files:
        df=pd.read_csv(path+'/unique_miRNA/'+file, 
                      sep='\t',
                      engine='python',
                      usecols=[6,7])
        print(df.head(2))
        df.rename(columns={'6':'sequences'}, inplace=True)
        df['7']=df['7'].round(3)
        motif=file.split('_')[1] #name of file 
        #print(motif.split('-')[1])
        columnname=motif.split('-')[1] #nmae of column 
        
        if 'x1' in motif and motif in file:
            df.rename(columns={'7':columnname}, inplace=True)
            x1=x1.append(df, ignore_index=True)
            
        elif 'xins' in motif and motif in file:
            df.rename(columns={'7':columnname}, inplace=True)
            xins=xins.append(df, ignore_index=True)
    
    col={'unc1':'sum', 'unc2': 'sum', 'unc3':'sum', 
         'sm31':'sum', 'sm32':'sum', 'sm33':'sum'}
    #,'sequences':'first'}
    grouped_xins=xins.groupby(xins['sequences']).aggregate(col)
    grouped_x1=x1.groupby(xins['sequences']).aggregate(col)
    print(grouped_xins[0:10])
    
    #replace 0 with NAN 
    
    grouped_xins=grouped_xins.replace(0, np.nan)
    grouped_x1=grouped_x1.replace(0, np.nan)
    print(grouped_x1[0:10])
    
    #drop rows with at more than one NAN 
    
    grouped_xins=grouped_xins.dropna(subset=['unc1',
                                             'unc2', 'unc3'], thresh=2)
    grouped_xins=grouped_xins.dropna(subset=['sm31',
                                             'sm32', 'sm33'], thresh=2)
    print(grouped_xins[0:10])
    
    grouped_x1=grouped_x1.dropna(subset=['unc1', 
                                         'unc2', 'unc3'], thresh=2)
    grouped_x1=grouped_x1.dropna(subset=['sm31', 
                                         'sm32', 'sm33'], thresh=2)
    print(grouped_xins[0:10])    
    
    grouped_x1.to_csv(
        '{}/edgeR_dataframe/known_miRNAs_for_edgeR_x1.csv'.format(
        path), sep='\t', 
                 header=True, 
                 index=True)
    
    grouped_xins.to_csv(
        '{}/edgeR_dataframe/known_miRNAs_for_edgeR_xins.csv'.format(
        path), sep='\t', 
                 header=True, 
                 index=True)
            
    
    return grouped_xins, grouped_x1
get_known_unique()    

#get all miRNAs (known and novel)

 def get_unique():
    """
    create dataframes 
    1-> sequence unc1 unc2 unc3 sm31 sm32 sm33 
    2-> leave only sequences, which are present at 
    least in two reps, discard rest 
    """
    unique_path=output_folder(path+'/edgeR_dataframe')
    files= get_files_2(path+'/unique_miRNA', 'uni*bed2')
    
    #unique_path=output_folder('/Users/annaliznar/Downloads/DELETEMEE')
    #files= get_files_2('/Users/annaliznar/Downloads/', 'unique_*')
    #print(files)
    
    
    x1=pd.DataFrame(columns=['sequences', 'unc1', 'unc2',
                             'unc3', 'sm31', 'sm32', 'sm33'])
    
    xins=pd.DataFrame(columns=['sequences', 'unc1', 'unc2',
                               'unc3', 'sm31', 'sm32', 'sm33'])
    merged_xins=pd.DataFrame(columns=['sequences', 'unc1', 'unc2',
                                      'unc3', 'sm31', 'sm32', 'sm33'])
    #x2=x2.set_index('sequences')

    
    for file in files:
        df=pd.read_csv(path+'/unique_miRNA/'+file, 
                      sep='\t',
                      engine='python',
                      usecols=[6,7])
        df.rename(columns={'6':'sequences'}, inplace=True)
        df['7']=df['7'].round(3)
        motif=file.split('_')[1] #name of file 
        #print(motif.split('-')[1])
        columnname=motif.split('-')[1] #nmae of column 
        
        if 'x1' in motif and motif in file:
            df.rename(columns={'7':columnname}, inplace=True)
            x1=x1.append(df, ignore_index=True)
            
        elif 'xins' in motif and motif in file:
            df.rename(columns={'7':columnname}, inplace=True)
            xins=xins.append(df, ignore_index=True)
    #print(len(x1['sequences']))  
    
    #x1=x1.drop_duplicates(subset='sequences',keep='first')
    #xins=xins.drop_duplicates(subset='sequences',keep='first')
    #print(len(x1['sequences'])) 
    #print(x1.head(5))
    #print(xins.head(5))
    
    
    #print(len(xins['sequences']))   'sequences':'first', 
    col={'unc1':'sum', 'unc2': 'sum', 'unc3':'sum', 
         'sm31':'sum', 'sm32':'sum', 'sm33':'sum'}#,'sequences':'first'}
    grouped_xins=xins.groupby(xins['sequences']).aggregate(col) 
    #.apply(pd.DataFrame)
    grouped_x1=x1.groupby(xins['sequences']).aggregate(col)
    print(grouped_xins[0:10])
    
    #replace 0 with NAN 
    
    grouped_xins=grouped_xins.replace(0, np.nan)
    grouped_x1=grouped_x1.replace(0, np.nan)
    print(grouped_x1[0:10])
    
    #drop rows with at more than one NAN 
    
    grouped_xins=grouped_xins.dropna(subset=['unc1', 
                                             'unc2', 'unc3'], thresh=2)
    grouped_xins=grouped_xins.dropna(subset=['sm31', 
                                             'sm32', 'sm33'], thresh=2)
    print(grouped_xins[0:10])
    
    grouped_x1=grouped_x1.dropna(subset=['unc1', 
                                         'unc2', 'unc3'], thresh=2)
    grouped_x1=grouped_x1.dropna(subset=['sm31', 
                                         'sm32', 'sm33'], thresh=2)
    print(grouped_xins[0:10])    
    
    grouped_x1.to_csv(
        '{}/edgeR_dataframe/for_edgeR_x1.csv'.format(
        path), sep='\t', 
                 header=True, 
                 index=True)
    
    grouped_xins.to_csv(
        '{}/edgeR_dataframe/for_edgeR_xins.csv'.format(
        path), sep='\t', 
                 header=True, 
                 index=True)
    return grouped_xins, grouped_x1
get_unique()


# In[ ]:


#calc norm factor with total reads 
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob


import matplotlib.pyplot as plt 


os.chdir('//data/output_sm3_KD/grep_miRNA/')
path=os.getcwd()

def get_files(motif)-> list:
    """
    load all files names into a list
    """
    #global path
    
    os.chdir('/data/output_sm3_KD/grep_miRNA/')
    path=os.getcwd()
    files=glob(motif)
    return(files)

def get_files_2(path_2, motif)-> list:
    """
    load all files names into a list
    """
    os.chdir(path_2)
    path_2=os.getcwd()
    files=glob(motif)
    return(files)

def output_folder(path):
    """
    creates new output folder 
    """
    if not os.path.exists(path):
        os.makedirs(path)
    unique_output=path
    return unique_output

def unique_df_pandas(files):
    """
    load data into df with pandas 
    """
    unique_output = output_folder(path
                                  +'/unique_miRNA')
    print(unique_output)
    for file in files:
        #sub_csv=csv[11]
        df=pd.read_csv(path+'/'+file,sep='\t',  
                       engine= 'python', 
                       header=None)
        #print(df.head(5))
        unique = df.drop_duplicates(subset=[6], 
                                    keep='first')
        #print(unique.head(5))
        
        if not os.path.exists(unique_output+
                              'unique_'+file.split(
                                  'weighted_counts_sort_m_')[1]):
            unique.to_csv('{}/unique_{}'.format(
                 unique_output, 
                 file.split('weighted_counts_sort_m_')[1]), 
                 sep='\t', 
                 header=True, 
                 index=False)
            #pass

files=get_files('weighted*')        
unique_df_pandas(files) 

def plot_miRNA(data, name,path):
    """
    barplot
    """
    summary=data
    p=summary.plot.bar(y='sum', 
                       x='name' ,
                       title=name, legend=False)
    p.set_ylabel('Total reads', fontsize=14)
    p.set_xlabel('miRNAs', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    fig_name=name.replace(' ', '_')
    plt.savefig(path+'/{}_sm3.pdf'.format(fig_name),
                bbox_inches='tight')
    plt.show()
    plt.close()
    
def total_count():
    """
    count the total count of miRNA
    """
    unique_path = output_folder(path+'/total_count')
    files=get_files_2(path+'/unique_miRNA', 'unique*')
    summary_n=pd.DataFrame( columns=['name','sum'])
    summary_k=pd.DataFrame( columns=['name','sum'])
    for file in files:
        
        if 'novel' in file:
            df_n=pd.read_csv(path+'/unique_miRNA/'+file,
                             sep='\t', 
                             engine= 'python',
                             header=None)
            summary_n=summary_n.append({'name':file.split('.bed')[0],
                                       'sum':df_n[7].sum()}, 
                                       ignore_index=True)
        elif 'known' in file:
            df_k=pd.read_csv(path+'/unique_miRNA/'+file, 
                             sep='\t', 
                             engine= 'python', 
                             header=None)
            summary_k=summary_k.append({'name':file.split('.bed')[0],
                                       'sum':df_k[7].sum()}, 
                                       ignore_index=True)

    summary_n['sum']  = summary_n['sum'].astype(float)
    summary_n['sum']  = summary_n['sum'].round(3)
    summary_n['name'] = summary_n['name'].replace('unique_',
                                                  '', regex=True)
    summary_n['name'] = summary_n['name'].replace('_mature',
                                                  '', regex=True)
    summary_n['name'] = summary_n['name'].replace('_novel',
                                                  '', regex=True)
    summary_n=summary_n.sort_values('name', ascending=False)
    summary_n=summary_n.reset_index(drop=True)
    
    summary_k['sum']  = summary_k['sum'].astype(float)
    summary_k['sum']  = summary_k['sum'].round(3)    
    summary_k['name'] = summary_k['name'].replace('unique_',
                                                  '', regex=True)
    summary_k['name'] = summary_k['name'].replace('_mature',
                                                  '', regex=True)
    summary_k['name'] = summary_k['name'].replace('_known',
                                                  '', regex=True)
    summary_k=summary_k.sort_values('name', ascending=False)
    summary_k=summary_k.reset_index(drop=True)
    
    summary_k.to_csv('{}/total_reads_known_miRNAs_sm3.csv'.format(
        unique_path),
                 sep='\t', 
                 header=True, 
                 index=False)
    print(summary_n)

    plot_miRNA(summary_n, "Total reads of novel miRNAs", 
               unique_path)
    summary_n.to_csv('{}/total_reads_novel_miRNAs_sm3.csv'.format(
        unique_path), 
                 sep='\t', 
                 header=True, 
                 index=False)
    
    print(summary_k)
    
    plot_miRNA(summary_k, 'Total reads of known miRNAs',
               unique_path)

    
total_count()

def mean_miRNA():
    """
    count mean for all replicates
    """
    unique_path = output_folder(path+'/mean_repl')
    files = get_files_2(path+'/unique_miRNA', 'unique*')
    summary=pd.DataFrame( columns=['name','sum'])
    for file in files:
        df=pd.read_csv(path+'/unique_miRNA/'+file,
                             sep='\t', 
                             engine= 'python',
                             header=None)
        summary=summary.append({'name':file.split('.bed')[0],
                                       'sum':df[7].sum()}, 
                                       ignore_index=True)
    summary['sum']  = summary['sum'].round(3)
    summary['name'] = summary['name'].replace('unique_',
                                              '', regex=True)
    summary['name'] = summary['name'].replace(
        {'_novel_mature': '',
         '_known_mature':''}, 
          regex=True)
    summary=summary.groupby('name').mean().reset_index()
    plot_miRNA(summary, 'Mean of novel and known miRNAs', 
               unique_path)
    summary.to_csv('{}/mean_miRNAs_sm3.csv'.format(
        unique_path), 
                 sep='\t', 
                 header=False, 
                 index=False)
    
    return summary
        
mean_miRNA()

def calc_norm_fac_known(): #new -> only knwon miRNAs
    """
    calculate normalization factor 
    """
    unique_path = output_folder(path+'/norm_fact_total_miRNAs')
    #files = get_files_2(path+'/mean_repl', '*csv')
    files = get_files_2(path+'/total_count', '*known_miRNAs_sm3.csv')

    for file in files:
        df=pd.read_csv(path+'/total_count/'+file,
                             sep='\t', 
                             engine= 'python')#,
                             #header=None)
    value_unc_x1=df.iloc[n,1]
    #print(value_unc_x1)
    value_kd_x1=df.iloc[n+3,1]
    #print(value_kd_x1)
    value_unc_xins=df.iloc[n+6,1]
    #print(value_unc_xins)
    value_kd_xins=df.iloc[n+9,1]
    #print(value_kd_xins)
    
    df['fac_unc_x1']=value_unc_x1/df.loc[:,'sum']
    #print(df)
    df['fac_unc_x1']=df.loc[0:5,'fac_unc_x1']
    #print(df)
    df['fac_kd_x1']=value_kd_x1/df.loc[:,'sum']
    #print(df)
    df['fac_kd_x1']=df.loc[0:5,'fac_kd_x1']
    #print(df)
    df['fac_unc_xins']=value_unc_xins/df.loc[:,'sum']
    #print(df)
    df['fac_unc_xins']=df.loc[6:11,'fac_unc_xins']
    #print(df)
    df['fac_kd_xins']=value_kd_xins/df.loc[:,'sum']
    #print(df)
    df['fac_kd_xins']=df.loc[6:11,'fac_kd_xins']
    #print(df)
    #factor['factor']=factor['factor1'].combine_first(factor['factor2'])
    df['unc_factor']=df['fac_unc_x1'].combine_first(df['fac_unc_xins'])
    df['kd_factor']=df['fac_kd_x1'].combine_first(df['fac_kd_xins'])
    #print(df)
    df=df.drop(['sum','fac_unc_x1','fac_kd_x1',
                'fac_unc_xins','fac_kd_xins'], axis=1)
    print(df)
    #x1
    x1=df.loc[6:11, :]
    print(x1)
    #xins
    xins=df.loc[0:5,:]
    print(xins)
    
    #df.to_csv('{}/norm_fac_by_miRNA_librarysize.csv'.format(unique_path), 
     #            sep='\t', 
      #           index=False)
    p=df.plot.bar(y=['unc_factor', 'kd_factor'], 
                    x='name' ,
                    title='Normalization Factors by miRNA library size', 
                    legend=True,
                    rot=28)
    p.set_ylabel('Normalization Factors', fontsize=14)
    p.set_xlabel('sample', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    
    #plt.savefig(unique_path+'/norm_fact_sm3_both_unc_kd.pdf',
     #           bbox_inches='tight')
    #plt.show()
    plt.close() 
    
    k=df.plot.bar(y='unc_factor', 
                    x='name' ,
    title='Normalization Factors by miRNA library size', 
                    legend=False,
                    rot=28,
                 ylim=(0,1.3))
    k.set_ylabel('Normalization Factors', fontsize=14)
    k.set_xlabel('sample', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    #plt.figure(figsize=(10,4))
    #plt.bar(x=1,width=2, align='center', height=0.3)
    for i, v in enumerate(df['unc_factor']):
        #print(i, v)
        #[p.text(v, i, '{:.2f}'.format(v))
        k.text(s='{:.1f}'.format(v), x=i, y=v+0.05, color="k",
           verticalalignment="center", 
           horizontalalignment="center", 
           size=14)
    
    #plt.savefig(unique_path+'/norm_fact_sm3_unc.pdf',
     #           bbox_inches='tight')
    #plt.show()
    plt.close() 
    
    
    k=x1.plot.bar(y='unc_factor', 
                    x='name' ,
    title='Normalization Factors by miRNA library size in X1', 
                    legend=False,
                    rot=28,
                 ylim=(0,1.3))
    k.set_ylabel('Normalization Factors', fontsize=14)
    k.set_xlabel('sample', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    #plt.figure(figsize=(10,4))
    #plt.bar(x=1,width=2, align='center', height=0.3)
    for i, v in enumerate(x1['unc_factor']):
        #print(i, v)
        #[p.text(v, i, '{:.2f}'.format(v))
        k.text(s='{:.2f}'.format(v), x=i, y=v+0.05, color="k",
           verticalalignment="center", 
           horizontalalignment="center", 
           size=14)
    
    #plt.savefig(unique_path+'/norm_fact_sm3_unc.pdf',
     #           bbox_inches='tight')
    plt.show()
    plt.close() 
    
    k=xins.plot.bar(y='unc_factor', 
                    x='name' ,
                    title='Normalization Factors by miRNA library size in Xins', 
                    legend=False,
                    rot=28,
                 ylim=(0,1.3))
    k.set_ylabel('Normalization Factors', fontsize=14)
    k.set_xlabel('sample', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    #plt.figure(figsize=(10,4))
    #plt.bar(x=1,width=2, align='center', height=0.3)
    for i, v in enumerate(xins['unc_factor']):
        #print(i, v)
        #[p.text(v, i, '{:.2f}'.format(v))
        k.text(s='{:.2f}'.format(v), x=i, y=v+0.05, color="k",
           verticalalignment="center", 
           horizontalalignment="center", 
           size=14)
    
    #plt.savefig(unique_path+'/norm_fact_sm3_unc.pdf',
     #           bbox_inches='tight')
    plt.show()
    plt.close() 
    
calc_norm_fac_known()

#from edge r

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob


import matplotlib.pyplot as plt 


def bar_plot_norm_factor():
    """
    calculate normalization factor 
    """
    path='/data1/AnLi_FM_MA/data/output_sm3_KD/grep_miRNA/'
    path = output_folder(path+'/norm_fact_known_miRNAs')
    #path='/data1/AnLi_FM_MA/data/output_sm3_KD/grep_miRNA/'
    
    df1=pd.DataFrame(columns=['sample', 'factor'], index=range(1,7), 
                    data={'sample': ['unc1_x1','unc2_x1','unc3_x1',
                                    'sm31_x1','sm32_x1','sm33_x1'],
                                     
                         'factor': [0.9339488,1.2987569,1.9585996,
                                         0.6989799,0.7723046,0.7797405]})
                                         
    df2=pd.DataFrame(columns=['sample', 'factor'], index=range(1,7), 
                     data={'sample': ['unc1_xins','unc2_xins','unc3_xins',
                                    'sm31_xins','sm32_xins','sm33_xins'],
                            'factor':[0.8932141,0.9148167,0.8994626,
                                         1.0698291,1.2064300,1.0541700]})
    print(df1)
    
    factor=df1
    
    p=factor.plot.bar(y='factor', 
                    x='sample' ,
                    title='Normalization Factors known miRNAs edgeR', 
                    legend=False,
                    rot=20, ylim=(0,4))
    p.set_ylabel('Normalization Factors', fontsize=14)
    p.set_xlabel('sample', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    for i, v in enumerate(factor['factor']):
        #print(i, v)
        #[p.text(v, i, '{:.2f}'.format(v))
        p.text(s='{:.2f}'.format(v), x=i, y=v+0.3, color="k",
           verticalalignment="top", 
           horizontalalignment="center", 
           size=14)
    #fig_name=name.replace(' ', '_')
    plt.savefig(path+'/norm_fact_sm3_X1_edgeRRR.pdf',
                bbox_inches='tight')
    plt.show()
    plt.close()
    factor.to_csv('{}/norm_f_X1_sm3.csv'.format(path), 
                 sep='\t',  
                 index=False)
    
    ##############################################################
    factor=df2
    
    p=factor.plot.bar(y='factor', 
                    x='sample' ,
                    title='Normalization Factors known miRNAs edgeR', 
                    legend=False,
                    rot=20, ylim=(0,1.5))
    p.set_ylabel('Normalization Factors', fontsize=14)
    p.set_xlabel('sample', fontsize=14)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    for i, v in enumerate(factor['factor']):
        #print(i, v)
        #[p.text(v, i, '{:.2f}'.format(v))
        p.text(s='{:.2f}'.format(v), x=i, y=v+0.13, color="k",
           verticalalignment="top", 
           horizontalalignment="center", 
           size=14)
    #fig_name=name.replace(' ', '_')
    plt.savefig(path+'/norm_fact_sm3_XINS_edgeRRR.pdf',
                bbox_inches='tight')
    plt.show()
    plt.close() 
    
    
    
    factor.to_csv('{}/norm_f_XINS_sm3.csv'.format(path), 
                 sep='\t',  
                 index=False)
                  
bar_plot_norm_factor()  


# In[ ]:


import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_files_2(path_2, motif)-> list:
    """
    load all files names into a list
    """
    os.chdir(path_2)
    path_2=os.getcwd()
    files=glob(motif)
    return(files)
#global
path='/data/output_sm3_KD/grep_miRNA/exon_bed_files_for_venn/'

def output_folder(path_t):
    """
    creates new output folder 
    """
    if not os.path.exists(path_t):
        os.makedirs(path_t)
    os.chdir(path_t)
    path_norm=os.getcwd()
    return path_norm

#x1
unc1_x1= 1.0463347231012188
unc2_x1= 1.1927652214650772
unc3_x1= 1.0
sm31_x1= 0.4343830505625732
sm32_x1=0.5353575507025394
sm33_x1=0.4076869396741633
#xins
unc1_xins= 0.7748175781274539
unc2_xins= 0.9897655890200868 
unc3_xins= 1.0
sm31_xins= 0.43128541413492033
sm32_xins= 0.5288086443406546
sm33_xins= 0.43957917427670845


def normalize_files():
    """
     load bed2 files 
     create new column with weighted and normalized counts
    """
    path_nor=output_folder('/data1/AnLi_FM_MA/data/output_sm3_KD/'+
                           '/normalized_total_miRNA')
    files= get_files_2(path, 'all_inter*bed')
    
    print (files)
    
    for file in files:            
        for chunk in pd.read_csv(path+file, sep='\t', engine='python', 
                                 header=None, 
                                 usecols=[5,6,12,13,14,17], 
                                 chunksize=100000):
        #df=pd.read_csv(path+file, sep='\t', engine='python',
        #header=None, usecols=[5,6,12,13,14,17])
            df=chunk
            
            if 'unc_sm31_X1' in file:
                df['norm']=df[17].multiply(unc1_x1)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                             index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                              index=False,
                             mode='a')

            elif 'unc_sm32_X1' in file:
                df['norm']=df[17].multiply(unc2_x1)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None,
                             index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                              index=False,
                             mode='a')
            elif 'unc_sm33_X1' in file:
                df['norm']=df[17].multiply(unc3_x1)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None,index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                              index=False,
                             mode='a')

            elif 'unc_sm31_Xins' in file:
                df['norm']=df[17].multiply(unc1_xins)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):                    
                    df.to_csv('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1]), 
                     sep='\t',  
                     header=None, 
                             index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,index=False,
                             mode='a')
            elif 'unc_sm32_Xins' in file:
                df['norm']=df[17].multiply(unc2_xins)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None, 
                             index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                             mode='a', index=False)
                    
                    
            elif 'unc_sm33_Xins' in file:
                df['norm']=df[17].multiply(unc3_xins)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):                    
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None, index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                             mode='a', index=False)

            elif 'mapped_sm31_X1' in file:
                df['norm']=df[17].multiply(sm31_x1)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):                    
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None, index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,
                             mode='a', index=False)                    
                    
            elif 'mapped_sm32_X1' in file:
                df['norm']=df[17].multiply(sm32_x1)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):  
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None, index=False)  
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None, index=False,
                             mode='a')
            elif 'mapped_sm33_X1' in file:
                df['norm']=df[17].multiply(sm33_x1)
                
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):                  
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None, index=False)  
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None, index=False,
                             mode='a')                


            elif 'mapped_sm31_Xins' in file:
                df['norm']=df[17].multiply(sm31_xins)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):   
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  index=False,
                         header=None) 
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None,index=False,
                             mode='a') 
            elif 'mapped_sm32_Xins' in file:
                df['norm']=df[17].multiply(sm32_xins)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  index=False,
                         header=None)   
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  index=False,
                                 header=None,
                             mode='a') 
            elif 'mapped_sm33_Xins' in file:
                df['norm']=df[17].multiply(sm33_xins)
                if not os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                        path_nor, file.split(
                            'all_intersected_sorted_mapped_')[1]), 
                         sep='\t',  
                         header=None, index=False)
                if os.path.exists('{}/norm_{}.csv'.format(
                    path_nor, file.split(
                        'all_intersected_sorted_mapped_')[1])):
                    df.to_csv('{}/norm_{}.csv'.format(
                                path_nor, file.split(
                                    'all_intersected_sorted_mapped_')[1]), 
                                 sep='\t',  
                                 header=None, index=False,
                             mode='a') 

            
normalize_files()     


# In[ ]:


import sys,os
import csv
import pandas as pd
import numpy as np
from glob import glob

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def get_files_2(path_2, motif)-> list:
    """
    load all files names into a list
    """
    os.chdir(path_2)
    path_2=os.getcwd()
    files=glob(motif)
    return(files)

def output_folder(path_t):
    """
    creates new output folder 
    """
    if not os.path.exists(path_t):
        os.makedirs(path_t)
    os.chdir(path_t)
    path_norm=os.getcwd()
    return path_norm

col_13=['schMed4_repeatMaskerDNA_', 
        'schMed4_repeatMaskerLINE_', 'schMed4_repeatMaskerLTR_', 
       'schMed4_repeatMaskerUnknown_', 'maker']
col_14=['exon', 'five_prime_UTR', 'three_prime_UTR']
output_sense={} #dict
output_antisense={} #dict

for Te in col_13:
    if Te not in output_sense:
        output_sense[Te] = 0
    if Te not in output_antisense:
        output_antisense[Te] = 0
        
for kind in col_14:
    if kind not in output_sense:
        output_sense[kind] = 0
    if kind not in output_antisense:
        output_antisense[kind] = 0
        
def summary():
    """
    get normalized files, open with pandas 
    create df for sense and antisense reads
    find each annotated element in col 4 
    count the sum from weighted count 
    """
    #path='/data1/AnLi_FM_MA/data/output_sm3_KD/grep_miRNA/normalized/'
    path='/data1/AnLi_FM_MA/data/output_sm3_KD/'+
    '/normalized_total_miRNA'
    intersect=get_files_2(path, 'norm*.bed.csv')
    path_sum=output_folder(
        '/data1/AnLi_FM_MA/data/output_sm3_KD/grep_miRNA/summary_for_barplot')
    
    print(intersect)
    
    summary_string=['summary_sm3']

    with open(path_sum+'/maker_new_sm3_summary.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(summary_string)

    #i=intersect[0]
    #path='/data1/AnLi_FM_MA/data/output_sm3_KD/normalized/'
    
    for i in intersect: 
        print(i)
        for chunk in pd.read_csv(path+'/'+i,chunksize=10000, sep= '\t', 
                                 engine= 'python', header=None):
            
            print(chunk.head(2))

            sense=pd.DataFrame(chunk.loc[chunk[0]==chunk[2]]) #sense
            #print (sense.head(5))
            antisense=pd.DataFrame(chunk.loc[chunk[0]!=(chunk[2])]) #antisense
            #break
            #------------------------------------------------------------------
            maker_sense=pd.DataFrame(sense.loc[sense[3]=='maker'])
            maker_sense=maker_sense.reset_index()
            maker_antisense=pd.DataFrame(antisense.loc[antisense[3]=='maker'])
            maker_antisense=maker_antisense.reset_index()
            #-------------------------------------------------------------------
            for Te in col_13:
                for key, val in output_sense.items():
                        if Te == key:#output_sense.keys():
                            #print(Te)
                            te=pd.DataFrame(sense.loc[sense[3]==Te]) 
                            #df nut mit bestimmnten TE
                            sum_17=te[6].sum()
                            #print(sum_17)
                            output_sense[Te]+=(round(sum_17, 4)) 
                            #output_sense.value round(val, 8)
            #----------------------------------------------------------
                for key, val in output_antisense.items():
                        #print(key)
                        if Te == key:#output_sense.keys():
                            #print(Te)
                            te=pd.DataFrame(antisense.loc[antisense[3]==Te]) 
                            #df nut mit bestimmnten TE
                            sum_17=te[6].sum()
                            #print(sum_17)
                            output_antisense[Te]+=(round(sum_17, 4)) 
                            #output_sense.value
            #------------------------------------------------------------
            for kind in col_14:
                #print(kind)
                for key, value in output_sense.items():
                    if kind == key:
                        #print(kind)
                        Kind=pd.DataFrame(maker_sense.loc[
                            maker_sense[4].astype(str)==kind])
                        #print(Kind)
                        sum_17_14=Kind[6].sum() #weighted count norm [7]
                        output_sense[kind]+=(round(sum_17_14, 4)) 
                #-------------------------------------------------------
                for key, value in output_antisense.items():
                    if kind == key:
                        Kind=pd.DataFrame(maker_antisense.loc[(
                            maker_antisense[4].astype(str)==kind)])
                        sum_17_14=Kind[6].sum()
                        output_antisense[kind]+=(round(sum_17_14, 4))  
            #break
        with open(path_sum+'/maker_new_sm3_summary.csv', 'a') as csv_file:
                writer = csv.writer(csv_file)
                writer.writerow([i])
                writer.writerow(['sense'])
                for key, value in output_sense.items():
                    print(key, value)
                    writer.writerow([key, value])
                writer.writerow(['antisense'])
                for key, value in output_antisense.items():

                    writer.writerow([key, value])
        for key, value in output_sense.items():
            output_sense[key]=0
        for key, value in output_antisense.items():
            output_antisense[key]=0
        
        
summary()


# In[ ]:




