#!/usr/bin/env python
# coding: utf-8
#written by Anna Liznar in jupyter
# In[ ]:


import subprocess, os, sys
import pandas as pd 
import re
from glob import glob 
from zipfile import ZipFile


# In[ ]:


def total_number_reads():
    """
    go into different fastqc outputs to find the total number of sequenced reads
    1. go into Folder called 1-fastqc in CLASH-> Analysis-> 1-fastqc 
    2. there are 11 different folders as zip -> go into that and open folder 
        -> open fastqc_data.txt 
    3. go into txt file and find in lines Filename with string and Total Sequences with regrex  
        Filename 
        Total Sequences 0123456789
    4. write into dictionary
        dict itself key filename and value the other sub-dict
        dict in dict (sub-dict)-> total sequences and the amount of it as value
    """
    os.chdir('/data1/eCLASH/')
    path=os.getcwd()
    folders=glob('*.zip')
    dic_table={}
    
    for folder in folders:
        #unzip
        #print(folder)
        with ZipFile(folder, 'r') as zip:
            #print(zip)
            #zip.printdir()
            files=zip.namelist()
            for i in files:
                if i.endswith('data.txt'):
                    #zip.extract(i)
                    with zip.open(i) as f:
                        for line in f:
                            #print(line)
                            if 'Total Sequences' in str(line):
                                #print(line)
                                motif=re.findall(r'[0-9]{6,10}', str(line))
                                #print(motif)
                                dic_table[zip.filename.split(
                                    '_fastqc.zip')[0]]={
                                    'Total number of sequenced reads ':motif[0]}
                                
            

    return dic_table

total_number_reads()


# In[ ]:


def show_amount():
    """
    """
    dic_=total_number_reads()
    df=pd.DataFrame.from_dict(dic_)
    df=df.apply(pd.to_numeric)
    ax = df.plot.bar( rot=0)
    ax
    
    #plt.savefig('/data1/eCLASH/Stacked_total amount.pdf', bbox_inches='tight')
    return df
    
show_amount()    


# In[ ]:


def trimmed_reads():
    """
    go into different fastqc outputs to find the total number of sequenced reads
    1. go into Folder called 2-trimmomatic in 
        CLASH-> Analysis-> 2-trimmomatic -> report
    2. there are 11 different folders as zip -> go into that and open folder 
        -> open fastqc_data.txt 
    3. go into txt file and find in lines Filename with string and Total Sequences with regrex  
        Filename blablabla
        Total Sequences 0123456789
    4. update into dictionary
        dict itself key filename and value the other sub-dict
        dict in dict (sub-dict)-> total sequences and the amount of it as value
    """
    os.chdir('/data1/eCLASH/report/')
    path=os.getcwd()
    folders=glob('*.zip')
    dic_table=total_number_reads()
    
    for folder in folders:
        #unzip
        #print(folder)
        os.chdir(path)
        with ZipFile(folder, 'r') as zip:
            #print(zip)
            files=zip.namelist() #lists all items
            for i in files:
                if i.endswith('data.txt'):
                    #print(i) -> path to file and file 
                    with zip.open(i) as f:
                        for line in f:
                            if 'Total Sequences' in str(line):
                                #print(line)
                                motif=re.findall(r'[0-9]{6,10}', str(line))
                                #print(motif)
                                kkf=zip.filename.split(
                                    '_trimmomatic_fastqc.zip')[0]
                                dic_table[kkf.split('processed_')[1]].update(
                                    {'Total number of trimmed reads':motif[0]})
                                #update nested dictionary 
        
    return dic_table
trimmed_reads()            


# In[ ]:


def collapsed_reads():
    """
    to update the ditionary with collapsed reads 'collapsed_reads' opens 
    statistic output and finds name and number of sequences
    
    regular expressions are needed
    """
    os.chdir('/data1/eCLASH')
    path=os.getcwd()
    file=glob('col*.txt')[0]
    dic_table=trimmed_reads()
    
    #motif_filename=r'Nr[0,1]{1}[8,9,0,1,2]{1}_[a-z]{4,5}_[1-3]{1}'
    motif_seq=r'[0-9]{7,8}'
    
    os.chdir('/data1/eCLASH')
    with open(file, 'r') as f:
        for line in f:
            #print(line)
            namee=line.split('collapser_processed_')[1]
            name=namee.split('_trimmomatic')[0]
            #name=re.findall(motif_filename, line)
            #print(name)
            no=re.findall(motif_seq, line)
            #print(name[0], no[0])
            
            dic_table[name].update({
                'Total number of reads after removing PCR duplicates':no[0]})
            
    return dic_table
collapsed_reads()    


# In[ ]:


def sortme_reads():
    """
    to update the dictionary with reads after rRNA extraction 'sortme_reads' opens 
    statistic output and finds name and number of sequences
    
    regular expressions are needed
    """
    os.chdir('/data1/eCLASH')
    path=os.getcwd()
    file=glob('sort*.txt')[0]
    dic_table=collapsed_reads()
    
    motif_filename=r'[0-9]{1,2}_[0-9]{1}[a-z]{1}'
    motif_seq=r'[0-9]{2,9}'
    
    os.chdir('/data1/eCLASH')
    with open(file, 'r') as f:
        for line in f:

            if 'Total reads = ' in line:
                #print(line)
                namee=line.split('processed_')[1]
                name=namee.split('_trimmomatic_rRNA')[0]
                #name=re.findall(motif_filename, line)
                
            elif 'Total reads passing E-value threshold' in line:
                no=re.findall(motif_seq, line)
                #print(name[0], no[0])
                dic_table[name].update({'Total number of rRNA':no[0]})
                
            elif 'Total reads failing E-value threshold' in line:
                noother=re.findall(motif_seq, line)
                dic_table[name].update({'Total number of reads other than rRNA':noother[0]})
            
            
    return dic_table
sortme_reads()   


# In[ ]:


def STAR_reads():
    """
    """
    os.chdir('/data1/eCLASH/')
    path=os.getcwd()
    file=glob('STAR*.txt')[0]
    #dic_table=sortme_reads()
    dic_STAR={}
    
    #motif_filename=r'[0-9]{1,2}_[0-9]{1}[a-z]{1}'
    motif_seq=r'[0-9]{3,8}'
    
    os.chdir('/data1/eCLASH/')
    with open(file, 'r') as f:
        for line in f:
            #print(line)
            #processed_Nr08_ctrl_1_.Log.final.out
            #namee=line.split('processed_')[1]
            name=line.split('_.Log.final.out')[0]
            #name=re.findall(motif_filename, line)
            no=re.findall(motif_seq, line)
            #print(no)
            word=re.findall(r'[a-z]{7}_[a-z]{6,10}_[a-z]{4,6}', line) #20 #27 #29 #24
            word_0=str(word[0])
            
            if name not in dic_STAR.keys():
                dic_STAR[name]={word_0:no[0]}
            else:
                dic_STAR[name].update({word_0:no[0]})
            
            
    return dic_STAR#, dic_table      
    #return dic_table
STAR_reads()   


# In[ ]:


#make adataframe out of dictionaries 
def dic_to_df1():
    """
    save dictionary from above functions to a dataframe 
    """
    dic_STAR=STAR_reads()
    
    reorder=['Total_Amount', 'genomic_mapper_reads', 
             '% genomic mapper reads',
             'genomic_unique_mapper',
             '% unique mapper reads',
         'genomic_multiple_mapper','% multimapper reads',
             'genomic_unmappable_reads', '% unmappable reads']

    col_reorder=['Nr08_ctrl_1','Nr09_ctrl_2','Nr10_piRNA_1','Nr11_piRNA_2','Nr12_piRNA_3']
    df_STAR=pd.DataFrame(dic_STAR)
    df_STAR=df_STAR.apply(pd.to_numeric)

    df_STAR.loc['Total_Amount']=df_STAR.loc[
        'genomic_mapper_reads']+df_STAR.loc[
        'genomic_unmappable_reads']
        
    df_STAR.loc['% genomic mapper reads']=df_STAR.loc[
        'genomic_mapper_reads']/df_STAR.loc[
        'Total_Amount']*100
    df_STAR.loc['% unique mapper reads']=df_STAR.loc[
        'genomic_unique_mapper']/df_STAR.loc[
        'Total_Amount']*100
    
    df_STAR.loc['% multimapper reads']=df_STAR.loc[
        'genomic_multiple_mapper']/df_STAR.loc['Total_Amount']*100
    
    df_STAR.loc['% unmappable reads']=df_STAR.loc[
        'genomic_unmappable_reads']/df_STAR.loc['Total_Amount']*100

    pd.options.display.float_format = '{:.2f}'.format
    df_STAR=df_STAR.reindex(reorder) #reorder index
    df_STAR=df_STAR[col_reorder] #change order columns

    
    #df_table.to_csv('/data1/eCLASH/preliminary_analysis.csv', sep='\t', 
   # header=True, index=True)
    
    #df_STAR.to_csv('/data1/eCLASH/1_UMI_Overview.csv', sep='\t', header=True, index=True)
    
    return df_STAR

dic_to_df1()


# In[ ]:


def calc_chimeric():
    """
    take pairs.aligned.tab and calc amount
    """
    path='/data1/eCLASH/'
    results_path='4-Chimera_Alex/results_not_on_genes/'
    reps=['ctrl_1', 'ctrl_2', 'piRNA_1', 'piRNA_2', 'piRNA_3']
    path_list=[]
    for rep in reps:
        path_list.append(results_path+rep+'/pairs.aligned.tab')
    for i in path_list:
        cd='cat {} | wc -l > {}_chimera.txt'.format(path+i, path+i.split('/pairs.aligned.tab')[0])
        print(cd)
        subprocess.check_call(cd, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
calc_chimeric()    


# In[ ]:


def make_chimeric_df():
    """
    read fragments.large.filtered.Aligned.bam and count chimeric reads (from all reps)
    make df
    """
    df_path='/data1/eCLASH/'
    p='/data1/eCLASH/1-Zamore_UMI/NEW/2-Cutadapt/results/'
    results_path=df_path+'4-Chimera_Alex/results_on_genes/'

    reps=['Nr08_ctrl_1_cut', 'Nr09_ctrl_2_cut', 'Nr10_piRNA_1_cut', 'Nr11_piRNA_2_cut', 'Nr12_piRNA_3_cut']
    path_list=[]
    for rep in reps:
        path_list.append(p+rep+'/fragments.large.filtered.Aligned.bam')
    #print(path_list)
    for i in path_list:
        cd="samtools view '{}' | cut -f1 | sed 's/-.*//' - | sort | uniq | wc -l > {}_possiblechimeras.txt".format(
        i, i.split('/fragments.large.filtered.Aligned.bam')[0])
        print(cd)
        subprocess.check_call(cd, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        #break
        
make_chimeric_df()


# In[ ]:


def append_chimeras():
    """
    open Overview_mapped.csv 
    make two df upper band(u) and lower band(d)
    include chimeric reads in Overview of u
    save as csv
    """
    df_path='/data1/eCLASH/'
    txt_file_path=df_path+'4-Chimera_Alex/results_not_on_genes/'
    
    names=['Nr08_ctrl_1', 'Nr09_ctrl_2', 'Nr10_piRNA_1', 'Nr11_piRNA_2', 'Nr12_piRNA_3']
    filename='_chimera.txt'
    dic_chimera={}
    os.chdir(txt_file_path)
    for name in names:
        with open(txt_file_path+name+filename, 'r') as f:
            line=f.readline().rstrip() #rstrip() looses the newline character
            dic_chimera[name]={'chimeric read':line}

    df_chim=pd.DataFrame.from_dict(dic_chimera)
    print(df_chim)
    df_chim=df_chim.apply(pd.to_numeric)
    df_chim['Nr10_piRNA_1']=(df_chim['Nr10_piRNA_1'])
    df_chim['Nr12_piRNA_3']=(df_chim['Nr12_piRNA_3'])
    df_chim['Nr11_piRNA_2']=(df_chim['Nr11_piRNA_2'])
    df=pd.read_csv('/data1/eCLASH/1_UMI_Overview_FÜRDIEMA_TOTAL_AMOUNT.csv', sep='\t', index_col= [0])

    df_upper=df[['Nr08_ctrl_1', 'Nr09_ctrl_2', 'Nr10_piRNA_1', 'Nr11_piRNA_2', 'Nr12_piRNA_3']]

    df_upper=df_upper.append(df_chim)

    df_upper.loc['genomic_unmappable_reads']=df_upper.loc[
        'genomic_unmappable_reads']-df_upper.loc['chimeric read']
    
    df_upper.loc['% chimeric read']=df_upper.loc[
                                            'chimeric read']*100/df_upper.loc[
                                                                  'Total_Amount']
    df_upper.loc['% unmappable reads']=df_upper.loc['% unmappable reads']-df_upper.loc['% chimeric read']
    
    #df_upper.to_csv('/data1/eCLASH/Overview_chimeric.csv',
                   # sep='\t', header=True, index=True)
    
    
    
    return df_upper
    
append_chimeras()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


# In[ ]:


def make_stacked_bp():
    """
    from df make stacked bar plot 
    """
    #os.chdir('/data1/eCLASH/')
    #os.chdir('/data1/eCLASH/')
    p=os.getcwd()

    #file=glob('1_UMI_Overview_mapped_TOTAL_AMOUNT.csv')
    file=glob('Overview_chimeric.csv')

    df=pd.read_csv(p+'/'+file[0], sep='\t', index_col=[0])

    #df=df.drop(['Total_Amount','% unique mapper reads', '% multimapper reads', 'chimeric read',
         #       '% chimeric read', '% unmappable reads'])
    df=df.drop(['Total_Amount','% unique mapper reads', '% multimapper reads', '% chimeric read',
                 '% unmappable reads', '% genomic mapper reads'])
    #df=df.rename(index=str, columns={'1_1u':'1', '3_2u':'2', 
                                    # '5_3u':'3', '7_4u':'4','9_5u':'5'})
    df=df/1000

    df_t=df.transpose()
    df_t=df_t.rename(index=str, columns={'genomic_mapper_reads':'Mapper reads', 
                                         'genomic_unique_mapper_reads':'Unique mapper reads',
                                        'genomic_multiple_mapper_reads': 'Multimapper reads',
                                        'genomic_unmappable_reads':'Unmappable reads',
                                        'chimeric reads': 'Chimeric reads'})

    #fig = stacked_barplot(df, rotation=45, legend_loc='best')
    #fig.show()
    pal = sns.color_palette("Spectral")
    pal=['firebrick', 'lightcoral', 'green', 'navy','mediumslateblue']
    #pal=sns.color_palette('Dark2_r')
    p=df_t.plot.bar(stacked=True,  color=pal, alpha=0.4, rot=35, title='Mapping events')
    p.set_ylabel("Amount of mappings *10⁴", fontsize=12)
    p.set_xlabel("Replicate", fontsize=12)
    axes = plt.gca()
    axes.set_ylim([0,180])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('./UMI_Stacked_bp_chim_mappingevents.pdf', 
                #bbox_inches='tight')
    plt.show()
    plt.close
make_stacked_bp()


# In[ ]:


def make_stacked_bp_perc():
    """
    from df make stacked bar plot 
    """
    os.chdir('/data1/eCLASH/')#4-Chimera_Alex/results_not_on_genes')
    p=os.getcwd()
    #print(p)

    #file=glob('1_UMI_Overview_mapped_TOTAL_AMOUNT.csv')
    file=glob('Overview_chimeric.csv')
    #print(file)
    df=pd.read_csv(p+'/'+file[0], sep='\t', index_col=[0])
    #print(df)
    #Nr08_ctrl_1','Nr09_ctrl_2','Nr10_piRNA_1','Nr11_piRNA_2','Nr12_piRNA_3
    df=df.drop(['Total_Amount','genomic_mapper_reads', 'genomic_unique_mapper',
                '% genomic mapper reads',
                'genomic_multiple_mapper', 'genomic_unmappable_reads', 'chimeric read'])
    df=df.rename(index=str, columns={'Nr08_ctrl_1':'ctrl_1', 
                                     'Nr09_ctrl_2':'ctrl_2', 
                                     'Nr10_piRNA_1':'eCLASH_1', 
                                     'Nr11_piRNA_2':'eCLASH_2',
                                     'Nr12_piRNA_3':'eCLASH_3'})
    #df=df/1000000

    df_t=df.transpose()
    #df_t=df_t.rename(index=str, columns={'genomic_mapper_reads':'Mapper reads', 
                                     #    'genomic_unique_mapper_reads':'Unique mapper reads',
                                     #   'genomic_multiple_mapper_reads': 'Multimapper reads',
                                     #   'genomic_unmappable_reads':'Unmappable reads',
                                      #  'chimeric reads': 'Chimeric reads'})
    print(df_t)
    df_aa=df_t.loc[['eCLASH_1', 'eCLASH_2', 'eCLASH_3']]
    df_a=pd.DataFrame(df_aa.agg('mean'))
   # print((df_a))
    df_a=df_a.rename(index=str, columns={0:'1-3'})
    df_a=df_a.transpose()
    print((df_a))
    
    #fig = stacked_barplot(df, rotation=45, legend_loc='best')
    #fig.show()
    pal = sns.color_palette("Dark2_r")
    p=df_t.plot.bar(stacked=True,  color=pal, alpha=0.4, rot=0, title='Mapping events')
    p.set_ylabel("Distribution of mapping events [%]", fontsize=12)
    p.set_xlabel("Replicate", fontsize=12)
    axes = plt.gca()
    axes.set_ylim([0,110])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('./UMI_Stacked.pdf', bbox_inches='tight')
    plt.show()
    plt.close
    
    
    #
    pal = sns.color_palette("Dark2_r")
    p=df_a.plot.bar(stacked=True,  color=pal, alpha=0.4, rot=0, title='Mapping events', width=0.3)
    p.set_ylabel("Distribution of mapping events [%]", fontsize=12)
    p.set_xlabel("Replicate", fontsize=12)
    axes = plt.gca()
    axes.set_ylim([0,110])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    fig = plt.gcf()
    fig.set_size_inches(1.5, 4)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    
    #plt.savefig('./Stacked.pdf', bbox_inches='tight')
    plt.show()
    plt.close
make_stacked_bp_perc()

