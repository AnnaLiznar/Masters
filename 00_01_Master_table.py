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
        Filename blablabla
        Total Sequences 0123456789
    4. write into dictionary
        dict itself key filename and value the other sub-dict
        dict in dict (sub-dict)-> total sequences and the amount of it as value
    """
    os.chdir('/data1/Master/')
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
    os.chdir('/data1/Master/CLASH/Analysis/2-trimmomatic/report/')
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
                                dic_table[zip.filename.split(
                                    '_trimmomatic_fastqc.zip')[0]].update(
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
    os.chdir('/data1/Master/CLASH/Analysis/')
    path=os.getcwd()
    file=glob('col*.txt')[0]
    dic_table=trimmed_reads()
    
    motif_filename=r'[0-9]{1,2}_[0-9]{1}[a-z]{1}'
    motif_seq=r'[0-9]{7,8}'
    
    os.chdir('/data1/Master/CLASH/Analysis/')
    with open(file, 'r') as f:
        for line in f:
            name=re.findall(motif_filename, line)
            no=re.findall(motif_seq, line)
            
            dic_table[name[0]].update({
                'Total number of reads after removing PCR duplicates':no[0]})
            
    return dic_table
collapsed_reads()    


# In[ ]:


def sortme_reads():
    """
    to update the ditionary with reads after rRNA extraction 'sortme_reads' opens 
    statistic output and finds name and number of sequences
    
    regular expressions are needed
    """
    os.chdir('/data1/Master/')
    path=os.getcwd()
    file=glob('sort*.txt')[0]
    dic_table=collapsed_reads()
    
    motif_filename=r'[0-9]{1,2}_[0-9]{1}[a-z]{1}'
    motif_seq=r'[0-9]{5,8}'
    
    os.chdir('/data1/Master/CLASH/Analysis/')
    with open(file, 'r') as f:
        for line in f:
            #print(line)
            if 'Total reads = ' in line:
                #print(line)
                name=re.findall(motif_filename, line)
            elif 'Total reads passing E-value threshold' in line:
                no=re.findall(motif_seq, line)
                #print(name[0], no[0])

                dic_table[name[0]].update({'Total number of rRNA':no[0]})
            elif 'Total reads failing E-value threshold' in line:
                noother=re.findall(motif_seq, line)
                dic_table[name[0]].update({'Total number of reads other than rRNA':noother[0]})
            
            
    return dic_table
sortme_reads()   


# In[ ]:


def STAR_reads():
    """
    to update the ditionary with reads after rRNA extraction 'sortme_reads' opens 
    statistic output and finds name and number of sequences
    
    regular expressions are needed
    """
    os.chdir('/data1/Master/CLASH/')
    path=os.getcwd()
    file=glob('STAR*.txt')[0]
    #dic_table=sortme_reads()
    dic_STAR={}
    
    motif_filename=r'[0-9]{1,2}_[0-9]{1}[a-z]{1}'
    motif_seq=r'[0-9]{6,8}'
    
    os.chdir('/data1/Master/CLASH/')
    with open(file, 'r') as f:
        for line in f:
            #print(line)
            name=re.findall(motif_filename, line)
            no=re.findall(motif_seq, line)
            word=re.findall(r'\w{20,29}', line) #20 #27 #29 #24
            word_0=str(word[0])
            #print(word_0)
            #print(name, no)
            
            #dic_table[name[0]].update({word_0:no[0]})
            
            if name[0] not in dic_STAR.keys():
                dic_STAR[name[0]]={word_0:no[0]}
            else:
                dic_STAR[name[0]].update({word_0:no[0]})
            
            
    return dic_STAR#, dic_table      
    #return dic_table
STAR_reads()   


# In[ ]:


#make adataframe out of dictionaries 
def dic_to_df():
    """
    save dictionary from above functions to a dataframe 
    """
    dic_STAR=STAR_reads()
    dic_table=sortme_reads()
    
    reorder_index=['Total number of sequenced reads ', 'Total number of trimmed reads', 
                   'Total number of reads after removing PCR duplicates', 
                   'Total number of rRNA', '% rRNA',
                   'Total number of reads other than rRNA', '% reads other than rRNA']
    
    reorder=['Total_Amount', 'genomic_mapper_reads', 'genomic_unique_mapper_reads',
             '% unique mapper reads',
         'genomic_multiple_mapper_reads','% multimapper reads',
             'genomic_unmappable_reads', '% unmappable reads']
    
    df=pd.DataFrame(dic_table)
    df=df.apply(pd.to_numeric)
    #df=df.astype(int)
    
    df.loc['% rRNA']=df.loc['Total number of rRNA']/df.loc[
        'Total number of reads after removing PCR duplicates']*100
    df.loc['% reads other than rRNA']=df.loc['Total number of reads other than rRNA']/df.loc[
        'Total number of reads after removing PCR duplicates']*100
    #print(df)
    pd.options.display.float_format = '{:.2f}'.format
    df_table=df.reindex(reorder_index)
    columns_table=list(df_table.columns.values)
    #print(columns_table)
    
    col_reorder=['1_1u','2_1d','3_2u','4_2d','5_3u','6_3d','7_4u','8_4d','9_5u', '10_5d','11_8c']
    df_table=df_table[col_reorder] #change order of header/columns
    
    df_STAR=pd.DataFrame(dic_STAR)
    df_STAR=df_STAR.apply(pd.to_numeric)
    df_STAR.loc['Total_Amount']=df_STAR.loc[
        'genomic_mapper_reads']+df_STAR.loc[
        'genomic_unmappable_reads']
    
    df_STAR.loc['% unique mapper reads']=df_STAR.loc[
        'genomic_unique_mapper_reads']/df_STAR.loc[
        'Total_Amount']*100
    
    df_STAR.loc['% multimapper reads']=df_STAR.loc[
        'genomic_multiple_mapper_reads']/df_STAR.loc['Total_Amount']*100
    
    df_STAR.loc['% unmappable reads']=df_STAR.loc[
        'genomic_unmappable_reads']/df_STAR.loc['Total_Amount']*100
    print(df_STAR)
    pd.options.display.float_format = '{:.2f}'.format
    df_STAR=df_STAR.reindex(reorder) #reorder index
    df_STAR=df_STAR[col_reorder] #change order columns
    
    #print(os.getcwd())
    
    #df_table.to_csv('./Overview.csv', sep='\t', header=True, index=True)
    
    df_STAR.to_csv('./Overview.csv', sep='\t', header=True, index=True)
    
    return df_table, df_STAR

dic_to_df()


# In[ ]:


def calc_chimeric():
    """
    take pairs.aligned.tab and calc amount
    """
    path='/data1/Master/CLASH/Analysis'
    results_path='/results/'
    reps=['1_1u', '3_2u', '5_3u', '7_4u', '9_5u']
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
    df_path='/data1/Master/CLASH/Analysis'
    results_path=df_path+'/results/'
    reps=['1_1u', '3_2u', '5_3u', '7_4u', '9_5u']
    path_list=[]
    for rep in reps:
        path_list.append(results_path+rep+'/fragments.large.filtered.Aligned.bam')
    #print(path_list)
    for i in path_list:
        cd="samtools view {} | cut -f1 | sed 's/-.*//' - | sort | uniq | wc -l > {}_chimera.txt".format(
        i, i.split('/fragments.large.filtered.Aligned.bam')[0])
        print(cd)
        subprocess.check_call(cd, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              shell=True)
        
make_chimeric_df()


# In[ ]:


def append_chimeras():
    """
    open Overview_mapped.csv 
    make two df upper band(u) and lower band(d)
    include chimeric reads in Overview of u
    save as csv
    """
    Overview_mapped_path='/data1/Master/CLASH/Analysis/'
    txt_file_path=Overview_mapped_path+'/results/'
    names=['1_1u', '3_2u', '5_3u', '7_4u', '9_5u']
    filename='_chimera.txt'
    dic_chimera={}
    os.chdir(txt_file_path)
    for name in names:
        with open(txt_file_path+name+filename, 'r') as f:
            line=f.readline().rstrip() #rstrip() looses the newline character

            dic_chimera[name]={'chimeric read':line}

    df_chim=pd.DataFrame.from_dict(dic_chimera)

    df=pd.read_csv(Overview_mapped_path+'Overview_mapped_TOTAL_AMOUNT.csv', sep='\t', index_col= [0])

    df_upper=df[['1_1u', '3_2u', '5_3u', '7_4u', '9_5u']]

    df_upper=df_upper.append(df_chim)

    df_upper=df_upper.apply(pd.to_numeric)

    df_upper.loc['genomic_unmappable_reads']=df_upper.loc[
        'genomic_unmappable_reads']-df_upper.loc['chimeric read']
    #print(df_upper)
    
    df_upper.loc['% chimeric read']=df_upper.loc[
                                            'chimeric read']*100/df_upper.loc[
                                                                  'Total_Amount']
    df_upper.loc['% unmappable reads']=df_upper.loc['% unmappable reads']-df_upper.loc['% chimeric read']
    #print(df_upper)
    
    #df_upper.to_csv('/data1/Master/CLASH/Analysis/Overview_upper_chimeric.csv',
                   # sep='\t', header=True, index=True)
    
    
    
    return df_upper
    
append_chimeras()


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def make_stacked_bp():
    """
    from df make stacked bar plot 
    """
    os.chdir('/data1/Master/CLASH/Analysis/')
    p=os.getcwd()
    #print(p)
    file=glob('Overview.csv')
    print(file)
    df=pd.read_csv(p+'/'+file[0], sep='\t', index_col=[0])
    print(df)
    df=df.drop(['Total_Amount','% unique mapper reads', '% multimapper reads', 
                '% chimeric read', '% unmappable reads'])
    df=df.rename(index=str, columns={'1_1u':'1', '3_2u':'2', 
                                     '5_3u':'3', '7_4u':'4','9_5u':'5'})
    df=df/1000000
    print(df)
    df_t=df.transpose()
    df_t=df_t.rename(index=str, columns={'genomic_mapper_reads':'Mapper reads', 
                                         'genomic_unique_mapper_reads':'Unique mapper reads',
                                        'genomic_multiple_mapper_reads': 'Multimapper reads',
                                        'genomic_unmappable_reads':'Unmappable reads',
                                        'chimeric reads': 'Chimeric reads'})

    pal = sns.color_palette("tab20")
    p=df_t.plot.bar(stacked=True,  color=pal, alpha=0.4, rot=0, title='Mapping events')
    p.set_ylabel("Amount of mappings *10⁶", fontsize=12)
    p.set_xlabel("Replicate", fontsize=12)
    axes = plt.gca()
    axes.set_ylim([0,11.5])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('./Stacked.pdf', bbox_inches='tight')
    plt.show()
    plt.close
make_stacked_bp()


# In[ ]:


def make_stacked_bp_perc():
    """
    from df make stacked bar plot 
    """
    os.chdir('/data1/Master/CLASH/Analysis/')
    p=os.getcwd()
    #print(p)
    file=glob('Overview_upper_chimeric.csv')

    df=pd.read_csv(p+'/'+file[0], sep='\t', index_col=[0])

    
    df=df.drop(['Total_Amount','genomic_mapper_reads', 'genomic_unique_mapper_reads',
                'genomic_multiple_mapper_reads', 'genomic_unmappable_reads', 'chimeric read'])
    df=df.rename(index=str, columns={'1_1u':'1', '3_2u':'2', '5_3u':'3', '7_4u':'4','9_5u':'5'})
    #df=df/1000000

    df_t=df.transpose()
    #df_t=df_t.rename(index=str, columns={'genomic_mapper_reads':'Mapper reads', 
                                     #    'genomic_unique_mapper_reads':'Unique mapper reads',
                                     #   'genomic_multiple_mapper_reads': 'Multimapper reads',
                                     #   'genomic_unmappable_reads':'Unmappable reads',
                                      #  'chimeric reads': 'Chimeric reads'})
    print(df_t)
    df_aa=df_t.loc[['2', '3', '4']]
    df_a=pd.DataFrame(df_aa.agg('mean'))

    df_a=df_a.rename(index=str, columns={0:'1-5'})
    df_a=df_a.transpose()
    print((df_a))

    pal = sns.color_palette("Dark2_r")
    p=df_t.plot.bar(stacked=True,  color=pal, alpha=0.4, rot=0, title='Mapping events')
    p.set_ylabel("Distribution of mapping events [%]", fontsize=12)
    p.set_xlabel("Replicate", fontsize=12)
    axes = plt.gca()
    axes.set_ylim([0,110])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    #plt.savefig('./Stacked.pdf', bbox_inches='tight')
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
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)

    fig = plt.gcf()
    fig.set_size_inches(1.5, 4)
    #plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), fontsize=12)
    
    #plt.savefig('./Stacked.pdf', bbox_inches='tight')
    plt.show()
    plt.close
make_stacked_bp_perc()


# In[ ]:




