#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys,os, subprocess
import csv
import pandas as pd
import numpy as np
from glob import glob
import re

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[ ]:


def get_files_2(path_2, motif)-> list:
    """
    load all files names into a list
    """
    os.chdir(path_2)
    path_2=os.getcwd()
    files=glob(motif)
    return(files)

#get_files_2('/data1/AnLi_FM_MA/data/output_sm3_KD/normalized/', 'norm*.csv')


# In[ ]:


def output_folder(path_t):
    """
    creates new output folder 
    """
    if not os.path.exists(path_t):
        os.makedirs(path_t)
    os.chdir(path_t)
    path_outp=os.getcwd()
    return path_outp


# In[ ]:


def get_csv(path):
    """
    open csv file 
    """
    df=pd.read_csv(path, sep=',', index_col=0, engine='python')
    return df


# In[ ]:


def find_elements():
    """
    open files with shell commands 
    input: files -> intersected files with exons, DNA, LTR, LINE and Unkn Tes and introns
    output: txt files with amount of chimeric sequences 
    """
    bedtools_path='/data1/eCLASH/6-bedtools/'
    dirs=['Cluster', 'DNA', 'Exon', 'Intron', 'LINE', 'LTR', 'No_annot', 'Unk']
    motif='grep*.bed'
    output_file_path='/data1/eCLASH/6-bedtools/Overview_targets'
    path_outp=output_folder(output_file_path) #creates new folder 
    
    for d in dirs:
        files=get_files_2(bedtools_path+d, motif)
        print(files)
        for file in files:
            ff=file.split('_fragments.large.filtered.Aligned.bed')[0]
            #print('this is the path: {}, this is the file: {}'.format(bedtools_path+d, file))
            if d=='No_annot':
                cmd="cat {} | cut -f4 | sed 's/-.*//' - | sort -u | wc -l > {}.txt".format(
                bedtools_path+d+'/'+file, path_outp+'/'+ff.split('grep_chimera_')[1])
                print(cmd)
                subprocess.check_call(cmd, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, 
                                 shell=True)
                #if no annot dont compare $6 and $12 
            else:
                cmd="cat {} | awk '($6==$12) {}' - | cut -f4 | sed 's/-.*//' - | sort -u | wc -l > {}.txt".format(
                    bedtools_path+d+'/'+file, 
                    '{print $0}',
                    path_outp+'/'+ff.split('grep_chimera_')[1])
                print(cmd)
                subprocess.check_call(cmd, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, 
                                 shell=True)
            
        
find_elements()


# In[ ]:


def get_total_chimeric_reads():
    """
    open with samtools view the bam file (large filtered aligned) includes all 'large' fragments of 
    chimeric reads
    and grep the ones in pairs.aligned.tab to make sure real chimeras are only counted
    count those 
    """
    names=['ctrl_1', 'ctrl_2', 'piRNA_1', 'piRNA_2', 'piRNA_3']
    general_path='/data1/eCLASH/4-Chimera_Alex/results_not_on_genes/amount_chimeras/'
    filename='fragments.large.filtered.Aligned.bam'
    
    chim_dic={}
    
    os.chdir(general_path)
    
    for name in names:
        file=general_path+name+'/'+filename
        txt=general_path+name+'/'+name+'_chimeric_names.txt'
        #print(file)
        #print(txt)
        
        c="samtools view {} | grep -Fwf {} - | cut -f1 | sort -u | wc -l > {}.txt".format(
            file, 
            txt,
            general_path+name+'_amount_chimeras')
        print(c)
        subprocess.check_call(c, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.PIPE, 
                                 shell=True)
        #break
        
get_total_chimeric_reads()
#samtools view 'bam' | grep -Fwf '//results_all/1_1u/1_1u_input_names.txt' - |cut -f1 | sort -u | wc -l 
#/data1/Master/CLASH/Analysis/6-chimeraSearch_clipseqtool/results_all/1_1u/1_1u_input_names.txt


# In[ ]:


#get chiemeras 

def append_chimeras():
    """
    open Overview_mapped.csv 
    make two df upper band(u) and lower band(d)
    include chimeric reads in Overview of u
    save as csv
    """
    Overview_mapped_path='/data1/eCLASH/'
    txt_file_path=Overview_mapped_path+'4-Chimera_Alex/results_not_on_genes/amount_chimeras/'
    names=['Nr08_ctrl_1', 'Nr09_ctrl_2', 'Nr10_piRNA_1', 'Nr11_piRNA_2', 'Nr12_piRNA_3']
    filename='_chimera.txt'
    dic_chimera={}
    os.chdir(txt_file_path)
    for name in names:
        with open(txt_file_path+'/'+name+filename, 'r') as f:
            line=f.readline().rstrip() #rstrip() looses the newline character
            #print(line)
            dic_chimera[name]={'chimeric read':line}
    #print(dic_chimera)
    df_chim=pd.DataFrame.from_dict(dic_chimera)
    print(df_chim)
    overview_path='/data1/eCLASH/1_UMI_Overview_mapped_TOTAL_AMOUNT.csv'
    df=pd.read_csv(overview_path, sep='\t', index_col=0)
    #print(df)

    #print(df_upper)
    df=df.append(df_chim)

    df=df.apply(pd.to_numeric)
    print(df)

    df.loc['% chimeric read']=df.loc[
                                            'chimeric read']*100/df.loc[
                                                                  'genomic_unmappable_reads']
    
    return df
    
append_chimeras()


# In[ ]:


def append_to_overview():
    """
    open txt files in 7-bedtools/Overview 
    append amount to Overview 
    """
    df=append_chimeras()
    df=df.rename(index=str, columns={'Nr08_ctrl_1':'ctrl_1', 
                                         'Nr09_ctrl_2':'ctrl_2',
                                        'Nr10_piRNA_1': 'piRNA_1',
                                        'Nr11_piRNA_2':'piRNA_2',
                                        'Nr12_piRNA_3': 'piRNA_3'})
    path='/data1/eCLASH/6-bedtools/'
    files=get_files_2(path+'Overview_targets', '*txt')
    
    ########## get df 

    dirs=['Cluster', 'DNA', 'Exon', 'Intron', 'LINE', 'LTR', 'No_annot', 'Unk']
    #print(df)
    #open txt files 
    dirs=sorted(dirs)
    files=sorted(files)
    dic_elemets={}

    for file in files:
            motif=re.findall(r'[a-z]{2,7}_[a-z]{5}_', file, re.IGNORECASE)
            #print(motif)
            motif=motif[0]
            mot=re.findall(r'[a-z]{2,7}', file, re.IGNORECASE)
            mo=mot[0]
            name=file.split(motif)[1]
            with open(path+'Overview_targets/'+file, 'r') as f:
                line=f.readline().rstrip() #rstrip() looses the newline character
                #print(line)
                if name.split('.txt')[0] not in dic_elemets.keys():
                    dic_elemets[name.split('.txt')[0]]={mo:line}
                else:
                    dic_elemets[name.split('.txt')[0]].update({mo:line})
    #print(dic_elemets)
    df_e=pd.DataFrame(dic_elemets)
    #print(df_e)
    df_e=df_e.rename(index= {'No': 'No annotation', 'Unk':'Unknown'}) #rename index
    
    df_all=df.append(df_e)
    
    df_all.to_csv(path+'/Overview_ALL_ELEMENTS.csv', header=True, index=True, sep='\t')
    
    return(df_all)
    
append_to_overview()


# In[ ]:


import sys, os 
import csv 
import pandas as pd 
from glob import glob 
import matplotlib.pyplot as plt
import numpy as np
import re

from pandas.tools.plotting import table


# In[ ]:


def get_df():
    """
    open csv file as df with pandas
    """
    path='/data1/eCLASH/6-bedtools/'
    os.chdir(path)
    p=os.getcwd()
    #print(p)
    file=glob('*ALL_E*')
    #print(file)
    df=pd.read_csv(path+'/'+file[0], sep='\t', index_col=0)
    df=df.iloc[[8,10,11,12,13, 14, 15, 16, 17],:]
    #print(df.index.tolist())
    return(df)
get_df()


# In[ ]:


def perc():
    """
    count percentages for plots
    """
    df=get_df()
    df.loc['% Cluster']=df.loc['Cluster']*100/df.loc['chimeric read']
    df.loc['% DNA']=df.loc['DNA']*100/df.loc['chimeric read']
    df.loc['% Exon']=df.loc['Exon']*100/df.loc['chimeric read']
    df.loc['% Intron']=df.loc['Intron']*100/df.loc['chimeric read']
    df.loc['% LINE']=df.loc['LINE']*100/df.loc['chimeric read']
    df.loc['% LTR']=df.loc['LTR']*100/df.loc['chimeric read']
    df.loc['% No annotation']=df.loc['No annotation']*100/df.loc['chimeric read']
    df.loc['% Unknown']=df.loc['Unknown']*100/df.loc['chimeric read']
    pd.options.display.float_format = '{:.2f}'.format
    #print(df)
    return (df)
perc()


# In[ ]:


def make_piechart():
    """
    take df from get_df
    plot a pie chart
    """
    #df=get_df()
    #df=df.iloc[[1,2,3,4,5,6,7,8],:]
    df=perc()
    df=df.iloc[[9,10,11,12,13,14,15,16],2:5]
    #pd.options.display.float_format = '{:.2f}'.format
    
    #print(cols)
    #plot=df.plot.pie(df)
    df.index=df.index.str.replace('%','', regex=True)
    df['all replicates']=df[['piRNA_1','piRNA_2','piRNA_3']].mean(axis=1)
    cols=df.columns.tolist()
    print(df)
    for rep in cols:
        
    
        plt.figure(figsize=(16,6))

        ax1=plt.subplot(121, aspect='equal')
        df.plot(kind='pie', 
               y=rep, ax=ax1,
               startangle=90,
               shadow=False,
               legend=False,
               fontsize=12,
               label=rep)



        # plot table
        ax2 = plt.subplot(122)
        plt.axis('off')
        #tbl = table(ax2, df, loc='center')
        #tbl.auto_set_font_size(False)
        #tbl.set_fontsize(14)
        #plt.show()
        
        values = df[rep]
        labels = ax2
        explode = (0, 0, 0.1, 0, 0, 0, 0, 0)
        p=plt.pie(values, labels= ['{:.2f} %'.format(x) for x in values],explode=explode, 
                counterclock=False,
               startangle=180, textprops={'fontsize': 14})
        plt.title('Distribution of annotated elements in {}'.format(rep))
        plt.legend(labels= df.index.tolist(), bbox_to_anchor=(1.07,1.05), loc="upper left") 
        #plt.savefig('/data1/eCLASH/6-bedtools/Overview_targets/PIECHART_{}.pdf'.format(
            #rep), bbox_inches='tight')
        plt.show()
        plt.close()
make_piechart()


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
                'genomic_multiple_mapper', 'genomic_unmappable_reads', 'chimeric read'])
    df=df.rename(index=str, columns={'Nr08_ctrl_1':'ctrl_1', 
                                     'Nr09_ctrl_2':'ctrl_2', 
                                     'Nr10_piRNA_1':'eCLASH_1', 
                                     'Nr11_piRNA_2':'eCLASH_2',
                                     'Nr12_piRNA_3':'eCLASH_3'})
    #df=df/1000000
    #print(df)
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
    #plt.savefig('./UMI_Stacked_bp_mappingandCHIMERIC_PERCENTAGE_all_reps.pdf', bbox_inches='tight')
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
    
    #plt.savefig('./Stacked_bp_mappingevents_PERCENTAGE.pdf', bbox_inches='tight')
    plt.show()
    plt.close
make_stacked_bp_perc()


# In[ ]:




