#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys


# # Location

# Folder with list of test genomes and reference genomes

# In this folder, there are folders for each test genome

# In[ ]:


import os
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio import Entrez
import pandas as pd  
import numpy as np
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation


# In[31]:


def read_list_acc_num(name):
    genomes_list = [x.strip() for x in open(name, 'r').readlines()]
    return genomes_list


# In[34]:


def retrieve_seq (acc_num): 
    gb_file = "{}/{}.gb".format(acc_num, acc_num)
    entry = SeqIO.read(open(gb_file,"r"), "genbank")
    return entry


# In[35]:


def database (genomes_list):
    d = {}
    for genome in genomes_list:
        d[genome] = str(retrieve_seq(genome).seq)
    return d


# In[36]:


def save_fasta_database (item, database):
    f = open('{}/Reference_db_without_{}.fasta'.format(item, item), "w")
    for genome, seq in database.items():
        if genome != item:
            f.write(">{}\n".format(genome))
            f.write(seq)
            f.write('\n')


# In[37]:


def save_fasta (genome, name):
    f = open(name, "w")
    f.write(">%s\n" %genome.name)
    f.write(str(genome.seq))
    f.write('\n')


# In[38]:


# input - .fastq reads, output - dictionary with DeepSim identifier and my identifier and read itself

# filename - 'acc_num.fastq'
def open_reads (filename, filtering):
    my_identifiers = [] #reads and my identifiers
    ids_deepsim = []
    reads = []
    for i, record in enumerate(SeqIO.parse(filename, "fastq")):
        if filtering != None:
            if len(str(record.seq)) >= filtering:
                my_identifiers.append('{}-{}'.format(filename.split('/')[-1][:-6], (i+1)))
                reads.append(str(record.seq))
                ids_deepsim.append(record.id)     

        
    id_dictionary = dict(zip(ids_deepsim, my_identifiers))
    reads_dictionary = dict(zip(my_identifiers, reads))
    return id_dictionary, reads_dictionary 


# In[39]:


# Save feads in the format for BLAST, input - reads dictionary
def save_reads(reads_dictionary, item):
    f = open("{}/Reads_{}.fasta".format(item, item), "w")
    for i in range(len(reads_dictionary)):
        f.write(">{}\n".format(list(reads_dictionary.keys())[i]))
        f.write(list(reads_dictionary.values())[i])
        f.write('\n')


# In[44]:


def indexes_of_sequence (start, end):
    return pd.Index(np.arange(start, end))


# In[46]:


def intersection(idx_subject, idx_coding_region):
    intersec = list(idx_subject.intersection(idx_coding_region))
    return intersec


# In[47]:


def index_transition(seq, stop):
    i=0
    count = 0
    while count < stop:
        if seq[i]!='-':
            count+=1
        i+=1
    return i


# In[48]:


def alignment_start_end(intersec, s_start, subject_seq):
    overlap_start = intersec[0] # intersection start 0-based
    overlap_end = intersec[-1] # intersection end 
    alignment_start = overlap_start - s_start # where the gene starts in the alignment
    alignment_end = index_transition(subject_seq, overlap_end - s_start+1) # 0-based
    
    return alignment_start, alignment_end


# In[49]:


def aligned_sequence(seq_with_gaps, alignment_start, alignment_end):
    clean_seq = seq_with_gaps[alignment_start:alignment_end].replace('-','')
    return clean_seq


# # Processing the result 

# In[52]:


def set_read_number_as_index(prediction, result, reads_prediction):
    for id_ in reads_prediction:
        result = result.append({'Read number': id_}, ignore_index=True)
    result = result.set_index('Read number')
    return result


# In[53]:


def split_dataframe (dataframe, strand, read_num):
    df = dataframe.loc[dataframe['Read number'] == read_num].sort_values(by=['Start'])
    df = df.loc[df['Strand'] == strand]
    return df


# In[54]:


def initialization_result (result, read_num):
    headers = result.columns.values
    for column_name in headers:
        result.loc[[read_num], [column_name]] = 0
    return result


# In[55]:


def start_end(dataframe, i):
    start = int(dataframe.iloc[i, 2].split(',')[0][1:])
    end = int(dataframe.iloc[i, 2].split(',')[1][:-1])
    return start, end


# In[56]:


def number_of_genes(result, read_num, df_gt, df_pred, strand):
    result.loc[[read_num], ['Number of genes GT {}'.format(str(strand))]] =  df_gt.shape[0]
    result.loc[[read_num], ['Number of predicted genes {}'.format(str(strand))]] =  df_pred.shape[0]
    return result


# In[58]:


def difference (list_gt, list_pred):
    return list(set(list_gt) - set(list_pred) + list(list_pred)- set(list_gt))


# In[64]:


def filling_in_result(result, read_num, df_gt, df_pred, strand):
    
    result = number_of_genes(result, read_num, df_gt, df_pred, strand)
    result = coverage_support(result, read_num, df_gt, df_pred, strand)
    result, indexes_of_not_counted_genes = split(result, read_num, df_gt, df_pred, strand)
    result, indexes_of_not_counted_genes = merge(result, read_num, df_gt, df_pred, strand, indexes_of_not_counted_genes)
    result = edge_error(result, read_num, df_gt, df_pred, strand, indexes_of_not_counted_genes)
    
    return result


# # Process output

# In[66]:



# In[67]:


def average_result(result):
    genome = item
    d = {'Genome': [genome], 'Number of genes GT 1': [0], 'Number of genes GT -1': [0],
     'Number of predicted genes 1': [0], 'Number of predicted genes -1': [0],
     'Coverage 1': [0], 'Coverage -1': [0],
     'Support 1': [0], 'Support -1': [0],
     'Merge 1': [0], 'Merge -1': [0],
     'Split 1': [0], 'Split -1': [0],
     'Error start 1': [0], 'Error start -1': [0],
     'Error end 1': [0], 'Error end -1': [0],
     'Read length': [0]}
    average = pd.DataFrame(data=d)
    average = average.set_index('Genome')
    headers = average.columns.values
    for column_name in headers:
        try:
            average.loc[[genome], [column_name]] = round(np.nanmean(result[[column_name]].values), 3)
        except:
            average.loc[[genome], [column_name]] = None       
    return average


# In[68]:


def save_result(dataframe, name):
    dataframe.to_csv("{}/{}_{}.csv".format(item, name, item))


# # Main code 

# In[1]:


test_list = read_list_acc_num('Test.acc_lst')
for item in test_list:
#     reads_file = 'Test_data/{}/{}.fastq'.format(item, item)
#     mapping_file = 'Test_data/{}/{}.paf'.format(item, item)
#     id_dictionary, reads_dictionary = open_reads (reads_file, 2000)
    result = pd.read_csv('{}/Initial_result_{}.csv'.format(item, item), index_col = 'Read number')
    average = average_result(result)
    save_result(average, 'Average_result')

    print('Average result {} Done'.format(item))

average = pd.read_csv('{}/Average_result_{}.csv'.format(test_list[0], test_list[0]), index_col = 'Genome')
for item in test_list[1:]:
    df = pd.read_csv('{}/Average_result_{}.csv'.format(item, item), index_col = 'Genome')
    average = pd.concat([average, df])
    
f = open('Average_result.txt', 'w')
f.write(average.to_latex())
average.to_csv("Average_result.csv")