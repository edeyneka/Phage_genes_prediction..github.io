#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys


# # Location

# Folder with list of test genomes and reference genomes

# In this folder, there are folders for each test genome

# In[ ]:


intersec_threshold = float(sys.argv[1])
merge_threshold = float(sys.argv[2])
split_threshold = float(sys.argv[3])


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


# In[57]:

def indexes_of_sequence (start, end):
    return pd.Index(np.arange(start, end))

def edge_error(result, read_num, df_gt, df_pred, strand, indexes_of_not_counted_genes):
#     intersec_threshold = 0.2
    for i in range(df_gt.shape[0]):

        start, end = start_end(df_gt, i)
        idx_gt = indexes_of_sequence (start, end) 
        gt = list(idx_gt)
        
        for j in indexes_of_not_counted_genes:
            start_pred, end_pred = start_end(df_pred, j)
            idx_pred = indexes_of_sequence (start_pred, end_pred) 
            pred = list(idx_pred)
        
            intersec = intersection(idx_gt, idx_pred)
            if len(intersec) != 0 and len(intersec)/len(pred) >= intersec_threshold and len(intersec)/len(gt) >= intersec_threshold:
                result.loc[[read_num], ['Error start {}'.format(str(strand))]] += int(np.abs(start - start_pred))
                result.loc[[read_num], ['Error end {}'.format(str(strand))]] += int(np.abs(end - end_pred))
                
    return result


# In[58]:


def difference (list_gt, list_pred):
    return list(set(list_gt) - set(list_pred) + list(list_pred)- set(list_gt))


# In[59]:


def coverage_support(result, read_num, df_gt, df_pred, strand):
    gt = []
    pred = []
    for i in range(df_gt.shape[0]):
        
        start, end = start_end(df_gt, i)
        gt += list(indexes_of_sequence (start, end))  # 0-based system

    for j in range(df_pred.shape[0]):

        start_pred, end_pred = start_end(df_pred, j)
        pred += list (indexes_of_sequence (start_pred, end_pred))
    try:
        result.loc[[read_num], ['Coverage {}'.format(str(strand))]] = len(list(set(gt).intersection(set(pred))))/len(gt)
        result.loc[[read_num], ['Support {}'.format(str(strand))]] = len(list(set(gt).intersection(set(pred))))/len(pred)
    except:
        pass
            
    return result        


# In[7]:


def merge (result, read_num, df_gt, df_pred, strand, indexes_of_not_counted_genes):
#     merge_threshold = 0.2 # threshold can be changed
    i_start = 1e10
    i_end = 0
    for j in range(df_pred.shape[0]):
        merge = 0
        start_pred, end_pred = start_end(df_pred, j)
        idx_pred = indexes_of_sequence (start_pred, end_pred) 
        pred = list(idx_pred)
        
        for i in range(df_gt.shape[0]):

            start, end = start_end(df_gt, i)
            idx_gt = indexes_of_sequence (start, end) 
            gt = list(idx_gt)
            
            intersec = intersection(idx_gt, idx_pred)
            if len(intersec) != 0:
                if i < i_start:
                    i_start = i
                if i > i_end:
                    i_end = i
                if len(intersec)/len(pred) >= merge_threshold and len(intersec)/len(gt) >= merge_threshold:
                    result.loc[[read_num], ['Merge {}'.format(str(strand))]] += 1
                    merge += 1

                    
        if merge == 1:
            result.loc[[read_num], ['Merge {}'.format(str(strand))]] -= 1
            
        if merge != 1 and merge!=0 and i_start != 1e10: 
            result.loc[[read_num], ['Merge {}'.format(str(strand))]] -= 1
            start, _ = start_end(df_gt, i_start)
            _, end = start_end(df_gt, i_end)
            try:
                indexes_of_not_counted_genes.remove(i_start)
            except:
                pass
            try:
                indexes_of_not_counted_genes.remove(i_end)
            except:
                pass
            result.loc[[read_num], ['Error start {}'.format(str(strand))]] += int(np.abs(start - start_pred))
            result.loc[[read_num], ['Error end {}'.format(str(strand))]] += int(np.abs(end - end_pred))

    return result, indexes_of_not_counted_genes


# In[8]:


def split (result, read_num, df_gt, df_pred, strand):
#     split_threshold = 0.2 # threshold can be changed
    
    j_start = 1e10 # assignment because if there is 
    j_end = 0
    indexes_of_not_counted_genes = list(range(df_pred.shape[0]))
    for i in range(df_gt.shape[0]):

        split = 0
        start, end = start_end(df_gt, i)
        idx_gt = indexes_of_sequence (start, end) 
        gt = list(idx_gt)
        
        for j in range(df_pred.shape[0]):
            start_pred, end_pred = start_end(df_pred, j)
            idx_pred = indexes_of_sequence (start_pred, end_pred) 
            pred = list(idx_pred)
        
            intersec = intersection(idx_gt, idx_pred)
            if len(intersec) != 0:
                if j < j_start:
                    j_start = j
                if j > j_end:
                    j_end = j
                if len(intersec)/len(gt) >= split_threshold and len(intersec)/len(gt) >= split_threshold:
                    result.loc[[read_num], ['Split {}'.format(str(strand))]] += 1
                    split += 1

        
        if split == 1:
            result.loc[[read_num], ['Split {}'.format(str(strand))]] -= 1
            
        if split != 1 and split !=0 and j_start != 1e10:   
            result.loc[[read_num], ['Split {}'.format(str(strand))]] -= 1

            start_pred, _ = start_end(df_pred, j_start)
            _, end_pred = start_end(df_pred, j_end)
            try:
                indexes_of_not_counted_genes.remove(j_start)
            except:
                pass
            try:
                indexes_of_not_counted_genes.remove(j_end)
            except:
                pass

            result.loc[[read_num], ['Error start {}'.format(str(strand))]] += int(np.abs(start - start_pred))
            result.loc[[read_num], ['Error end {}'.format(str(strand))]] += int(np.abs(end - end_pred))
            
    if j_start == 1e10:
        result.loc[[read_num], ['Error start {}'.format(str(strand))]] = None
        result.loc[[read_num], ['Error end {}'.format(str(strand))]] = None
        
    return result, indexes_of_not_counted_genes


# In[64]:


def filling_in_result(result, read_num, df_gt, df_pred, strand):
    
    result = number_of_genes(result, read_num, df_gt, df_pred, strand)
    result = coverage_support(result, read_num, df_gt, df_pred, strand)
    result, indexes_of_not_counted_genes = split(result, read_num, df_gt, df_pred, strand)
    result, indexes_of_not_counted_genes = merge(result, read_num, df_gt, df_pred, strand, indexes_of_not_counted_genes)
    result = edge_error(result, read_num, df_gt, df_pred, strand, indexes_of_not_counted_genes)
    
    return result


# In[65]:


def result_dataframe(ground_truth, prediction):
    d = {'Read number': [], 'Number of genes GT 1': [], 'Number of genes GT -1': [],
     'Number of predicted genes 1': [], 'Number of predicted genes -1': [],
     'Coverage 1': [], 'Coverage -1': [],
     'Support 1': [], 'Support -1': [],
     'Merge 1': [], 'Merge -1': [],
     'Split 1': [], 'Split -1': [],
     'Error start 1': [], 'Error start -1': [],
     'Error end 1': [], 'Error end -1': [],
     'Read length': []}
    result = pd.DataFrame(data=d)
    
    # extracting prediction for a particular read
    _, idx = np.unique(prediction[['Read number']].values, return_index=True)
    # list of reads for which there is a prediction
    reads_prediction = [x[0] for x in prediction[['Read number']].values[np.sort(idx)]]
    
    
    result = set_read_number_as_index(prediction, result, reads_prediction)

    for read_num in reads_prediction:

        # copying read length in result dataframe

        result = initialization_result (result, read_num)
        result.loc[[read_num], ['Read length']] = prediction.loc[prediction['Read number'] == read_num][['Read length']].values[0][0]
        
        df_gt_for = split_dataframe (ground_truth, 1, read_num)
        df_pred_for = split_dataframe (prediction, 1, read_num)
        
        df_gt_rev = split_dataframe (ground_truth, -1, read_num)
        df_pred_rev = split_dataframe (prediction, -1, read_num)
        
        result = filling_in_result(result, read_num, df_gt_for, df_pred_for, 1)
        result = filling_in_result(result, read_num, df_gt_rev, df_pred_rev, -1)
        
    return result


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
    ground_truth = pd.read_csv('{}/Ground_truth_{}.csv'.format(item, item), index_col = 'Unnamed: 0')
    prediction = pd.read_csv('{}/Prediction_{}.csv'.format(item, item), index_col = 'Unnamed: 0')
    result = result_dataframe(ground_truth, prediction)
    save_result(result, 'Initial_result')
    
    print('Result {} Done'.format(item))

