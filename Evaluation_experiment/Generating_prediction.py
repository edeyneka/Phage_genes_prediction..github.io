#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys


# # Location

# Folder with list of test genomes and reference genomes

# In this folder, there are folders for each test genome

# In[ ]:


location = sys.argv[1]
filtering = int(sys.argv[2])
merge_overlap_parameter = float(sys.argv[3])
gene_threshold = int(sys.argv[4])


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
def open_reads (filename):
    ids = []
    reads = []
    for i, record in enumerate(SeqIO.parse(filename, "fasta")):
        reads.append(str(record.seq))
        ids.append(record.id)     

        
    reads_dictionary = dict(zip(ids, reads))
    return reads_dictionary

# In[39]:


# Save feads in the format for BLAST, input - reads dictionary
def save_reads(reads_dictionary, item):
    f = open("{}/Reads_{}.fasta".format(item, item), "w")
    for i in range(len(reads_dictionary)):
        f.write(">{}\n".format(list(reads_dictionary.keys())[i]))
        f.write(list(reads_dictionary.values())[i])
        f.write('\n')


# # Alignment

# ### Inputs: alignment and ground truth 

# In[42]:


def parse_alignment(alignment_file):
    file = open(alignment_file, "r").readlines()
    alignment = []
    for string in file:
        if '#' not in string:
            alignment.append(string.split('\t'))

    for string in file:
        if '# Fields' in string:
            columns = string[10:-1].split(', ')  
    try:
        df = pd.DataFrame(alignment, columns = columns)
        return df
    except:
        return None
    


# In[43]:


def alignment_to_zero_based(df):
    for i in range(df.shape[0]):
        df.iloc[i]['q. start'] = str(int(df.iloc[i]['q. start'])-1) # 0-based
        if int(df.iloc[i]['s. start']) < int(df.iloc[i]['s. end']):
            df.iloc[i]['s. start'] = str(int(df.iloc[i]['s. start'])-1) # 0-based
        else: 
            df.iloc[i]['s. end'] = str(int(df.iloc[i]['s. end'])-1) # 0-based
        df.iloc[i]['query seq'] = df.iloc[i]['query seq'].replace('\n','')
        df.iloc[i]['subject seq'] = df.iloc[i]['subject seq'].replace('\n','')
        
    return df

# In[44]:


def indexes_of_sequence (start, end):
    return pd.Index(np.arange(start, end))

def retrieve_seq_db (acc_num): 
    gb_file = "/tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/Reference_data/{}/{}.gb".format(acc_num, acc_num)
    entry = SeqIO.read(open(gb_file,"r"), "genbank")
    return entry

# In[45]:


def initialization_prediction(df_q, i):
    subject_acc = df_q.iloc[i]['subject acc.'] 
    query_acc = df_q.iloc[i]['query acc.']
    
    subject_entry = retrieve_seq_db (subject_acc)
    
    subject_strand = int(df_q.iloc[i]['sbjct frame'])
    
    s_start = int(df_q.iloc[i]['s. start'])
    s_end = int(df_q.iloc[i]['s. end'])
    
    q_start = int(df_q.iloc[i]['q. start'])
    q_end = int(df_q.iloc[i]['q. end'])
    
    if s_end < s_start:
        idx_subject = indexes_of_sequence (s_end,s_start) # indexes that correspond to the alignment in the reference 
    else:
        idx_subject = indexes_of_sequence (s_start,s_end) # indexes that correspond to the alignment in the reference 
    
    subject_seq = df_q.iloc[i]['subject seq']
    query_seq = df_q.iloc[i]['query seq']
    

    read = reads_dictionary[query_acc]
    
    return subject_entry, idx_subject, subject_strand, s_start, subject_seq, query_acc, query_seq

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


# In[3]:


def gene_within_mapping(subject_entry, idx_subject, subject_strand, s_start, subject_seq, query_acc, query_seq, read, prediction):
    for j in range(len(subject_entry.features)): # iterate over features
        # here we take into accout the strand and coding region
        if subject_entry.features[j].type=='CDS':
            #now we check if the mapped region is within the coding region
            idx_coding_region = indexes_of_sequence(subject_entry.features[j].location.start.position, 
                                          subject_entry.features[j].location.end.position) # indexes of features
            intersec = intersection(idx_subject, idx_coding_region)
            if len(intersec) != 0 : # if they intersect with indexes in alignment
                alignment_start, alignment_end = alignment_start_end(intersec, s_start, subject_seq)
                query_gene = aligned_sequence(query_seq, alignment_start, alignment_end)
                subject_gene = aligned_sequence(subject_seq, alignment_start, alignment_end)
                gene_location = [read.find(query_gene), read.find(query_gene)+len(query_gene)]
                if (gene_location[1] - gene_location[0] + 1) >= 30:
                    if subject_entry.features[j].location.strand == subject_strand:
                        strand = 1
                    else:
                        strand = -1
                    prediction = prediction.append({'Read number': query_acc, 'Strand': strand, 
                                                'Location of gene in the read': '['+str(gene_location[0]) +',' 
                                                                     +str(gene_location[1])+')', 'Read length': len(read), 
                                                    'Start':gene_location[0], 'End': gene_location[1], 
                                                    'Function': subject_entry.features[j].qualifiers['product'][0]}, ignore_index=True) # 0-based system
    return prediction

# In[ ]:


def split_dataframe (dataframe, strand, read_num):
    df = dataframe.loc[dataframe['Read number'] == read_num].sort_values(by=['Start'])
    df = df.loc[df['Strand'] == strand]
    return df


# In[4]:


def merge_overlaps(prediction, strand, read_num):
    df = split_dataframe (prediction, strand, read_num)
    
    for i in list(df.index):
        start_1 = df.loc[i]['Start']
        end_1 = df.loc[i]['End']
        idx_1 = indexes_of_sequence (start_1, end_1) 
        idx_1_list = list(idx_1)
        for j in list(df.index):
            if i != j:
                start_2 = df.loc[j]['Start']
                end_2 = df.loc[j]['End']
                idx_2 = indexes_of_sequence (start_2, end_2) 
                idx_2_list = list(idx_2)
                intersec = intersection(idx_1, idx_2)
                if len(intersec) != 0 and len(intersec)/len(idx_1_list) >= merge_overlap_parameter:
                    if len(idx_1_list) > len(idx_2_list):
                        try:
                            prediction = prediction.drop([j])
                        except:
                            pass
                    else:
                        try:
                            prediction = prediction.drop([i])
                        except:
                            pass      
    return prediction


# In[5]:


def filtering_df(prediction): 
    _, idx = np.unique(prediction[['Read number']].values, return_index=True)
    reads_prediction = [x[0] for x in prediction[['Read number']].values[np.sort(idx)]]
    for read_num in reads_prediction:
        prediction = merge_overlaps(prediction, 1, read_num)
        prediction = merge_overlaps(prediction, -1, read_num)
    return prediction


# In[6]:


def prediction_dataframe(alignment_file):
    d = {'Read number': [], 'Strand': [], 'Location of gene in the read': [], 'Read length': [], 'Start':[], 'End': [], 'Function': []}
    prediction = pd.DataFrame(data=d)
    alignment_dataframe = parse_alignment(alignment_file)
    try:
        if alignment_dataframe == None:
            return prediction
    except: 
        alignment_dataframe = alignment_to_zero_based(alignment_dataframe)



        #extracting alignment for particular read
        _, idx = np.unique(alignment_dataframe[['query acc.']].values, return_index=True)

        for query in [x[0] for x in alignment_dataframe[['query acc.']].values[np.sort(idx)]]: # if the list is not complete, then there are no hits for this read
            df_q = alignment_dataframe.loc[alignment_dataframe['query acc.'] == query] # dataframe for one query acc.

            for i in range(df_q.shape[0]): # iterate over hits within one read
                subject_entry, idx_subject, subject_strand, s_start, subject_seq, query_acc, query_seq = initialization_prediction(df_q, i)

                prediction = gene_within_mapping(subject_entry, idx_subject, subject_strand, s_start, subject_seq, query_acc, query_seq, reads_dictionary[query_acc], prediction)

        prediction = filtering_df(prediction)       
        return prediction


# In[69]:


def save_result(dataframe, name):
    dataframe.to_csv("{}/{}_{}.csv".format(item, name, item))


# # Main code 

# In[1]:


test_list = read_list_acc_num('Test.acc_lst')
for item in test_list:

    reads_file = '{}/Reads_{}.fasta'.format(item, item)
    reads_dictionary = open_reads (reads_file)
    prediction = prediction_dataframe('{}/Alignment_result_{}.txt'.format(item, item))
    save_result(prediction, 'Prediction')
    
    print('Prediction {} Done'.format(item))

