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


# ### Input: two lists of accesion numbers. First - for reference database, second - for test 

# Save reference genomes and make database

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


# # Parsing mapping file

# In[1]:


def mapping(filename, id_dictionary):
    df = pd.read_table(filename, header=None)
    d = {'Read number': [], 'Strand': [], 'Location of gene in the read': [], 'Read length':[], 'Start': [], 'End': [], 'Function':[] }
    ground_truth = pd.DataFrame(data=d)

    for i in range(df.shape[0]):
        if id_dictionary.get(df.iloc[i,0]) != None:
            start_for = df.iloc[i,7]
            end_for = df.iloc[i,8]

            start_rev = df.iloc[i,6]-df.iloc[i,8]
            end_rev = df.iloc[i,6] - (df.iloc[i,7]-1)
            subject_acc = df.iloc[i,5] 


            #retrieving subject genome
            subject_genome = retrieve_seq (subject_acc)

            #location of mapped regions
            idx_ref_for = pd.Index(np.arange(start_for, end_for))
            idx_ref_rev = pd.Index(np.arange(start_rev, end_rev))
            for j in range(len(subject_genome.features)):
                if subject_genome.features[j].type=='CDS':
                    # for the forward strand
                    ground_truth = gene_locations(df.iloc[i,0], ground_truth, subject_genome.features[j], idx_ref_for, start_for, end_for)

                    # for the reverse strand
                    #ground_truth = gene_locations(df.iloc[i,0], ground_truth, subject_genome.features[j], idx_ref_rev, -1, start_rev, end_rev)
    return ground_truth


# In[2]:


def gene_locations(read_id, ground_truth, subject_genome_feature, idx_ref, start, end): # subject_genome_features[j]
    idx = pd.Index(np.arange(subject_genome_feature.location.start.position, 
                             subject_genome_feature.location.end.position))
    strand = subject_genome_feature.strand
    overlap = list(idx_ref.intersection(idx))
    if len(overlap) != 0:
        overlap_start = overlap[0]
        overlap_end = overlap[-1]
        ground_truth = ground_truth.append({'Read number': id_dictionary.get(read_id), 'Strand': strand, 
                                'Location of gene in the read': '['+str(overlap_start - start) +',' 
                                                                 +str(overlap_end - start+1)+')', 
                                            'Start':overlap_start - start, 'End':overlap_end - start+1,
                                            'Read length': len(reads_dictionary[id_dictionary.get(read_id)]),
                                           'Function':subject_genome_feature.qualifiers['product'][0]} , ignore_index=True) # 0-based system
    return ground_truth


# In[ ]:


def save_result(dataframe, name):
    dataframe.to_csv("{}/{}_{}.csv".format(item, name, item))


# # Main code 

# In[1]:


test_list = read_list_acc_num('Test.acc_lst')
for item in test_list:
    reads_file = '{}/{}.fastq'.format(item, item)
    mapping_file = '{}/{}.paf'.format(item, item)
    id_dictionary, reads_dictionary = open_reads (reads_file, filtering)
    ground_truth = mapping(mapping_file, id_dictionary)
    save_result(ground_truth, 'Ground_truth')
    
    print('Ground truth {} Done'.format(item))

