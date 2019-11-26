#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys


# # Location

# Folder with list of test genomes and reference genomes

# In this folder, there are folders for each test genome

# In[ ]:


location = sys.argv[1]


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


# In[ ]:


def retrieve_seq_db (acc_num): 
    gb_file = "/tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/Reference_data/{}/{}.gb".format(acc_num, acc_num)
    entry = SeqIO.read(open(gb_file,"r"), "genbank")
    return entry


# In[35]:


def database (genomes_list):
    d = {}
    for genome in genomes_list:
        d[genome] = str(retrieve_seq_db (genome).seq)
    return d


# In[36]:


def save_fasta_database (item, database):
    f = open('{}/Reference_db_without_{}.fasta'.format(item, item), "w")
    for genome, seq in database.items():
        if genome != item:
            f.write(">{}\n".format(genome))
            f.write(seq)
            f.write('\n')


# Create database

# In[ ]:


reference_list = read_list_acc_num('Reference.acc_lst')
db = database (reference_list)


# Save test genomes and create database without a test genome

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


# In[ ]:


test_list = read_list_acc_num('Test.acc_lst')
for item in test_list:
    save_fasta_database (item, db)
    cmd1 = 'makeblastdb -in /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/{}/{}/Reference_db_without_{}.fasta -dbtype nucl -out /scratch/ekaterina/{}/{}/Reference_db_without_{}.db'.format(location, item,  item, location, item, item)
    os.system(cmd1)

