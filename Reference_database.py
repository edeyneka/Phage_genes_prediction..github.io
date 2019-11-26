#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


# In[3]:


def read_list_acc_num(name):
    genomes_list = [x.strip() for x in open(name, 'r').readlines()]
    return genomes_list


# In[ ]:


def save_genebank(acc_num):
    Entrez.email = "E.Deyneka@student.tudelft.nl"
    handle = Entrez.efetch(db="nuccore", id=acc_num, rettype="gb", retmode="text")
    entry = SeqIO.read(handle, "genbank")
    output_file = open('/tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/Reference_data/{}/{}.gb'.format(acc_num, acc_num), 'w')
    SeqIO.write(entry, output_file, 'genbank')


# In[ ]:


reference_list = read_list_acc_num('Reference.acc_lst')


# In[ ]:


for item in reference_list:
    cmd = 'mkdir /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/Reference_data/{}'.format(item)
    os.system(cmd)
    try:
        save_genebank(item)
        print('{} Done'.format(item))
    except:
        print('Error {}'.format(item))


# In[ ]:




