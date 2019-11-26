#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys


# # Read number

# In[ ]:


read_num = sys.argv[1] # number of read (with ac num) for which I want to generate a picture


# In[ ]:


from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.pairwise2 import format_alignment
from Bio import pairwise2
from Bio import Entrez
import pandas as pd  
import numpy as np
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt


# In[ ]:


def open_reads (filename):
    ids = []
    reads = []
    for i, record in enumerate(SeqIO.parse(filename, "fasta")):
        reads.append(str(record.seq))
        ids.append(record.id)     

        
    reads_dictionary = dict(zip(ids, reads))
    return reads_dictionary 


# In[ ]:

def start_end(dataframe, i):
    start = int(dataframe.iloc[i, 2].split(',')[0][1:])
    end = int(dataframe.iloc[i, 2].split(',')[1][:-1])
    return start, end


class MyCustomTranslator(BiopythonTranslator):

    def compute_feature_color(self, feature):
        if feature.type == "Prediction":
            return "#3498db"
        if feature.type == "source":
            return '#bdc3c7'
        else:
            return "#9b59b6"
        
    def compute_feature_label(self, feature):
        if feature.qualifiers != 'Read' and feature.qualifiers != 'Pred':
            return '{}'.format(feature.qualifiers)
        else:
            return None


# In[1]:


def visualization(read_num):
    item = read_num.split('-')[0]
    reads_file = '{}/Reads_{}.fasta'.format(item, item)

    reads_dictionary = open_reads (reads_file)
    prediction = pd.read_csv('{}/Prediction_{}.csv'.format(item, item), index_col = 'Unnamed: 0')
    ground_truth = pd.read_csv('{}/Ground_truth_{}.csv'.format(item, item), index_col = 'Unnamed: 0')
    sequence_string = reads_dictionary['{}'.format(read_num)]
    sequence_object = Seq(sequence_string, IUPAC.unambiguous_dna)


    record = SeqRecord(sequence_object,
                       id='{}'.format(read_num))
    feature = SeqFeature(FeatureLocation(start=0, end=len(sequence_string), strand=1), type='source', qualifiers='Read')
    record.features.append(feature)

    df_pred = prediction.loc[prediction['Read number'] == '{}'.format(read_num)]
    for i in range(df_pred.shape[0]):
        start_pred, end_pred = start_end(df_pred, i)
        strand_pred = df_pred.iloc[i]['Strand']
        protein = df_pred.iloc[i]['Function']
        feature = SeqFeature(FeatureLocation(start=start_pred, end=end_pred, strand=strand_pred), type='Prediction', qualifiers=protein)
        record.features.append(feature)

    df = ground_truth.loc[ground_truth['Read number'] == '{}'.format(read_num)]    
    for i in range(df.shape[0]):
        start, end = start_end(df, i)
        strand = df.iloc[i]['Strand']
        protein = df.iloc[i]['Function']
        feature = SeqFeature(FeatureLocation(start=start, end=end, strand=strand), type='CDS', qualifiers=protein)
        record.features.append(feature)
    
    graphic_record = MyCustomTranslator().translate_record(record)
    ax, _ = graphic_record.plot(figure_width=20)
    plt.title('Ground truth genes (purple) and predicted genes (blue) in read {}'.format(read_num))    
    plt.savefig('{}/{}.png'.format(item, read_num))    


# In[ ]:


visualization(read_num)

