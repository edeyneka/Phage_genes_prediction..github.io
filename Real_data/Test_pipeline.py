#!/usr/bin/env python
# coding: utf-8

# In[30]:


import sys


# # Location

# Folder with list of test genomes and reference genomes

# In this folder, there are .bash file and folders for each test genome

# In[ ]:


location = sys.argv[1] # folder = number of experiment or other information
item = sys.argv[2] # filemane with reads
merge_overlap_parameter = float(sys.argv[3])
gene_threshold = int(sys.argv[4])
intersec_threshold = sys.argv[5]
merge_threshold = sys.argv[6]
split_threshold = sys.argv[7]


# In[ ]:

from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

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
from dna_features_viewer import BiopythonTranslator



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
        d[genome] = str(retrieve_seq_db(genome).seq)
    return d


# In[36]:


def save_fasta_database (database):
    f = open('/tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/Reference_data/Reference_db_all_genomes.fasta', "w")
    for genome, seq in database.items():
        f.write(">{}\n".format(genome))
        f.write(seq)
        f.write('\n')


# Create database

# In[ ]:


#reference_list = read_list_acc_num('Reference.acc_lst')
#db = database (reference_list)
#print("DB")

# In[37]:


def save_fasta (genome, name):
    f = open(name, "w")
    f.write(">%s\n" %genome.name)
    f.write(str(genome.seq))
    f.write('\n')


# In[1]:


# input - .fastq reads, output - dictionary with DeepSim identifier and my identifier and read itself

# filename - 'acc_num.fastq'
def open_reads (filename):
    ids = []
    reads = []
    for i, record in enumerate(SeqIO.parse(filename, "fastq")):
        reads.append(str(record.seq))
        ids.append(record.id)     


    reads_dictionary = dict(zip(ids, reads))
    return reads_dictionary 


# In[39]:


# Save feads in the format for BLAST, input - reads dictionary
def save_reads(reads_dictionary, item):
    f = open("{}/Reads.fasta".format(item), "w")
    for i in range(len(reads_dictionary)):
        f.write(">{}\n".format(list(reads_dictionary.keys())[i]))
        f.write(list(reads_dictionary.values())[i])
        f.write('\n')


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
        df.iloc[i]['s. start'] = str(int(df.iloc[i]['s. start'])-1) # 0-based
        df.iloc[i]['query seq'] = df.iloc[i]['query seq'].replace('\n','')
        df.iloc[i]['subject seq'] = df.iloc[i]['subject seq'].replace('\n','')
        
    return df


# In[44]:


def indexes_of_sequence (start, end):
    return pd.Index(np.arange(start, end))


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
        result.loc[[read_num], ['Error start {}'.format(str(strand))]] = 'NaN'
        result.loc[[read_num], ['Error end {}'.format(str(strand))]] = 'NaN'
        
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


def save_result(dataframe, item, name):
    dataframe.to_csv("{}/{}.csv".format(item, name))


# Creating database with all genomes

# In[ ]:


#save_fasta_database (db)
#cmd0 = 'makeblastdb -in /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/Reference_data/Reference_db_all_genomes.fasta -dbtype nucl -out /scratch/ekaterina/Reference_db_all_genomes.db'
#os.system(cmd0)


# Generating Alignment file

# In[ ]:

cmd3 = 'mkdir /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/{}/{}'.format(location, item)
os.system(cmd3)
reads_file = '{}.fastq'.format(item)
reads_dictionary = open_reads (reads_file)
save_reads(reads_dictionary, item)
# cmd2 = 'mkdir /scratch/ekaterina/{}'.format(location)
# os.system(cmd2)
cmd4 = 'blastn -db /scratch/ekaterina/Reference_db_all_genomes.db -query {}/Reads.fasta -out {}/Alignment_result.txt -strand both -outfmt "7 qacc sacc evalue qstart qend sstart send bitscore score pident nident qframe sframe qcovs qseq sseq"'.format(item, item)
os.system(cmd4)


# In[ ]:


prediction = prediction_dataframe('{}/Alignment_result.txt'.format(item))
save_result(prediction, item, 'Prediction')
print('Prediction {} Done'.format(item))


# Visualization

# In[ ]:


class MyCustomTranslator(BiopythonTranslator):

    def compute_feature_color(self, feature):
        if feature.type == "Prediction":
            return "#3498db"
        elif feature.type == "source":
            return '#bdc3c7'
        else:
            return "#ccccff"
        
    def compute_feature_label(self, feature):
        if feature.qualifiers != 'Read':
            return '{}'.format(feature.qualifiers)
        else:
            return None


# In[1]:


def visualization(item, read_num, seq):
    
    sequence_string = seq
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

    graphic_record = MyCustomTranslator().translate_record(record)
    ax, _ = graphic_record.plot(figure_width=20)
    plt.title('Predicted genes (blue) in read {}'.format(read_num))    
    plt.savefig('{}/{}.png'.format(item, read_num))    


# In[ ]:


for read_num, seq in reads_dictionary.items():
    visualization(item, read_num, seq)


# In[ ]:




