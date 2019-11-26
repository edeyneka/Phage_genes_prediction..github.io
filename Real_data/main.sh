#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=04:00:00
#SBATCH --mem=16384
#SBATCH --mail-type=END

export PATH=/tudelft.net/staff-umbrella/nb2161/tools/ncbi-blast-2.9.0+/bin:$PATH
export PATH=/tudelft.net/staff-umbrella/nb2161/tools/minimap2-2.16_x64-linux:$PATH
export PATH=/tudelft.net/staff-umbrella/nb2161/tools/miniasm-0.3:$PATH
export PATH=/tudelft.net/staff-umbrella/nb2161/tools/MUMmer3.23:$PATH
export PATH=/tudelft.net/staff-umbrella/nb2161/tools/bwa-0.7.17:$PATH
export PATH=/tudelft.net/staff-umbrella/nb2161/tools/samtools-1.9:$PATH
######Include your code below ###
source ~/python3virt/virtualenv3/bin/activate
location=Real_data
item=File_4_three_reads
intersec_threshold=0.2
merge_threshold=0.2
split_threshold=0.2
merge_overlap_parameter=0.95
gene_threshold=30
python3 Test_pipeline.py $location $item $merge_overlap_parameter $gene_threshold $intersec_threshold $merge_threshold $split_threshold
###### 
