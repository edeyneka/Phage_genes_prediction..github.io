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
location=
coverage=5
filtering=2000
intersec_threshold=0.2
merge_threshold=0.2
split_threshold=0.2
merge_overlap_parameter=0.95
gene_threshold=30
#python3 Create_test_data.py $location
#python3 Create_reference_database.py $location
#deactivate
#export PATH=$PATH:/home/nfs/edeyneka/anaconda2/bin
#file="/tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/$location/Test.acc_lst"
#cd ~/DeepSimulator
#while IFS= read -r line
#do
	#mkdir /scratch/ekaterina/DeepSim
	#./deep_simulator.sh -K $coverage -i /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/$location/"$line"/"$line".fasta -o /scratch/ekaterina/DeepSim
	#cp /scratch/ekaterina/DeepSim/test.fastq /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/$location/"$line"/"$line".fastq
	#cp /scratch/ekaterina/DeepSim/mapping.paf /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/$location/"$line"/"$line".paf
	#rm -rf /scratch/ekaterina/DeepSim
	#echo $line
#done <"$file"
#source ~/python3virt/virtualenv3/bin/activate
#cd /tudelft.net/staff-umbrella/abeellab/ekaterina/Experiment/$location
#python3 Generating_ground_truth.py $location $filtering
#python3 Blast_reads.py $location
#python3 Generating_prediction.py $location $filtering $merge_overlap_parameter $gene_threshold
python3 Result.py $intersec_threshold $merge_threshold $split_threshold
python3 Average_result.py
###### 
