#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10 
#SBATCH --mem=25G
#SBATCH --time=4-0:00:00
#SBATCH --job-name=medaka
#SBATCH -o /nesi/nobackup/vuw03922/stdout/medaka.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/medaka.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules and activate conda env
module purge && module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/vuw03922/medaka

PROJECT=$1
assembly=$2 #e.g. blue-45G_filtered_h50 etc
sample=$3 #e.g. antipodean_sa52* or gibsons_d36*
filter=$4 #e.g. filtered_50 or filtered_q8_h50 etc
#NPROC=$(nproc)#prints the number of processing units available to the current process

draft=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/1_raw/$assembly/assembly.fasta
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/2_polished/$assembly'_medaka_1'

mkdir -p $out_dir/tmp_fastq
cp /nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$sample'_lib5'/fastq/$sample*$filter.fastq.gz $out_dir/tmp_fastq
cp /nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$sample'_lib6'/fastq/$sample*$filter.fastq.gz $out_dir/tmp_fastq
cat $out_dir/tmp_fastq/*.fastq.gz | gunzip -c > $out_dir/tmp_fastq/combined.reads.fastq

#copy assembly file to tmp dir as advised https://github.com/nanoporetech/medaka/issues/476#issuecomment-1830214360
cp $draft $out_dir/tmp_fastq

#run medaka
medaka_consensus -d $out_dir/tmp_fastq/assembly.fasta \
	-i $out_dir/tmp_fastq/combined.reads.fastq \
	-o $out_dir \
	-m r1041_e82_400bps_sup_v4.1.0
#	-t ${NPROC} \

#rm -r $out_dir/tmp_fastq  

conda deactivate


#Cluster: mahuika
#Job ID: 41857432
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:   34.7%  1-17:38:12 of 5-00:00:00 time limit
#CPU Efficiency:  20.8%  1-19:22:07 of 8-16:11:00 core-walltime
#Mem Efficiency:  53.3%  13.32 GB of 25.00 GB

