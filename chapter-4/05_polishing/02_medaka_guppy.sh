#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=80G
#SBATCH --time=1-12:00:00
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

draft=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/2_polished/$assembly/consensus.fasta
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/2_polished/$assembly'_medaka_2'

mkdir -p $out_dir/tmp_fastq
#antipodean
#cp /nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$sample'_lib1'/fastq/$sample*$filter.fastq.gz $out_dir/tmp_fastq
#cp /nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$sample'_lib2'/fastq/$sample*$filter.fastq.gz $out_dir/tmp_fastq

#gibson
cp /nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$sample'_lib3'/fastq/$sample*$filter.fastq.gz $out_dir/tmp_fastq
cp /nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$sample'_lib4'/fastq/$sample*$filter.fastq.gz $out_dir/tmp_fastq

cat $out_dir/tmp_fastq/*.fastq.gz | gunzip -c > $out_dir/tmp_fastq/combined.reads.fastq

#copy assembly file to tmp dir as advised https://github.com/nanoporetech/medaka/issues/476#issuecomment-1830214$
cp $draft $out_dir/tmp_fastq

#run medaka
medaka_consensus -d $out_dir/tmp_fastq/consensus.fasta \
        -i $out_dir/tmp_fastq/combined.reads.fastq \
        -o $out_dir \
        -m r104_e81_sup_g610
#       -t ${NPROC} \

#rm -r $out_dir/tmp_fastq

conda deactivate
  
#Cluster: mahuika
#Job ID: 41977094
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:   25.5%  1-00:26:47 of 4-00:00:00 time limit
#CPU Efficiency:  20.4%  1-00:54:59 of 5-02:13:55 core-walltime
#Mem Efficiency:  63.0%  62.98 GB of 100.00 GB

