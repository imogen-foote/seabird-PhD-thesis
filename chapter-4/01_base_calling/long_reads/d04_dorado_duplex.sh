#!/bin/bash -e
#SBATCH --gpus-per-node=A100:1
#SBATCH --partition=hgx
#SBATCH --cpus-per-task=4
#SBATCH --profile=task
#SBATCH --mem=200G
#SBATCH --time=3-12:00:00 
#SBATCH --job-name=dorado_duplex
#SBATCH -o /nesi/nobackup/vuw03922/stdout/dorado_duplex.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/dorado_duplex.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module purge
module load Dorado/0.4.1

#set up
PROJECT=$1
LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/pod5
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/duplex_bam

#run dorado basecaller
dorado duplex --device 'cuda:all' dna_r10.4.1_e8.2_400bps_sup@v4.1.0 $read_dir/ > $out_dir/$LIB.bam

#create sequencing summary file
dorado summary $out_dir/$LIB.bam > $out_dir/sequencing_summary_$LIB.txt


#Cluster: mahuika
#Job ID: 40556564
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   54.9%  2-17:53:46 of 5-00:00:00 time limit
#CPU Efficiency:  52.9%  2-21:42:05 of 5-11:47:32 core-walltime
#Mem Efficiency:  74.0%  147.98 GB of 200.00 GB


#Cluster: mahuika
#Job ID: 40556553
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   52.8%  2-15:21:11 of 5-00:00:00 time limit
#CPU Efficiency:  53.9%  2-20:21:19 of 5-06:42:22 core-walltime
#Mem Efficiency:  77.0%  154.05 GB of 200.00 GB

