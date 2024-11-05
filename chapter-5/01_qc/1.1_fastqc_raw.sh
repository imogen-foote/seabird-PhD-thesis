#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=12
#SBATCH --mem=10G
#SBATCH --time=0-5:00
#SBATCH --job-name=fastqc
#SBATCH -o /nesi/nobackup/vuw03922/stdout/fastqc.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/fastqc.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module purge
module load FastQC/0.12.1

#set environment
PROJECT=$1
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/illumina
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/fastqc/fastqc_raw
mkdir -p $output_dir

#run program
fastqc $read_dir/*.f*.gz  -t 12 --noextract --outdir $output_dir/


#Cluster: mahuika
#Job ID: 45380559
#State: COMPLETED
#Cores: 6
#Tasks: 1
#Nodes: 1
#Job Wall-time:   61.9%  03:05:43 of 05:00:00 time limit
#CPU Efficiency: 178.8%  1-09:12:05 of 18:34:18 core-walltime
#Mem Efficiency:  51.5%  5.15 GB of 10.00 GB

