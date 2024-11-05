#!/bin/bash 
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=0-12:00:00
#SBATCH --job-name=porechop
#SBATCH -o /nesi/nobackup/vuw03922/stdout/porechop.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/porechop.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#call modules
module purge
module load Porechop/0.2.4-gimkl-2020a-Python-3.8.2

#set up
PROJECT=$1
LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/duplex_fastq

porechop -i $read_dir/*.fastq.gz -o $read_dir/$LIB.trimmed_duplex.fastq.gz --threads 40


#Cluster: mahuika
#Job ID: 40339882
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:   79.2%  09:30:18 of 12:00:00 time limit
#CPU Efficiency: 113.2%  1-19:01:42 of 1-14:01:12 core-walltime
#Mem Efficiency:  90.0%  144.02 GB of 160.00 GB

