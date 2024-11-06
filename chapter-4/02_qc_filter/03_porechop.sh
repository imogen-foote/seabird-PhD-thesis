#!/bin/bash 
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=8
#SBATCH --mem=350G
#SBATCH --time=1-00:00:00
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
#LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/all_fastq

porechop -i $read_dir/$PROJECT.simplex.fastq.gz -o $read_dir/$PROJECT.simplex_trimmed.fastq.gz --threads 40
porechop -i $read_dir/$PROJECT.duplex.fastq.gz -o $read_dir/$PROJECT.duplex_trimmed.fastq.gz --threads 40

#Cluster: mahuika
#Job ID: 39383318
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:   69.4%  05:33:19 of 08:00:00 time limit
#CPU Efficiency: 116.7%  1-01:56:35 of 22:13:16 core-walltime
#Mem Efficiency:  71.3%  85.59 GB of 120.00 GB


#Cluster: mahuika
#Job ID: 39383319
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:   89.5%  07:09:43 of 08:00:00 time limit
#CPU Efficiency: 113.0%  1-08:22:02 of 1-04:38:52 core-walltime
#Mem Efficiency:  90.6%  108.69 GB of 120.00 GB 

#Cluster: mahuika
#Job ID: 39765114
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:   83.1%  08:18:38 of 10:00:00 time limit
#CPU Efficiency: 115.4%  1-14:22:19 of 1-09:14:32 core-walltime
#Mem Efficiency:  99.3%  139.01 GB of 140.00 GB


#Cluster: mahuika
#Job ID: 39774920
#State: TIMEOUT
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:  100.1%  10:00:22 of 10:00:00 time limit
#CPU Efficiency:  91.5%  1-12:37:39 of 1-16:01:28 core-walltime
#Mem Efficiency:  89.7%  143.47 GB of 160.00 GB

