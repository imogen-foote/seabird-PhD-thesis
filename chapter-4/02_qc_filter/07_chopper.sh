#!/bin/bash 
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=8
#SBATCH --mem=50MB
#SBATCH --time=0-2:00:00
#SBATCH --job-name=chopper
#SBATCH -o /nesi/nobackup/vuw03922/stdout/chopper.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/chopper.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module purge
module load chopper/0.5.0-GCC-11.3.0

#set up
PROJECT=$1
LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/fastq

#run chopper to filter low quality reads
#gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper --headcrop 50 --threads 8 | gzip > $read_dir/$LIB'_trimmed_clean_filtered_h50.fastq.gz'
gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 8 --headcrop 50 --threads 8 | gzip > $read_dir/$LIB'_trimmed_clean_filtered_q8_h50.fastq.gz' 
#gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 15 --headcrop 50 --threads 8 | gzip > $read_dir/$LIB'_trimmed_clean_filtered_q15_h50.fastq.gz'
#gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 7 --headcrop 50 --threads 8 | gzip > $read_dir/$LIB'_trimmed_clean_filtered_q7_h50.fastq.gz'
#gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 10 --threads 8 | gzip > $read_dir/$LIB'_trimmed_clean_filtered_q10.fastq.gz'
#gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 15 --threads 8 | gzip > $read_dir/$LIB'_trimmed_clean_filtered_q15.fastq.gz'

#gridion lib
#Cluster: mahuika
#Job ID: 40148582
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   42.6%  00:51:04 of 02:00:00 time limit
#CPU Efficiency: 138.5%  01:10:43 of 00:51:04 core-walltime
#Mem Efficiency:   0.6%  2.75 MB of 500.00 MB

#promethion lib
#Cluster: mahuika
#Job ID: 40746303
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   86.6%  12:59:40 of 15:00:00 time limit
#CPU Efficiency: 139.7%  18:08:55 of 12:59:40 core-walltime
#Mem Efficiency:  11.3%  5.66 MB of 50.00 MB

