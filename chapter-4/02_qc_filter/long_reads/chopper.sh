#!/bin/bash 
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-00:30:00
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

##run chopper to filter low quality reads
gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 10 | gzip > $read_dir/$LIB'_choppertest_q10.fastq.gz'
gunzip -c $read_dir/$LIB'_trimmed_clean.fastq.gz' | chopper -q 5 | gzip > $read_dir/$LIB'_choppertest_q5.fastq.gz'


