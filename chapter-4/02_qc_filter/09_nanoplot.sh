#!/bin/bash 
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-1:00:00
#SBATCH --job-name=nanoplot
#SBATCH -o /nesi/nobackup/vuw03922/stdout/nanoplot.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/nanoplot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load NanoPlot/1.41.0-gimkl-2022a-Python-3.10.5

#set environment
PROJECT=$1
LIB=$2

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/fastq
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/nanoplot/

#run program
#NanoPlot -t 2 --fastq $read_dir/$LIB'.fastq.gz' -o $output_dir/$LIB
#NanoPlot -t 2 --fastq $read_dir/$LIB'.trimmed.fastq.gz' -o $output_dir/$LIB'_trimmed'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean.fastq.gz' -o $output_dir/$LIB'_trimmed_clean'
NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q8_h50.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q8_h50' 
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q10.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q10'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q15.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q15'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q7_h50.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q7_h50'
