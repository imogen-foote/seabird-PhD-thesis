#!/bin/bash 
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0-2:0:00
#SBATCH --job-name=nanoplot
#SBATCH -o /nesi/nobackup/vuw03922/stdout/nanoplot.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/nanoplot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load Python/3.10.5-gimkl-2022a

#activate conda environment
source /nesi/project/vuw03922/nanoplot/bin/activate
export PYTHONNOUSERSITE=1 #makes sure local packages installed in home folder ~/.local/lib/pythonX.Y/site-packages/ (where X.Y is the Python version, e.g. 3.8) by pip install --user are excluded from your conda environments.

#set environment
PROJECT=$1
LIB=$2
#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/fastq
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/nanoplot/

#run program
#NanoPlot -t 2 --fastq $read_dir/$LIB'.trimmed.fastq.gz' -o $output_dir/$LIB'_trimmed'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean.fastq.gz' -o $output_dir/$LIB'_trimmed_clean'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_h50.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_h50'
NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q8_h50.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q8_h50'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q10.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q10'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q15.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q15'
#NanoPlot -t 2 --fastq $read_dir/$LIB'_trimmed_clean_filtered_q7_h50.fastq.gz' -o $output_dir/$LIB'_trimmed_clean_filtered_q7_h50'


deactivate

#promethion
#Cluster: mahuika
#Job ID: 40429989
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   44.1%  05:17:50 of 12:00:00 time limit
#CPU Efficiency:  99.3%  05:15:46 of 05:17:50 core-walltime
#Mem Efficiency:  52.7%  1.05 GB of 2.00 GB

#minion
#Cluster: mahuika
#Job ID: 40737467
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    8.4%  00:40:12 of 08:00:00 time limit
#CPU Efficiency:  92.7%  00:37:17 of 00:40:12 core-walltime
#Mem Efficiency:  50.2%  1.00 GB of 2.00 GB

