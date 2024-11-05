#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=1G
#SBATCH --time=0:10:00
#SBATCH --job-name=flye_mito
#SBATCH -o /nesi/nobackup/vuw03922/stdout/flye_mito.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/flye_mito.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/sequence/ont
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/mitogenome/flye_mito_only/1_raw
mkdir -p $output_dir

#run program
flye	--nano-hq $read_dir/mitochondrial_reads_$FILTER.fastq.gz \
	--out-dir $output_dir/$SAMPLE'_'$FILTER'_mito' \
	--genome-size 0.019m \
	--threads 10 \
	--read-error 0.03  \
	--scaffold \
	--no-alt-contigs 

#then visualise in bandage, figure out which contig is the mitogenome and extract on command line using
#awk -v RS='>' '/contig_10/ {print ">"$0}' assembly.fasta > blue-45G_mitogenome.fasta

#Cluster: mahuika
#Job ID: 42727577
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:    2.3%  00:04:05 of 03:00:00 time limit
#CPU Efficiency:  71.2%  00:14:32 of 00:20:25 core-walltime
#Mem Efficiency:  10.9%  448.27 MB of 4.00 GB

