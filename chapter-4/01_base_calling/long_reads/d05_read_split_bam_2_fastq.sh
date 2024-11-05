#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50MB
#SBATCH --time=05:00:00
#SBATCH --job-name=bam2fastq
#SBATCH -o /nesi/nobackup/vuw03922/stdout/bam2fastq.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/bam2fastq.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0

#set up
PROJECT=$1
LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/duplex_bam
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/fastq
mkdir $out_dir

#run conversion
samtools view $read_dir/$LIB.bam -d dx:1 | samtools fastq | bgzip > $out_dir/$LIB.duplex.fastq.gz
samtools view $read_dir/$LIB.bam -d dx:0 | samtools fastq | bgzip > $out_dir/$LIB.simplex.fastq.gz
samtools view $read_dir/$LIB.bam -d dx:-1 | samtools fastq | bgzip > $out_dir/$LIB.simplex_w_duplex_offspring.fastq.gz


#Cluster: mahuika
#Job ID: 40619770
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   40.7%  03:15:27 of 08:00:00 time limit
#CPU Efficiency: 113.5%  03:41:46 of 03:15:27 core-walltime
#Mem Efficiency:   1.0%  20.21 MB of 2.00 GB

