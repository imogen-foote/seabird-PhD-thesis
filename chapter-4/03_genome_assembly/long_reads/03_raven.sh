#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=20
#SBATCH --mem=70G
#SBATCH --time=0-15:0:00
#SBATCH --job-name=raven
#SBATCH -o /nesi/nobackup/vuw03922/stdout/raven.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/raven.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load Raven/1.5.0-GCC-9.2.0

#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/sequence/ont
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/raven/1_raw

mkdir -p $output_dir/$SAMPLE'_'$FILTER
cd $output_dir/$SAMPLE'_'$FILTER

#run program
raven -t 20 $read_dir/*'_combined_'$FILTER.fastq.gz > $output_dir/$SAMPLE'_'$FILTER/assembly.fasta


#Cluster: mahuika
#Job ID: 41873267
#State: COMPLETED
#Cores: 10
#Tasks: 1
#Nodes: 1
#Job Wall-time:   96.6%  14:29:07 of 15:00:00 time limit
#CPU Efficiency: 160.9%  9-17:02:26 of 6-00:51:10 core-walltime
#Mem Efficiency:  75.3%  52.72 GB of 70.00 GB

