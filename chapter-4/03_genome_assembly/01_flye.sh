#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=220G
#SBATCH --time=4-00:0:00
#SBATCH --job-name=flye
#SBATCH -o /nesi/nobackup/vuw03922/stdout/flye.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/flye.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#--qos=debug

#3d, 220G 10cpu

#load modules
module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/*/fastq
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/1_raw

#run program
flye	--nano-hq $read_dir/*$FILTER.fastq.gz \
	--out-dir $output_dir/$SAMPLE'_'$FILTER \
	--genome-size 1.2g \
	--threads 10 \
	--read-error 0.03  \
	--scaffold \
	--no-alt-contigs 
#	--resume



#full dataset, headcrop 50 only
#Cluster: mahuika
#Job ID: 40888253
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:   81.1%  2-10:22:16 of 3-00:00:00 time limit
#CPU Efficiency: 183.5%  22-07:24:36 of 12-03:51:20 core-walltime
#Mem Efficiency:  88.8%  266.55 GB of 300.00 GB


#filtered q8 + headcrop 50
#Cluster: mahuika
#Job ID: 41061800
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:   43.1%  1-07:00:23 of 3-00:00:00 time limit
#CPU Efficiency: 179.0%  11-13:31:28 of 6-11:01:55 core-walltime
#Mem Efficiency:  59.7%  179.04 GB of 300.00 GB

