#!/bin/bash -e
#SBATCH --gpus-per-node=A100:1
#SBATCH --partition=gpu,hgx
#SBATCH --cpus-per-task=4
#SBATCH --mem=11G
#SBATCH --time=0-15:00:00 
#SBATCH --job-name=dorado
#SBATCH -o /nesi/nobackup/vuw03922/stdout/dorado.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/dorado.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module purge
module load Dorado/0.3.4

#set up
PROJECT=$1
LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/pod5
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/bam

#run dorado basecaller
dorado basecaller --device 'cuda:all' dna_r10.4.1_e8.2_400bps_sup@v4.1.0 $read_dir/ > $out_dir/$LIB.calls.bam

#resume interrupted dorado run
#dorado basecaller --device 'cuda:all' dna_r10.4.1_e8.2_400bps_sup@v4.1.0 $read_dir/ --resume-from $out_dir/$LIB.calls.bam > $out_dir/calls.bam

#remove incomplete bam file (from interrupted run)
#rm $out_dir/$LIB.calls.bam
#mv $out_dir/calls.bam $out_dir/$LIB.calls.bam

#create sequencing summary file
dorado summary $out_dir/$LIB.calls.bam > $out_dir/sequencing_summary_$LIB.txt


#first run
#Cluster: mahuika
#Job ID: 39221238
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   99.2%  09:55:26 of 10:00:00 time limit
#CPU Efficiency:  66.1%  13:06:56 of 19:50:52 core-walltime
#Mem Efficiency:  86.4%  8.64 GB of 10.00 GB

#resumed run
#Cluster: mahuika
#Job ID: 39650463
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   19.4%  01:09:44 of 06:00:00 time limit #some libraries took extra 4 h
#CPU Efficiency:  79.6%  01:51:05 of 02:19:28 core-walltime
#Mem Efficiency:  78.1%  8.59 GB of 11.00 GB

