#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=20MB
#SBATCH --time=0-0:05:00
#SBATCH --job-name=make_yaml
#SBATCH -o /nesi/nobackup/vuw03922/stdout/make_yaml.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/make_yaml.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#set environment
PROJECT=$1
POP=$2
#read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/illumina
read_dir=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/raw_data/illumina
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/paleomix
mkdir -p $output_dir

for f in $read_dir/*_R1.fastq.gz ; do
base=$(basename ${f%_H*})
mkdir $output_dir/$base
cp /nesi/nobackup/vuw03922/scripts/Chapter4/02_read_alignment/blank_$POP.yaml $output_dir/$base/$base.yaml
sed -i "s/xxsamplexx/${base}/g" $output_dir/$base/$base.yaml
done 

#manually edited ANT_OR01, ANT_OR04 and ANT_SA52 as these had multiple libraries, plus gibsons ones


#Cluster: mahuika
#Job ID: 43198657
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    0.0%  00:00:01 of 01:00:00 time limit
#CPU Efficiency:   0.0%  00:00:00 of 00:00:01 core-walltime
#Mem Efficiency:   0.0%  0.00 MB of 10.00 GB
