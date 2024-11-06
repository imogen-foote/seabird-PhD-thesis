#!/bin/bash 
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-00:5:00
#SBATCH --job-name=pycoqc
#SBATCH -o /nesi/nobackup/vuw03922/stdout/pycoqc.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/pycoqc.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load Python/3.10.5-gimkl-2022a

#activate conda environment
source /nesi/project/vuw03922/ont_env/bin/activate
export PYTHONNOUSERSITE=1 #makes sure local packages installed in home folder ~/.local/lib/pythonX.Y/site-packages/ (where X.Y is the Python version, e.g. 3.8) by pip install --user are excluded from your conda environments.

#set environment
PROJECT=$1
LIB=$2

#set paths
#for guppy libraries
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/fastq
#for dorado libraries
#read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/bam

output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/pycoQC/

#run program
pycoQC -f $read_dir/sequencing_summary*.txt --min_pass_qual 8 -o $output_dir/$LIB.pycoqc_qual8.html
#pycoQC -f $read_dir/sequencing_summary*.txt --min_pass_qual 10 --min_pass_len 200 -o $output_dir/$LIB.pycoqc_qual10_len200.html
#pycoQC -f $read_dir/sequencing_summary*.txt --min_pass_qual 15 --min_pass_len 200 -o $output_dir/$LIB.pycoqc_qual15_len200.html



deactivate




#FOR GRIDION/GUPPY LIBRARIES
#Cluster: mahuika
#Job ID: 39437462
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    1.1%  00:00:41 of 01:00:00 time limit
#CPU Efficiency:  43.9%  00:00:18 of 00:00:41 core-walltime
#Mem Efficiency:  12.4%  253.12 MB of 2.00 GB

#FOR PROMETHION/DORADO LIBRARIES
#Cluster: mahuika
#Job ID: 39439002
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    9.5%  00:00:57 of 00:10:00 time limit
#CPU Efficiency:  89.5%  00:00:51 of 00:00:57 core-walltime
#Mem Efficiency:  34.5%  1.04 GB of 3.00 GB

