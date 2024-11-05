#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --job-name=bandage
#SBATCH -o /nesi/nobackup/vuw03922/stdout/bandage.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/bandage.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load conda env
module purge && module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/vuw03922/bandage


#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/1_raw/$SAMPLE'_'$FILTER
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/mitogenome/$SAMPLE'_'$FILTER
mkdir -p $output_dir

#download organelle
get_organelle_config.py -a animal_mt

#run program
get_organelle_from_assembly.py 	-F animal_mt \
				-g $read_dir/assembly_graph.gfa \
				-o $output_dir \
				-t 10 \
				--overwrite \
				--verbose \
				--max-depth 48


conda deactivate

#antipodean didn't extract just one fragment so had to save the slimmed assembly graph in bandage, identify the mitochondrial (circular) contig and save that as a separate fasta file
#both antipodean and gibson's assembled as 38kb fragment so mitogenome needs to be properly assembled separately, but following scripts remove this 38km misassembled fragment

#Cluster: mahuika
#Job ID: 42115809
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   19.9%  00:59:48 of 05:00:00 time limit
#CPU Efficiency: 184.1%  01:50:04 of 00:59:48 core-walltime
#Mem Efficiency:  59.3%  2.96 GB of 5.00 GB

