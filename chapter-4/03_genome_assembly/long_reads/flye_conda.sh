#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:10:00
#SBATCH --job-name=flye
#SBATCH -o /nesi/nobackup/vuw03922/stdout/flye.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/flye.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#activate conda environment
module purge
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
source activate /nesi/project/vuw03922/flye

#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/*/fastq
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/flye

#run program
flye	--nano-hq $read_dir/*$FILTER.fastq.gz \
	--out-dir $output_dir/$SAMPLE'_'$FILTER \
	--genome-size 1.2g \
	--threads 10 \
	--read-error 0.03  \
	--scaffold \
	--no-alt-contigs

conda deactivate

#$(date +%d_%m_%Y)
