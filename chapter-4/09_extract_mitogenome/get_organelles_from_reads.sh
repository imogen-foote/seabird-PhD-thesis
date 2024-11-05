#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --time=05:00:00
#SBATCH --job-name=getorganelles
#SBATCH -o /nesi/nobackup/vuw03922/stdout/getorganelles.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/getorganelles.%j.err
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
R1=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/sequence/illumina/R1/*_combined_R1.fastq.gz
R2=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/sequence/illumina/R2/*_combined_R2.fastq.gz
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/mitogenome/$SAMPLE'_'$FILTER'_from_reads'
mkdir -p $output_dir

#download organelle
get_organelle_config.py -a animal_mt

#run program
get_organelle_from_reads.py -1 $R1 -2 $R2 -o $output_dir -F animal_mt -t 10 --overwrite

conda deactivate

