#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=assembly-stats
#SBATCH -o /nfs/scratch/footeim/stdout/assembly-stats.%j.out
#SBATCH -e /nfs/scratch/footeim/stdout/assembly-stats.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#activate conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate assembly-stats

#set environment
PROJECT=$1
ASSEMBLY=$2
ASSEMBLER=$3

#set paths
read_dir=$SCRATCH/projects/$PROJECT/data/Chapter3/$ASSEMBLER/$ASSEMBLY
output_dir=$SCRATCH/projects/$PROJECT/output/Chapter3/${ASSEMBLER}_assembly/$ASSEMBLY/assembly-stats
mkdir $SCRATCH/projects/$PROJECT/output/Chapter3/${ASSEMBLER}_assembly/$ASSEMBLY/assembly-stats

#run program
assembly-stats $read_dir/*ssembly.fasta > $output_dir/${ASSEMBLY}_assembly-stats.txt

conda deactivate
