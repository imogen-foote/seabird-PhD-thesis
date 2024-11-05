#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=gaas
#SBATCH -o /nfs/scratch/footeim/stdout/gaas.%j.out
#SBATCH -e /nfs/scratch/footeim/stdout/gaas.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#activate conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda activate gaas

#set environment
PROJECT=$1
ASSEMBLY=$2
ASSEMBLER=$3

#set paths
read_dir=$SCRATCH/projects/$PROJECT/data/Chapter3/$ASSEMBLER/$ASSEMBLY
output_dir=$SCRATCH/projects/$PROJECT/output/Chapter3/${ASSEMBLER}_assembly/$ASSEMBLY/gaas_fasta_statistics
mkdir $SCRATCH/projects/$PROJECT/output/Chapter3/${ASSEMBLER}_assembly/$ASSEMBLY/gaas_fasta_statistics

#run program
gaas_fasta_statistics.pl -f $read_dir/*ssembly.fasta > $output_dir/gaas_fasta_statistics.txt

conda deactivate
