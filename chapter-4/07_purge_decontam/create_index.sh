#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=0-1:00:00
#SBATCH --job-name=index
#SBATCH -o /nesi/nobackup/vuw03922/stdout/index.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/index.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0

PROJECT=$1
ASSEMBLER=$2
ASSEMBLY=$3

read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/$ASSEMBLER/2_polished/$ASSEMBLY

#index
samtools faidx $read_dir/consensus.fasta
