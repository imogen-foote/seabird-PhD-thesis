#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=0-3:00
#SBATCH --job-name=multiqc
#SBATCH -o /nesi/nobackup/vuw03922/stdout/multiqc.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/fmultiqc.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module load MultiQC/1.13-gimkl-2022a-Python-3.10.5

PROJECT=$1

#set environment
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/fastqc/fastqc_raw
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/multiqc/multiqc_raw
mkdir -p $out_dir

#run program
multiqc $read_dir/ --outdir $out_dir

