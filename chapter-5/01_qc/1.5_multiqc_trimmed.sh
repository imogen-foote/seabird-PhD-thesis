#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=multiqc_trimmed
#SBATCH -o /nfs/scratch2/footeim/Ant_albatross/stdout/multiqc_trimmed.%j.out
#SBATCH -e /nfs/scratch2/footeim/Ant_albatross/stdout/multiqc_trimmed.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module load multiqc/1.7

#set environment
read_dir=/nfs/scratch2/footeim/Ant_albatross/output/fastqc_trimmed
out_dir=/nfs/scratch2/footeim/Ant_albatross/output/multiqc_trimmed

#run program
multiqc $read_dir/ --outdir $out_dir

