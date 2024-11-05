#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=3G
#SBATCH --partition=quicktest
#SBATCH --time=0-1:00
#SBATCH --job-name=fastqc_trimmed
#SBATCH -o /nfs/scratch2/footeim/Ant_albatross/stdout/fastqc_trimmed.%j.out
#SBATCH -e /nfs/scratch2/footeim/Ant_albatross/stdout/fastqc_trimmed.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module load fastqc/0.11.7

#set environment
sample=$1
read_dir=/nfs/scratch2/footeim/Ant_albatross/output/trimmed	
output_dir=/nfs/scratch2/footeim/Ant_albatross/output/fastqc_trimmed

#run program
fastqc $read_dir/$sample.*_p*.fq.gz  -t 2 --noextract --outdir $output_dir/

#mv $read_dir/$sample*_fastqc.zip  $output_dir
#mv $read_dir/$sample*_fastqc.html $output_dir
