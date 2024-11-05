#!/bin/bash
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=6G
#SBATCH --partition=bigmem
#SBATCH --time=2-1:00
#SBATCH --job-name=trimmomatic
#SBATCH -o /nfs/scratch/footeim/stdout/trim.%j.out
#SBATCH -e /nfs/scratch/footeim/stdout/trim.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
trimmomatic=/nfs/home/footeim/bin/trimmomatic-0.39.jar

#set environment
PROJECT=$1
read_dir=$SCRATCH/projects/$PROJECT/raw_data/illumina
output_dir=$SCRATCH/projects/$PROJECT/data/trimmomatic

#run trimmomatic
#trimmomatic can only run on one file at a time so need to create a for loop
for f in $read_dir/*_R1.fastq.gz ; do
base=$(basename ${f%_H*})
java -jar $trimmomatic PE \
-threads 12 -phred33 \
${f} $read_dir/${base}*_R2.fastq.gz \
$output_dir/${base}.forward_paired.fq.gz $output_dir/${base}.forward_unpaired.fq.gz \
$output_dir/${base}.reverse_paired.fq.gz $output_dir/${base}.reverse_unpaired.fq.gz \
ILLUMINACLIP:$SCRATCH/projects/$PROJECT/resources/adapt_seq/NexteraPE-PE.fa:2:30:10
done

