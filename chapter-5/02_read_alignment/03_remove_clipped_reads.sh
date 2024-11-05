#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH -a 1-43
#SBATCH --mem=30G
#SBATCH --time=1:30:00
#SBATCH --job-name=rm_clipped
#SBATCH -o /nesi/nobackup/vuw03922/stdout/rm_clipped.%A_%a.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/rm_clipped.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
 
PROJECT=$1
IND=$2 #blue-45G or d36
N=${SLURM_ARRAY_TASK_ID}


sample_list=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/paleomix/samples.txt
sample=$( head -n $N $sample_list | tail -n 1 | cut -f 1 )

read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/paleomix/${sample}


# create list of all clipped reads
samtools view $read_dir/$sample'_1.0.'$IND'_nuclear.bam' |awk '$6 ~ /H|S/ {print $1}' |sort -u > $read_dir/$sample'_clipped_reads_list.txt'

# filter clipped reads
samtools view -h $read_dir/$sample'_1.0.'$IND'_nuclear.bam' | fgrep -wvf $read_dir/$sample'_clipped_reads_list.txt' | samtools view -b - -o $read_dir/$sample'_1.0.'$IND'_nuclear_clipped.bam'

# create index for new bam file
samtools index $read_dir/$sample'_1.0.'$IND'_nuclear_clipped.bam'


#for row in $(seq 1 43); do
#	#extract sample name for current iteration
#	sample=$(head -n $row $sample_list | tail -n 1 | cut -f 1)
#
#	#define path to read dir for the current sample
#	read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/paleomix/${sample}
#
#	#create list of all clipped reads
#	samtools view $sample'_1.0_nuclear.bam' |awk '$6 ~ /H|S/ {print $1}' |sort -u > $sample'_clipped_reads_list.txt'
#
#	#filter clipped reads
#	samtools view -h $sample'_1.0_nuclear.bam' | fgrep -wvf $sample'_clipped_reads_list.txt' | samtools view -b - -o $sample'_1.0_nuclear_clipped.bam'
#
#	#create index for new bam file
#	samtools index $sample'_1.0_nuclear_clipped.bam'
#done
