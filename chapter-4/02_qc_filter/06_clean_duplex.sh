#!/bin/bash 
#SBATCH --cpus-per-task=8
#SBATCH --mem=5G
#SBATCH --time=0-12:00:00
#SBATCH --job-name=clean
#SBATCH -o /nesi/nobackup/vuw03922/stdout/clean.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/clean.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#call modules
module purge
module load Nextflow/23.04.1
module load Singularity/3.11.3


#set up
PROJECT=$1
LIB=$2
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB/duplex_fastq
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/

mkdir $read_dir/clean_work
cd $read_dir/clean_work

#run clean to remove DNA CS reads
nextflow run hoelzer/clean -r v1.0.0-beta.1 \
	--input_type nano \
	--input $read_dir/$LIB'.trimmed_duplex.fastq.gz' \
	--control dcs \
	--output $read_dir/clean_output \
	-profile singularity	\
	--cleanup_work_dir    

#--host eco,hsa 

mv $read_dir/clean_output/minimap2/$LIB/$LIB'.unmapped.fastq.gz' $read_dir/$LIB'_trimmed_clean_duplex.fastq.gz'
mv $read_dir/clean_output/Summary/multiqc_report.html $output_dir/multiQC/$LIB'_trimmed_clean_DUPLEX_multiqc.html'

