#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=50MB
#SBATCH --time=0:10:00
#SBATCH --job-name=concatfastq
#SBATCH -o /nesi/nobackup/vuw03922/stdout/concatfastq.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/concatfastq.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#set up
PROJECT=$1
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont

mkdir $read_dir/all_fastq
mkdir $read_dir/all_fastq/simplex_tmp
mkdir $read_dir/all_fastq/duplex_tmp

cp $read_dir/*/fastq/*simplex.fastq.gz $read_dir/all_fastq/simplex_tmp
cat $read_dir/all_fastq/simplex_tmp/*.fastq.gz > $read_dir/all_fastq/$PROJECT.simplex.fastq.gz

cp $read_dir/*/fastq/*duplex.fastq.gz $read_dir/all_fastq/duplex_tmp
cat $read_dir/all_fastq/duplex_tmp/*.fastq.gz > $read_dir/all_fastq/$PROJECT.duplex.fastq.gz

#rm -r $read_dir/all_fastq/simplex_tmp
#rm -r $read_dir/all_fastq/duplex_tmp 
