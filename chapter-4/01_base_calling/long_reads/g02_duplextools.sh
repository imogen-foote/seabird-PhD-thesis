#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G
#SBATCH --time=0-1:00
#SBATCH --job-name=duplextools
#SBATCH -o /nesi/nobackup/vuw03922/stdout/duplextools.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/duplextools.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load Python/3.10.5-gimkl-2022a
module load duplex-tools/0.2.20-gimkl-2022a-Python-3.10.5

#set environment
PROJECT=$1
LIB=$2

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/duplex_tools

#run program
#investigate sequencing summary output from Guppy to identify candidate pairs
duplex_tools pairs_from_summary $read_dir/fastq/sequencing_summary*.txt $output_dir/$LIB
#'basecall-to-basecall alignment filtering' - analyses basecalls of candidate reads, checking for similarity
duplex_tools filter_pairs $output_dir/$LIB/pair_ids.txt $read_dir/fastq
