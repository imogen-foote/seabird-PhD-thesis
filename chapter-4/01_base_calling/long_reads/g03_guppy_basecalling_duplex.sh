#!/bin/bash
#SBATCH --gpus-per-node=A100:1
#SBATCH --partition=gpu,hgx
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=0-15:00:00
#SBATCH --job-name=guppy_duplex
#SBATCH -o /nesi/nobackup/vuw03922/stdout/guppy_duplex.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/guppy_duplex.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load ont-guppy-gpu/6.5.7
module load SAMtools/1.16.1-GCC-11.3.0 

#set up
PROJECT=$1
LIB=$2
MODEL=$3 #either dna_r10.4_e8.1_sup.cfg or dna_r9.4.1_e8.1_sup.cfg

pair_ids=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/duplex_tools/$LIB/pair_ids_filtered.txt
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB

mkdir $read_dir/duplex_calls/TMP_fast5
cp $read_dir/fast5/fast5_*/*.fast5 $read_dir/duplex_fastq/TMP_fast5

mkdir $read_dir/duplex_fastq/TMP_fastq

### run guppy
guppy_basecaller_duplex 			\
    -i $read_dir/duplex_fastq/TMP_fast5		\
    -r 					\
	-s $read_dir/duplex_fastq/TMP_fastq	\
    --device auto 			\
	-c /opt/nesi/CS400_centos7_bdw/ont-guppy-gpu/6.5.7/data/$MODEL	\
    --duplex_pairing_mode from_pair_list	\
    --duplex_pairing_file $pair_ids	\
	--num_callers 4

#concatenate reads
cat $read_dir/duplex_fastq/TMP_fastq/pass/*fastq > $read_dir/duplex_fastq/$LIB.duplex.pass.fastq
bgzip $read_dir/duplex_fastq/$LIB.duplex.pass.fastq
cat $read_dir/duplex_fastq/TMP_fastq/fail/*fastq > $read_dir/duplex_fastq/$LIB.duplex.fail.fastq
bgzip $read_dir/duplex_fastq/$LIB.duplex.fail.fastq
cat $read_dir/duplex_fastq/$LIB.duplex.pass.fastq.gz $read_dir/duplex_fastq/$LIB.duplex.fail.fastq.gz > $read_dir/duplex_fastq/$LIB.duplex.fastq.gz


#rename summary files with duplex 
mv $read_dir/fastq/TMP_fastq/sequencing_summary.txt $read_dir/fastq/TMP_fastq/sequencing_summary_$LIB_duplex.txt
mv $read_dir/fastq/TMP_fastq/sequencing_telemetry.js $read_dir/fastq/TMP_fastq/sequencing_telemetry_$LIB_duplex.js

#copy files from tmp_fastq to fastq directory
cp $read_dir/duplex_fastq/TMP_fastq/* $read_dir/duplex_fastq
#then remove individual fastq files and tmp directories
rm $read_dir/duplex_fastq/*.fastq
rm -r $read_dir/duplex_fastq/TMP_fast5
rm -r $read_dir/duplex_fastq/TMP_fastq


#Cluster: mahuika
#Job ID: 39534227
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   13.5%  03:15:04 of 1-00:00:00 time limit
#CPU Efficiency:  52.3%  03:23:56 of 06:30:08 core-walltime
#Mem Efficiency:   6.9%  1.37 GB of 20.00 GB


#Cluster: mahuika
#Job ID: 39609267
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   48.4%  11:37:04 of 1-00:00:00 time limit
#CPU Efficiency:  98.7%  22:56:21 of 23:14:08 core-walltime
#Mem Efficiency:  23.1%  1.84 GB of 8.00 GB
