#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=0-0:5:00
#SBATCH --job-name=align
#SBATCH -o /nesi/nobackup/vuw03922/stdout/align.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/align.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load MAFFT/7.505-gimkl-2022a-with-extensions

#set environment
fasta=/nesi/nobackup/vuw03922/projects/BycatchID/data/fasta
outgroup=/nesi/nobackup/vuw03922/resources/reference_genomes/Phoebastria_albatrus_mitogenome/mitochondrial_sequence.fasta

#combine all sequences
cat $fasta/combined.fasta $outgroup  > $fasta/combined_w_outgroup.fasta

#align
mafft --auto $fasta/combined_w_outgroup.fasta > $fasta/combined_mafft_mitogenome_w_outgroup.fasta


#Cluster: mahuika
#Job ID: 43967788
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    2.0%  00:00:06 of 00:05:00 time limit
#CPU Efficiency:  83.3%  00:00:05 of 00:00:06 core-walltime
#Mem Efficiency:   0.0%  0.00 MB of 1.00 GB

