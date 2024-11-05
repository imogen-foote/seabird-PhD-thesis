#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=32
#SBATCH --mem=30G
#SBATCH --time=1-0:00:00
#SBATCH --job-name=BRAKER
#SBATCH -o /nesi/nobackup/vuw03922/stdout/BRAKER.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/BRAKER.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#3days, 32G, 8cpus

#load modules
module purge
module load Apptainer/1.2.2

#define variables
PROJECT=$1
SAMPLE=$2

export BRAKER_SIF=/nesi/nobackup/vuw03922/scripts/Chapter3/11_genome_annotation/braker3.sif

protein=/nesi/nobackup/vuw03922/resources/protein_db/combined_protein.fasta
genome=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/5_repeat_masking/$SAMPLE'_q8_h50_masked'/$SAMPLE'_softmasked_renamed.fasta'
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation/braker_extra_protein
mkdir $out_dir
cd $out_dir

#run braker
singularity exec $BRAKER_SIF braker.pl --genome=$genome --prot_seq=$protein --threads=32 --gff3 --species=$PROJECT'_extra' 


#Cluster: mahuika
#Job ID: 43368254
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:   79.0%  3-03:52:21 of 4-00:00:00 time limit
#CPU Efficiency: 167.3%  21-03:46:45 of 12-15:29:24 core-walltime
#Mem Efficiency:  45.8%  22.88 GB of 50.00 GB

