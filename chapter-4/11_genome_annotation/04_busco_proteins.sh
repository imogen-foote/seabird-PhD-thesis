#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=1500M
#SBATCH --time=0-0:30:00
#SBATCH --job-name=BUSCO
#SBATCH -o /nesi/nobackup/vuw03922/stdout/BUSCO.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/BUSCO.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load BUSCO/5.3.2-gimkl-2020a

#set variables
PROJECT=$1
SAMPLE=$2
FILTER=$3

#define directories
proteins=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation/braker_extra_protein/braker.aa
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter3/flye_assembly/6_gene_annotation/$SAMPLE'_'$FILTER/$SAMPLE'_'$FILTER'_genes_extra_protein.busco'

mkdir -p $out_dir

echo "Run BUSCO"
cd $out_dir

#run busco
busco 	-i $proteins \
	-o $SAMPLE'_'$FILTER'_genes_extra_protein_busco' \
	-l aves_odb10 \
	-m proteins \
	--cpu 10 
#	--restart


#manual: https://busco.ezlab.org/

#Cluster: mahuika
#Job ID: 43646580
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:    3.2%  00:19:12 of 10:00:00 time limit
#CPU Efficiency: 141.0%  02:15:22 of 01:36:00 core-walltime
#Mem Efficiency:   3.1%  946.18 MB of 30.00 GB
