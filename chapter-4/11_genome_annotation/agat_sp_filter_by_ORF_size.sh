#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=0:15:00
#SBATCH --job-name=agat_stats
#SBATCH -o /nesi/nobackup/vuw03922/stdout/agat_stats.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/agat_stats.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load AGAT/1.0.0-gimkl-2022a-Perl-5.34.1-R-4.2.1

#define variables
PROJECT=$1
SIZE=$2 #genome size in order to compute more statistics - Antip 1254625516, Gibsons 1263156490

#define paths
GFF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation/braker_extra_protein
AGAT_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation/braker_extra_protein/AGAT_stats

#filter by ORF size
agat_sp_filter_by_ORF_size.pl	--gff $GFF/braker.gff3 			\
				--size 50				\
				-o $GFF/braker_filtered_by_ORF_size.gff3


#run agat
agat_sp_functional_statistics.pl	--gff $GFF/braker_filtered_by_ORF_size.gff3 	\
					--gs $SIZE 					\
					-o $AGAT_DIR/filtered_by_ORF_size_report.txt

#Cluster: mahuika
#Job ID: 43550708
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    4.1%  00:02:26 of 01:00:00 time limit
#CPU Efficiency:  95.2%  00:02:19 of 00:02:26 core-walltime
#Mem Efficiency:  27.3%  1.36 GB of 5.00 GB

