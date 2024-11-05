#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G
#SBATCH --time=0-8:00:00
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
assembly_name=$2
assembler=$3

#define directories
assembly=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/$assembler/2_polished/$assembly_name/consensus.fasta
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter3/${assembler}_assembly/2_polished/$assembly_name/$assembly_name.busco

mkdir -p $out_dir

echo "Run BUSCO"
cd $out_dir

#run busco
busco 	-i $assembly \
	-o $assembly_name'_busco' \
	-l aves_odb10 \
	-m genome \
	--cpu 10 
#--restart


#manual: https://busco.ezlab.org/

#Cluster: mahuika
#Job ID: 40995492
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:   42.6%  05:06:37 of 12:00:00 time limit
#CPU Efficiency: 178.6%  1-21:37:35 of 1-01:33:05 core-walltime
#Mem Efficiency:  30.2%  15.10 GB of 50.00 GB

