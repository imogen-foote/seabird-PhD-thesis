#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --time=1-0:00:00
#SBATCH --job-name=racon
#SBATCH -o /nesi/nobackup/vuw03922/stdout/racon.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/racon.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load Racon/1.5.0-GCC-11.3.0
module load minimap2/2.24-GCC-11.3.0
module load HTSlib/1.18-GCC-11.3.0

#define variables
PROJECT=$1
SAMPLE=$2 #blue-45G or d36
FILTER=$3

fastq=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/sequence/ont/*combined_$FILTER.fastq.gz
bam=/nesi
assembly=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/raven/1_raw/$SAMPLE'_'$FILTER/assembly.fasta
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/raven/2_polished/$SAMPLE'_'$FILTER'_racon'
mkdir -p $out_dir
sam=$out_dir/$SAMPLE'_'$FILTER'_mapped.sam'

#map reads to assembly to get overlaps between sequences and targe
minimap2 -ax map-ont $assembly $fastq > $sam
bgzip $sam

#run racon
racon -t 10 -q 12 \
	$fastq \
	$sam'.gz' \
	$assembly > $out_dir/$SAMPLE'_'$FILTER'_racon.fasta'


#for minimap step only
#Cluster: mahuika
#Job ID: 41937658
#State: FAILED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:   27.9%  13:24:29 of 2-00:00:00 time limit
#CPU Efficiency:  47.6%  1-07:53:27 of 2-19:02:25 core-walltime
#Mem Efficiency:  24.6%  12.30 GB of 50.00 GB

#for racon step only
#Cluster: mahuika
#Job ID: 41955041
#State: COMPLETED
#Cores: 5
#Tasks: 1
#Nodes: 1
#Job Wall-time:    6.1%  02:54:29 of 2-00:00:00 time limit
#CPU Efficiency: 159.2%  23:08:41 of 14:32:25 core-walltime
#Mem Efficiency:  84.2%  210.55 GB of 250.00 GB

