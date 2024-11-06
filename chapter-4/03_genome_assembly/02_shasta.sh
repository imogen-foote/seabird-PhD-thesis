#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=64
#SBATCH --mem=300G
#SBATCH --time=3-0:0:00
#SBATCH --job-name=shasta
#SBATCH -o /nesi/nobackup/vuw03922/stdout/shasta.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/shasta.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#--qos=debug

#activate conda environment
module purge && module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/vuw03922/shasta

#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/*/fastq
output_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/shasta/1_raw

#create tmp directory to copy fastq.gz files into to unzip 
tmp=$output_dir/fastq
mkdir $tmp
cp $read_dir/*$FILTER.fastq.gz $tmp

#unzip
module load HTSlib/1.18-GCC-11.3.0
parallel bgzip -d {} ::: $tmp/*.fastq.gz

#run program
echo "beginning shasta assembly"
shasta --config Nanopore-R10-Slow-Nov2022 \
	--threads 64 \
	--input $tmp/*.fastq \
	--assemblyDirectory $output_dir/$SAMPLE'_'$FILTER \
	--Reads.minReadLength 100

rm -r $tmp

conda deactivate
