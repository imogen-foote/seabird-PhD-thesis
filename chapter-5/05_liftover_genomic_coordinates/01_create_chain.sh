#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=2-0:0:00
#SBATCH --job-name=create_chain
#SBATCH -o /nesi/nobackup/vuw03922/stdout/create_chain.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/create_chain.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load Nextflow/23.10.0
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

conda activate /nesi/project/vuw03922/nf-LO


#variables
antipodean=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked/blue-45G_softmasked_longest65.fasta
gibsons=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/d36_q8_h50_masked/d36_softmasked_longest65.fasta
out_dir=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter4/liftover
mkdir -p $out_dir

nextflow run evotools/nf-LO \
	--source $gibsons 	\
	--target $antipodean	\
	--distance near 	\
	--aligner lastz		\
	--outdir $out_dir/liftover_antipodean_gibsons_65scaff \
	-profile conda		\
	--max_cpus 2		\
	--max_memory 6.GB 

conda deactivate


#Cluster: mahuika
#Job ID: 43680580
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:   40.6%  19:28:52 of 2-00:00:00 time limit
#CPU Efficiency:  62.5%  2-00:40:11 of 3-05:55:28 core-walltime
#Mem Efficiency:  21.6%  3.45 GB of 16.00 GB

