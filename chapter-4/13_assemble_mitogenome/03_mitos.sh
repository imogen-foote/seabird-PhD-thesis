#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --job-name=mitos
#SBATCH -o /nesi/nobackup/vuw03922/stdout/mitos.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/mitos.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load conda env
module purge && module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1
conda activate /nesi/project/vuw03922/mitos


#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/mitogenome/flye_mito_only/$SAMPLE'_'$FILTER'_mito'
ref=/nesi/nobackup/vuw03922/resources/metazoan_mitos_ref
#mitos=/nesi/nobackup/vuw03922/scripts/Chapter3/extract_mitogenome/mitos2_wrapper.py
out_dir=$read_dir/annotated
mkdir $out_dir
cd $out_dir

#run program
#$mitos FASTA $read_dir/*.fasta


runmitos.py -i $read_dir/$SAMPLE'_mitogenome.fasta' \
	-o $out_dir \
	-c 2 \
	-R $ref \
	-r refseq63m

conda deactivate

