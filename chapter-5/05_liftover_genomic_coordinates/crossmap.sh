#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:05:00
#SBATCH --job-name=crossmap
#SBATCH -o /nesi/nobackup/vuw03922/stdout/crossmap.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/crossmap.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

conda activate /nesi/project/vuw03922/crossmap

#define paths
vcf=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter4/snp/GibsonsAlbatross_full
chain=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter4/liftover/liftover_antipodean_gibsons/chainnet/liftover.chain
reference=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked/blue-45G_softmasked_renamed.fasta
out_dir=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter4/crossmap_snp
mkdir -p $out_dir

#make tmp copy of vcf file and unzip 
mkdir $out_dir/tmp
cp $vcf/GibsonsAlbatross_full_qc.vcf.gz $out_dir/tmp/tmp_GibsonsAlbatross_full_qc.vcf.gz
gzip -d $out_dir/tmp/tmp_GibsonsAlbatross_full_qc.vcf.gz

#run crossmap
CrossMap.py vcf $chain	\
	$out_dir/tmp/tmp_GibsonsAlbatross_full_qc.vcf		\
	$reference					\
	$out_dir/GibsonsAlbatross_full_qc_crossmap.vcf	

# CrossMap.py vcf input_chain_file input_VCF_file ref_genome_file output_file

#remove tmp directory
rm -r $out_dir/tmp

conda deactivate
