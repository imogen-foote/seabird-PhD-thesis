#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=3G
#SBATCH --time=0-0:20:00
#SBATCH --job-name=merge_vcf
#SBATCH -o /nesi/nobackup/vuw03922/stdout/merge_vcf.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/merge_vcf.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#define variables
PROJECT=$1
SET=$PROJECT'_'$2

#load modules
module purge
module load BCFtools/1.19-GCC-11.3.0
module load R/4.3.2-foss-2023a

antipodean=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter4/snp/AntipodeanAlbatross_full/AntipodeanAlbatross_full_qc.vcf.gz
gibsons=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter4/liftover_snp/GibsonsAlbatross_liftover_full_qc.vcf.gz
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET
mkdir -p $out_dir

#tab delimited population file
POP=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/${SET}_pop_info.tsv

#custom scripts
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
vcf2R=$R/vcf2Rinput.R
gds2plink=$R/gds2plink.R


index vcf file
bcftools index $gibsons

#merge
bcftools merge $antipodean $gibsons		\
	-o $out_dir/${SET}_ALLsnps_qc.vcf.gz

#identify overlapping snps
bcftools isec -n=2 -Oz					\
	-o $out_dir/${SET}_overlapping_snps.txt	\
	$antipodean $gibsons

#filter to only retain overlapping snps
bcftools view -T $out_dir/${SET}_overlapping_snps.txt -Oz 	\
	-o $out_dir/${SET}_overlapping_qc.vcf.gz		\
	$out_dir/${SET}_ALLsnps_qc.vcf.gz

#create gds
Rscript $vcf2R --gzvcf $out_dir/${SET}_qc.vcf.gz	\
	--snprelate_out $out_dir/${SET}_qc			\
	--genlight_out $out_dir/${SET}_qc

#create plink bed file
Rscript $gds2plink --gds_file $out_dir/${SET}_qc.gds  \
                        --out $out_dir/${SET}_qc	\
                        --pop_file $POP



#Cluster: mahuika
#Job ID: 44064159
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   24.4%  00:14:38 of 01:00:00 time limit
#CPU Efficiency:  99.3%  00:14:32 of 00:14:38 core-walltime
#Mem Efficiency:  34.1%  1.71 GB of 5.00 GB

