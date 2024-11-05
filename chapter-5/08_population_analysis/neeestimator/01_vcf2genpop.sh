#!/bin/bash
#SWATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0-1:0:00
#SBATCH --job-name=vcf2gpop
#SBATCH -o /nesi/nobackup/vuw03922/stdout/vcf2gpop.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/vcf2gpop.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load R/4.3.1-gimkl-2022a

#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
vcf2genpop=$R/vcf2genpop.R

###resources
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/${PROJECT}'_pop_info.tsv'

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
GPOP=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/genepop/$SET
mkdir -p $GPOP

#convert to gds
Rscript $vcf2genpop --gzvcf $VCF'_40chr.vcf.gz'	\
		--genpop_out $GPOP/${SET}'_neutral'		\
		--pop_file $pop_file




#Cluster: mahuika
#Job ID: 44099217
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   26.3%  00:01:19 of 00:05:00 time limit
#CPU Efficiency:  92.4%  00:01:13 of 00:01:19 core-walltime
#Mem Efficiency:  12.1%  124.09 MB of 1.00 GB

