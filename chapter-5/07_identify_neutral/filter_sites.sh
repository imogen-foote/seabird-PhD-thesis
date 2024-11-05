#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:05:00
#SBATCH --job-name=filter_sites
#SBATCH -o /nesi/nobackup/vuw03922/stdout/filter_sites.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/filter_sites.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load BCFtools/1.13-GCC-9.2.0
module load R/4.3.1-gimkl-2022a

#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
vcf2R=$R/vcf2Rinput.R

#set variables
PROJECT=$1
SET=${PROJECT}_$2
FILTER=$3

POP1="GIB"
POP2="DI"

#define paths
dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/${SET}/${SET}_${FILTER}
vcf_in=${dir}.vcf.gz
vcf_out=${dir}_${POP1}_${POP2}
SAMPLES=/nesi/nobackup/vuw03922/projects/AllAlbatross/resources/sample_info/${POP1}_${POP2}.txt

vcftools --gzvcf $vcf_in		\
	--keep $SAMPLES			\
	--recode-INFO-all --recode	\
	--out $vcf_out
mv ${vcf_out}.recode.vcf ${vcf_out}.vcf
bgzip ${vcf_out}.vcf
tabix vcf ${vcf_out}.vcf.gz

#convert to gds
Rscript $vcf2R --gzvcf ${vcf_out}.vcf.gz    \
                --snprelate_out ${dir}_${POP1}_${POP2}

#Cluster: mahuika
#Job ID: 45716516
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    6.0%  00:00:54 of 00:15:00 time limit
#CPU Efficiency:  87.0%  00:00:47 of 00:00:54 core-walltime
#Mem Efficiency:   1.0%  10.08 MB of 1.00 GB

