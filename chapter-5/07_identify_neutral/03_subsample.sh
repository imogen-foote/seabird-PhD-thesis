#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:5:00
#SBATCH --job-name=subsample_vcf
#SBATCH -o /nesi/nobackup/vuw03922/stdout/subsample_vcf.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/subsample_vcf.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load vcflib/1.0.1-GCC-9.2.0
module load BCFtools/1.13-GCC-9.2.0
module load HTSlib/1.19-GCC-11.3.0

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET

#randomly select ~ 11,000 snps from full dataset (11,000 / 60488 - 0.18 for sampling rate)
bcftools view $VCF'_neutral.vcf.gz' | vcfrandomsample -r 0.3 > $VCF'_neutral_subset.vcf'

bgzip $VCF'_neutral_subset.vcf'
