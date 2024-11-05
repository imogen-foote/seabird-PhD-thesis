#!/bin/bash
#SWATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:10:00
#SBATCH --job-name=remove_miss
#SBATCH -o /nesi/nobackup/vuw03922/stdout/remove_miss.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/remove_miss.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load HTSlib/1.13-GCCcore-9.2.0

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET

#filter
vcftools --gzvcf $VCF'_neutral.vcf.gz'		\
	--out $VCF'_neutral_nomiss'	\
	--max-missing 0.99			\
	--recode-INFO-all	--recode
mv $VCF'_neutral_nomiss.recode.vcf' $VCF'_neutral_nomiss.vcf'


bgzip $VCF'_neutral_nomiss.vcf'
