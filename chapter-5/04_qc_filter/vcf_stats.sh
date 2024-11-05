#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500M
#SBATCH --time=0-0:30:00
#SBATCH --job-name=vcfstats
#SBATCH -o /nesi/nobackup/vuw03922/stdout/vcfstats_%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/vcfstats_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load HTSlib/1.19-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET'_raw.vcf.gz'
OUT_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/vcfstats
mkdir -p $OUT_DIR

#calculate allele frequency
vcftools --gzvcf $VCF --freq2 --out $OUT_DIR/$SET --max-alleles 2

#calculate mean depth per individual
vcftools --gzvcf $VCF --depth --out $OUT_DIR/$SET

#calculate mean depth per site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT_DIR/$SET

#calculate quality score per site
vcftools --gzvcf $VCF --site-quality --out $OUT_DIR/$SET

#calculate proportion missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $OUT_DIR/$SET

#calculate proportion missing data per site
vcftools --gzvcf $VCF --missing-site --out $OUT_DIR/$SET

#calculate heteryozygosity and inbreeding coefficient per individual
#note: expected het might be overestimated if samples not all from sample pop due to Wahlund effect
vcftools --gzvcf $VCF --het --out $OUT_DIR/$SET

#Cluster: mahuika
#Job ID: 45156379
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   15.7%  00:18:52 of 02:00:00 time limit
#CPU Efficiency:  99.3%  00:18:44 of 00:18:52 core-walltime
#Mem Efficiency:   0.1%  3.80 MB of 5.00 GB

