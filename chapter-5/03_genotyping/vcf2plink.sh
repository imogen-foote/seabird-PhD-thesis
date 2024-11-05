#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=quicktest
#SBATCH --time=0-5:00
#SBATCH --job-name=vcf2plink
#SBATCH -o /nfs/scratch/footeim/stdout/vcf_plink.%j.out
#SBATCH -e /nfs/scratch/footeim/stdout/vcf_plink.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module load vcftools/0.1.16
module load GCC/11.2.0
module load plink/1.90

#variables
PROJECT=$1
SET_NEW=$PROJECT'_'$2
VCF_NEW=$SCRATCH/projects/$PROJECT/data/Chapter4/snp/$SET_NEW/$SET_NEW'_raw'


#convert vcf file to plink
vcftools --gzvcf $VCF_NEW.vcf.gz --plink --out $VCF_NEW.plink

#convert PLINK to BED
plink --file $VCF_NEW.plink --make-bed --noweb --out $VCF_NEW
