#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=0-2:00:00
#SBATCH --job-name=bcf_concat
#SBATCH -o /nesi/nobackup/vuw03922/stdout/bcf_concat.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/bcf_concat.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module load HTSlib/1.19-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load BCFtools/1.19-GCC-11.3.0

#variables
PROJECT=$1
SET_NEW=$PROJECT'_'$2

#path to new (merged) VCF file
VCF_NEW=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET_NEW/$SET_NEW'_raw'

#path where unmerged VCF files are located
TMP_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET_NEW/tmp

#merge vcf files in tmp dir
bcftools concat -Oz -o $VCF_NEW'.vcf.gz' $( ls -v $TMP_DIR/*'_raw_tmp2.vcf.gz' )
#$(ls ....) command substitution with the -v option lists and accesses files in numerical order so that scaffolds are ordered numerically in the resulting VCF file


#create indexes for raw vcf file
bgzip --reindex $VCF_NEW'.vcf.gz'
tabix -p vcf $VCF_NEW'.vcf.gz'

#remove tmp folder when you don't need it anymore
#rm -r $TMP_DIR


#Cluster: mahuika
#Job ID: 43612476
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    2.8%  00:08:28 of 05:00:00 time limit
#CPU Efficiency: 166.9%  00:14:08 of 00:08:28 core-walltime
#Mem Efficiency:   0.5%  20.97 MB of 4.00 GB

