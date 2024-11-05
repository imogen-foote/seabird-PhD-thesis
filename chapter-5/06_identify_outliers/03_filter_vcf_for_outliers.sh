#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --job-name=filtervcf
#SBATCH -o /nesi/nobackup/vuw03922/stdout/filtervcf.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/filtervcf.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###this script requires a file containing the overlapping outliers identified in pcadapt and outflank - this is done in R###

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module purge
module load HTSlib/1.19-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load R/4.3.1-gimkl-2022a

q=0.05
K=1

#define paths
vcf2R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts/vcf2Rinput.R
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_site_info.tsv'
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
outliers=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/outliers


#filter vcf file for outliers
if [ -e $outliers/${SET}_qc_q${q}_K${K}_21_overlapping_outliers.tsv ]
then
 #extract SNP info from outlier list
 #tail -n +2 $outliers/${SET}_qc_q${q}_K${K}_*_overlapping_outliers.tsv |  awk '{print $2}' | sed -e 's/:/\t/p' > $outliers/$SET'_all_outlier_LOC.tsv'
 #outputting outlier vcf file containing all snps
 vcftools --gzvcf $VCF'_qc.vcf.gz'                             				\
               --positions $outliers/${SET}_qc_q${q}_K${K}_21_overlapping_outliers.tsv	\
               --out $VCF'_outliers_21snps'                           				\
               --recode-INFO-all                               				\
               --recode
 mv $VCF'_outliers_21snps.recode.vcf' $VCF'_outliers_21snps.vcf'
 bgzip $VCF'_outliers_21snps.vcf'

 Rscript $vcf2R --gzvcf $VCF'_outliers_21snps.vcf.gz' 	\
               	--snprelate_out $VCF'_outliers_21snps'	\
		--genlight_out $VCF'_outliers_21snps'
fi


#Cluster: mahuika
#Job ID: 45212750
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   21.7%  00:01:05 of 00:05:00 time limit
#CPU Efficiency:  92.3%  00:01:00 of 00:01:05 core-walltime
#Mem Efficiency:   4.7%  241.99 MB of 5.00 GB

