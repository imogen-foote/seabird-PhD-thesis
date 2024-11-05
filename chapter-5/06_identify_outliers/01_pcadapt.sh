#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0:05:00
#SBATCH --job-name=pcadapt
#SBATCH -o /nesi/nobackup/vuw03922/stdout/pcadapt.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/pcadapt.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2
K=$3

###load packages
module purge
module load HTSlib/1.19-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
#module load BCFtools/1.19-GCC-11.3.0
module load R/4.3.1-gimkl-2022a
#module load PLINK/1.09b6.16


#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
pcadapt=$R/pcadapt_shell.R
vcf2R=$R/vcf2Rinput.R

###resources
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_site_info.tsv'

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
PCADAPT=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/outlier_analyses/pcadapt

### create directories
#mkdir $TMP
mkdir -p $PCADAPT

##################################### identify outlier loci  #############################
QVAL=0.05  #value to determine outliers (loci with p-values less than or equal to q are considered outliers
mkdir -p $PCADAPT/q${QVAL}
#K=$3 as above
Rscript $pcadapt --plink $VCF'_qc'			\
			--out $PCADAPT/q${QVAL}/$SET	\
			--K0 10				\
			--K $K				\
			--maf 0.05			\
			--q $QVAL			\
			--slw 10000			\
			--minNsnp 2			\
			--mode full

#filter vcf file for outliers
#if [ -e $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' ]
#then
# #extract SNP info from outlier list
# tail -n +2 $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' |  awk '{print $2}' | sed -e 's/:/\t/p' > $PCADAPT/$SET'_all_outlier_LOC.tsv'
# #outputting outlier vcf file containing all snps
# vcftools --gzvcf $VCF'_qc.vcf.gz'				\
#		--positions $PCADAPT/$SET'_all_outlier_LOC.tsv'	\
#		--out $VCF'_outliers'				\
#		--recode-INFO-all				\
#		--recode
# mv $VCF'_outliers.recode.vcf' $VCF'_outliers.vcf'
# bgzip $VCF'_outliers.vcf' 

# Rscript $vcf2R --gzvcf $VCF'_outliers.vcf.gz'	\
#		--snprelate_out $VCF'_outliers'
#fi			

