#!/bin/bash
#SWATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:5:00
#SBATCH --job-name=get_overlap_outliers
#SBATCH -o /nesi/nobackup/vuw03922/stdout/get_overlap_outliers.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/get_overlap_outliers.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2
K=$3 #K values used for pcadapt

###load packages
module load HTSlib/1.19-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load BCFtools/1.19-GCC-11.3.0
module load R/4.3.1-gimkl-2022a
module load PLINK/1.09b6.16

#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
vcf2R=$R/vcf2Rinput.R
gds2plink=$R/gds2plink.R

###resrouces
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_site_info.tsv'

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
TMP=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/tmp
mkdir -p $TMP
outlier=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/outlier_analyses
PCADAPT=$outlier/pcadapt
OUTFLANK=$outlier/outflank

###extract all SNPS identified as outliers by pcadapt and outflank and create new txt file containing the names of all unique values
#extract SNP names identified by pcadapt
awk 'NR>1 {print $2}' $PCADAPT/*K${K}_q0.05_all_outliers.tsv > ${VCF}_pcadapt_outliers.txt

#extract SNP names identified by outflank
awk 'NR>1 {print $1}' $OUTFLANK/*_all_outliers_*.tsv > ${VCF}_outflank_outliers.txt

#combine
cat ${VCF}_pcadapt_outliers.txt ${VCF}_outflank_outliers.txt | sort -u > ${VCF}_all_outliers.txt


################################# filter neutral loci ###################################
echo "Removing outliers, filtering by hwe, maf and thinning..."
if [ -e ${VCF}_all_outliers.txt ]; then
 echo "${VCF}_all_outliers.txt found."
 #remove outliers, filter by hwe, maf and thin
 vcftools 	--gzvcf $VCF'_qc.vcf.gz'				\
			--out $TMP/$SET'_tmp6'				\
			--exclude-positions ${VCF}_all_outliers.txt 	\
			--hwe 0.05 --maf 0.05 --thin 10000		\
			--recode-INFO-all --recode
else
 echo "{VCF}_all_outliers.txt not found, proceeding with filtering and thinning."
 #filter by hwe, maf and thin
 vcftools 	--gzvcf $VCF'_qc.vcf.gz'			\
			--out $TMP/$SET'_tmp6'			\
			--hwe 0.05 --maf 0.05 --thin 10000	\
			--recode-INFO-all --recode
fi 
bgzip -ci $TMP/$SET'_tmp6.recode.vcf' > $TMP/$SET'_tmp6.vcf.gz'

#create gds from vcf so that can convert to plink format for pruning step
echo "Creating tmp gds..."
Rscript $vcf2R --gzvcf $TMP/$SET'_tmp6.vcf.gz'	\
		--snprelate_out $TMP/$SET'_tmp6'

#creating plink files from gds
echo "Creating tmp plink files..."
Rscript $gds2plink --gds_file $TMP/$SET'_tmp6.gds'	\
		--out $TMP/$SET'_tmp6'			\
		--pop_file $pop_file

#filter linked sites
echo "Filtering linked sites..."
plink --bfile $TMP/$SET'_tmp6'		\
	--indep-pairwise 50 5 0.2	\
	--chr-set 65			\
	--out $TMP/$SET'_tmp6'
sed -i 's/:/\t/g' $TMP/$SET'_tmp6.prune.in' #replace colons with tabs to make file compatible with vcftools
#perform pairwise LD pruning by removing one variant from a pair of SNPs if correlation coeff is greater than 0.2 within a 50-SNP window, advancing 5 SNPs at a time saves a file containing SNP 
#positions to keep - independent SNPs

#remove snps in LD as identified by previous stp to extract independent SNPs
echo "Extracting independent SNPs..."
vcftools --gzvcf $TMP/$SET'_tmp6.vcf.gz'	\
	--out $VCF'_neutral'			\
	--positions $TMP/$SET'_tmp6.prune.in'	\
	--recode-INFO-all --recode
#uses output file from previous step of SNPs to keep 
#convert to vcf.gz
mv $VCF'_neutral.recode.vcf' $VCF'_neutral.vcf'
bgzip $VCF'_neutral.vcf'

#convert to gds
echo "Creating gds..."
Rscript $vcf2R --gzvcf $VCF'_neutral.vcf.gz' \
		--snprelate_out $VCF'_neutral'
#convert to plink
echo "Creating plink files..."
Rscript $gds2plink --gds_file $VCF'_neutral.gds'	\
		--out $VCF'_neutral'			\
		--pop_file $pop_file
#clean up
#rm -r $TMP/
