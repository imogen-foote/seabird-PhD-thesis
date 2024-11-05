#!/bin/bash
#SWATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:5:00
#SBATCH --job-name=neutral_selection
#SBATCH -o /nesi/nobackup/vuw03922/stdout/neutral_selection.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/neutral_selection.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2
K=$3

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
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$SET'_site_info.tsv'

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
TMP=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/tmp
PCADAPT=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/outlier_analyses/pcadapt

### create directories
mkdir $TMP
################################# filter neutral loci  ###################################
#parameters you used to for outlier detection using pcadapt to remove the correct outliers
QVAL=0.05
#K=$3 above
#remove outlier loci
if [ -e $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' ]
then
 tail -n +2 $PCADAPT/$SET*'_K'$K'_q'$QVAL'_all_outliers.tsv' |  awk '{print $2}' | sed -e 's/:/\t/p' > $PCADAPT/$SET'_all_outlier_LOC.tsv'
 vcftools --gzvcf $VCF'_qc.vcf.gz'                                \
          --out $TMP/$SET'_tmp6'                                  \
          --exclude-positions $PCADAPT/$SET'_all_outlier_LOC.tsv' \
          --hwe 0.05 --maf 0.05 --thin 10000                      \
          --recode-INFO-all --recode
else
 vcftools --gzvcf $VCF'_qc.vcf.gz'          \
          --out $TMP/$SET'_tmp6'            \
          --hwe 0.05 --maf 0.05 --thin 10000 \
          --recode-INFO-all --recode		
fi
bgzip -ci $TMP/$SET'_tmp6.recode.vcf' > $TMP/$SET'_tmp6.vcf.gz'


#get vcf.gz
Rscript $vcf2R --gzvcf $TMP/$SET'_tmp6.vcf.gz' 	\
               --snprelate_out $TMP/$SET'_tmp6'
#get gds
Rscript $gds2plink --gds_file $TMP/$SET'_tmp6.gds' \
                   --out $TMP/$SET'_tmp6'          \
                   --pop_file $pop_file

#filter linked sites
plink --bfile $TMP/$SET'_tmp6'		\
	--indep-pairwise 50 5 0.2	\
	--chr-set 65			\
	--out $TMP/$SET'_tmp6'
sed -i 's/:/\t/g' $TMP/$SET'_tmp6.prune.in' 	#replace colons with tabs to make file compatible with vcftools
#perform pairwise LD pruning by removing one variant from a pair of SNPs if correlation coeff is greater than 0.2 within a 50-SNP window, advancing 5 SNPs at a time
#saves a file containing SNP positions to keep - independent SNPs?

# extract independent SNPs
vcftools --gzvcf $TMP/$SET'_tmp6.vcf.gz'		\
		--out $VCF'_neutral'					\
		--positions $TMP/$SET'_tmp6.prune.in'	\
		--recode-INFO-all --recode
#uses output file from previous step of SNPs to keep

#convert to vcf.gz
mv $VCF'_neutral.recode.vcf' $VCF'_neutral.vcf'
bgzip $VCF'_neutral.vcf'
tabix -p vcf $VCF'_neutral.vcf.gz'

#convert to gds
Rscript $vcf2R --gzvcf $VCF'_neutral.vcf.gz'	\
		--snprelate_out $VCF'_neutral'

#convert to plink					
Rscript $gds2plink --gds_file $VCF'_neutral.gds'	\
		--out $VCF'_neutral'						\
		--pop_file $pop_file
#clean up
rm -r $TMP/


#Cluster: mahuika
#Job ID: 44099217
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   26.3%  00:01:19 of 00:05:00 time limit
#CPU Efficiency:  92.4%  00:01:13 of 00:01:19 core-walltime
#Mem Efficiency:  12.1%  124.09 MB of 1.00 GB

