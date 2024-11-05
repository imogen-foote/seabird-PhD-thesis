#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH -a 1-65
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0-0:30:00
#SBATCH --job-name=QC_filter
#SBATCH -o /nesi/nobackup/vuw03922/stdout/QC_filter.%A_%a.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/QC_filter.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


###!!!###                                       ###!!!####
# determine the number of scaffolds you want to genotype #
#               change --array accordingly               #
###!!!###                                       ###!!!####

###run input
PROJECT=$1
SET=$PROJECT'_'$2
N=${SLURM_ARRAY_TASK_ID}

###load packages
module load HTSlib/1.19-GCC-11.3.0
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
module load BCFtools/1.19-GCC-11.3.0
module load R/4.3.2-foss-2023a

#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
AB_script=$R/allelelic_imbalance_5.0.R

###resources
REF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/5_repeat_masking/*'_q8_h50_masked'/*'_softmasked_renamed.fasta'
#AB_exclude=$SCRATCH/projects/$PROJECT/resources/sample_info/high_coverage_samples.list #only if inital tests for allelic imbalance suggest removal of certain individuals.

###Obtain the scaffold name for the ith scaffold from your reference. A fai file should have been created when paoleomix indexed your genome.
LG=$( head -n $N $REF.fai | tail -n 1 | cut -f 1 )

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/tmp
mkdir -p $DIR
TMP=$DIR/$LG'_'$SET

##1## subsample vcf to select only current scaffold/linkage group 
bcftools view -Oz -o $TMP'_tmp1.vcf.gz' $VCF'_raw.vcf.gz' -r $LG
#-Oz save output as compressed vcf
#-r selects the region to subsample (current chromosome/scaffold/linkage group)

##2## filter genotypes with DP < 3 - doesn't remove, sets to missing
vcftools --gzvcf $TMP'_tmp1.vcf.gz'	\
	--out 	$TMP'_tmp2'				\
	--minDP 3						\
	--remove-indels					\
	--recode-INFO-all	--recode
mv $TMP'_tmp2.recode.vcf' $TMP'_tmp2.vcf'
bgzip $TMP'_tmp2.vcf'
#filters out snps in a sample if the depth of coverage for that snp in that samples is less than 3
#remove-indels removes all sites that contain an indel
#recode-INFO-all ensures all INFO fields associated with variants that pass the filtering are retained in output, --recode tells VCFtools to create a new VCF file containing only the variatns that meet the filtering critera (in this case minDP 3)

##3## basic filter parameters
vcftools --gzvcf $TMP'_tmp2.vcf.gz'				\
		--out $TMP'_tmp3'	--max-missing 0.95	\
		--min-alleles 2		--max-alleles 2		\
		--min-meanDP 8		--max-meanDP 30		\
		--minQ 600			--maf 0.01			\
		--recode-INFO-all	--recode
mv $TMP'_tmp3.recode.vcf' $TMP'_tmp3.vcf'
bgzip $TMP'_tmp3.vcf' 
# max-missing excludes sites with more than 5% mssing data
# min/max-alleles excludes sites that are not biallelic
# min/max meanDP excludes sites with mean depth less than 8 and greater than 30 - basically want to remove anything with super high depth cos repetitive regions etc. Different to above as filters based on the mean depth across all samples for a given variant.
# include only sites with Q>600 and minor allele frequency > 0.01


##4## test for allelic imbalance
#Select output from VCF (genotypes)
vcftools --gzvcf $TMP'_tmp3.vcf.gz'	\
		--out $TMP'_tmp4'			\
		--extract-FORMAT-info GT
#extract-FORMAT-info GT extracts all genotype entries from VCF file

#Select output from VCF (allelic depth)
vcftools --gzvcf $TMP'_tmp3.vcf.gz' \
		--out 	$TMP'_tmp4'			\
		--extract-FORMAT-info AD
#extract-FORMAT-info AD extracts all allelic depth entries from VCF file

#run binomial test to filter sites with allelic imbalance - could require high mem when many SNPs are to be analysed
Rscript $AB_script --GT_file  $TMP'_tmp4.GT.FORMAT'		\
					--AD_file  $TMP'_tmp4.AD.FORMAT'	\
					--out_file $TMP'_tmp4_qc'			\
					--conf.level 0.99					\
					--plots TRUE
#                   --remove $AB_exclude	### make a list of sample names you want to exclude

##5## remove sites suffering of allelic imbalance
vcftools --gzvcf $TMP'_tmp3.vcf.gz'		\
		--out $TMP'_tmp5_qc'			\
		--recode-INFO-all	--recode	\
		--exclude-positions $TMP'_tmp4_qc.exclude_pval0.01.list'
mv $TMP'_tmp5_qc.recode.vcf' $TMP'_tmp5_qc.vcf'
bgzip -fi $TMP'_tmp5_qc.vcf'
tabix -fp vcf $TMP'_tmp5_qc.vcf.gz'



#Cluster: mahuika
#Job ID: 43722191
#Array Job ID: 43721687_1
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    0.3%  00:08:04 of 2-00:00:00 time limit
#CPU Efficiency:  97.7%  00:07:53 of 00:08:04 core-walltime
#Mem Efficiency:  15.7%  1.26 GB of 8.00 GB

