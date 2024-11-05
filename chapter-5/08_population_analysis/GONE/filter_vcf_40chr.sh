#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --job-name=filter_vcf
#SBATCH -o /nesi/nobackup/vuw03922/stdout/filter_vcf.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/filter_vcf.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##the purpose of this script is to exclude contigs with only a few snps as GONE can't use these and they will throw an error

#load modules
module load BCFtools/1.19-GCC-11.3.0
module load R/4.3.1-gimkl-2022a
module load PLINK/1.09b6.16

#set variables
PROJECT=$1
SET=$PROJECT'_'$2
FILTER=$3

#define paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/${SET}_${FILTER}
#contigs=./contigs_to_include.txt

#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
vcf2R=$R/vcf2Rinput.R
gds2plink=$R/gds2plink.R

###resources
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_pop_info.tsv'


#filter to retain only first 40 chr
#bcftools view -Oz -r antipodean_1,antipodean_2,antipodean_3,antipodean_4,antipodean_5,antipodean_6,antipodean_7,antipodean_8,antipodean_9,antipodean_10,antipodean_11,antipodean_12,antipodean_13,antipodean_14,antipodean_15,antipodean_16,antipodean_17,antipodean_18,antipodean_19,antipodean_20,antipodean_21,antipodean_22,antipodean_23,antipodean_24,antipodean_25,antipodean_26,antipodean_27,antipodean_28,antipodean_29,antipodean_30,antipodean_31,antipodean_32,antipodean_33,antipodean_34,antipodean_35,antipodean_36,antipodean_37,antipodean_38,antipodean_39,antipodean_40 \
#	-o ${read_dir}_40chrom.vcf.gz	\
#	${read_dir}.vcf.gz
	
#bcftools view -Oz -r gibsons_1,gibsons_2,gibsons_3,gibsons_4,gibsons_5,gibsons_6,gibsons_7,gibsons_8,gibsons_9,gibsons_10,gibsons_11,gibsons_12,gibsons_13,gibsons_14,gibsons_15,gibsons_16,gibsons_17,gibsons_18,gibsons_19,gibsons_20,gibsons_21,gibsons_22,gibsons_23,gibsons_24,gibsons_25,gibsons_26,gibsons_27,gibsons_28,gibsons_29,gibsons_30,gibsons_31,gibsons_32,gibsons_33,gibsons_34,gibsons_35,gibsons_36,gibsons_37,gibsons_38,gibsons_39,gibsons_40 \
#	-o ${read_dir}_40chrom.vcf.gz	\
#	${read_dir}.vcf.gz	



#convert to gds
Rscript $vcf2R --gzvcf ${read_dir}_40chrom.vcf.gz \
               --snprelate_out ${read_dir}_40chrom

#convert to plink
Rscript $gds2plink --gds_file ${read_dir}_40chrom.gds \
                   --out ${read_dir}_40chrom \
                   --pop_file $pop_file

#make .ped and .map as well
plink --bfile ${read_dir}_40chrom --recode --out ${read_dir}_40chrom --chr-set 40
