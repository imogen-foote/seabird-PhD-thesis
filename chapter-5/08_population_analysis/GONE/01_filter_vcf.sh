#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:15:00
#SBATCH --job-name=filter_vcf
#SBATCH -o /nesi/nobackup/vuw03922/stdout/filter_vcf.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/filter_vcf.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

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


###filter to retain only antips
bcftools view -Oz -s ANT_DG46,ANT_DG47,ANT_DG48,ANT_DG49,ANT_DG50,ANT_MW16,ANT_MW17,ANT_MW18,ANT_MW19,ANT_MW21,ANT_MW22,ANT_MW24,ANT_NP28,ANT_NP30,ANT_NP31,ANT_NP32,ANT_NP33,ANT_NP34,ANT_NP35,ANT_NP36,ANT_NP39,ANT_NP40,ANT_NP41,ANT_NP42,ANT_NP43,ANT_NP44,ANT_OR01,ANT_OR02,ANT_OR03,ANT_OR04,ANT_OR05,ANT_OR06,ANT_OR07,ANT_OR08,ANT_OR10,ANT_OR11,ANT_OR13,ANT_OR14,ANT_OR15,ANT_SA51,ANT_SA52,ANT_SA53,ANT_SA54	\
	-o ${read_dir}_antipodean.vcf.gz	\
	${read_dir}.vcf.gz	

#convert to gds
Rscript $vcf2R --gzvcf ${read_dir}_antipodean.vcf.gz \
               --snprelate_out ${read_dir}_antipodean

#convert to plink
Rscript $gds2plink --gds_file ${read_dir}_antipodean.gds \
                   --out ${read_dir}_antipodean \
                   --pop_file $pop_file

#make .ped and .map as well
plink --bfile ${read_dir}_antipodean --recode --out ${read_dir}_antipodean --chr-set 65	

###filter to retain only gibsons
bcftools view -Oz -s GIB_AA01,GIB_AA04,GIB_AA05,GIB_AA07,GIB_AA08,GIB_AA09,GIB_AA13,GIB_AA15,GIB_AA16,GIB_AA17,GIB_AA19,GIB_AA22,GIB_AA24,GIB_AA25,GIB_AA26,GIB_AA27,GIB_AA29,GIB_AA30,GIB_AA31,GIB_DI48,GIB_DI49,GIB_DI50,GIB_DI51,GIB_DI52,GIB_DI53,GIB_DI54,GIB_DI55,GIB_MD32,GIB_MD33,GIB_MD34,GIB_MD35,GIB_MD36,GIB_MD37,GIB_MD38,GIB_MD39,GIB_MD40,GIB_MD41,GIB_MD42,GIB_MD43,GIB_MD44,GIB_MR45,GIB_MR46,GIB_MR47	\
	-o ${read_dir}_gibsons.vcf.gz	\
	${read_dir}.vcf.gz

#convert to gds
Rscript $vcf2R --gzvcf ${read_dir}_gibsons.vcf.gz \
               --snprelate_out ${read_dir}_gibsons

#convert to plink
Rscript $gds2plink --gds_file ${read_dir}_gibsons.gds \
                   --out ${read_dir}_gibsons \
                   --pop_file $pop_file

#make .ped and .map as well
plink --bfile ${read_dir}_gibsons --recode --out ${read_dir}_gibsons --chr-set 65	

