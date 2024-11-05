#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH -a 1-65
#SBATCH --mem=1G
#SBATCH --time=1-00:00:00
#SBATCH --job-name=mpileup
#SBATCH -o /nesi/nobackup/vuw03922/stdout/mpileup.%A_%a.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/mpileup.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###!!!###		       			  	 ###!!!####
# determine the number of scaffold you want to genotype            #
#               change --array accordingly		           #
###!!!###						  ###!!!####

#load modules
module load HTSlib/1.19-GCC-11.3.0
module load BCFtools/1.19-GCC-11.3.0

#variables
SCAFFOLD=${SLURM_ARRAY_TASK_ID}
PROJECT=$1
SAMPLE=$2
SET_NEW=$PROJECT'_'$3
TMP_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET_NEW/tmp
mkdir -p $TMP_DIR

#genome
REFERENCE=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/5_repeat_masking/$SAMPLE'_q8_h50_masked'/$SAMPLE'_softmasked_renamed.fasta'
FAI=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/5_repeat_masking/$SAMPLE'_q8_h50_masked'/$SAMPLE'_softmasked_renamed.fasta.fai'


BAMLIST=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_bam.list'
REGION=$( head -n $SCAFFOLD $FAI | tail -n 1 | cut -f 1 )

#genotype
bcftools mpileup -Ov									\
		-a 'FORMAT/AD,FORMAT/DP,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/AD'	\
		-f $REFERENCE								\
		-r $REGION								\
		-b $BAMLIST								|
bcftools call -Ov -mv > $TMP_DIR/$SCAFFOLD'_'$REGION'_'$SET_NEW'_raw_tmp1.vcf'
#do I need to remove -v flag to make it look at all sites rather than just variable sites

#update INFO fields
bcftools	+fill-tags 	$TMP_DIR/$SCAFFOLD'_'$REGION'_'$SET_NEW'_raw_tmp1.vcf'			\
		-Oz -o 		$TMP_DIR/$SCAFFOLD'_'$REGION'_'$SET_NEW'_raw_tmp2.vcf.gz'		\
		-- -t AC,AF,AN,MAF,NS,AC_Hom,AC_Het			
bgzip --reindex $TMP_DIR/$SCAFFOLD'_'$REGION'_'$SET_NEW'_raw_tmp2.vcf.gz'

rm $TMP_DIR/$SCAFFOLD'_'$REGION'_'$SET_NEW'_raw_tmp1.vcf'

