#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=0-0:30:00
#SBATCH --job-name=QC_merge
#SBATCH -o /nesi/nobackup/vuw03922/stdout/QC_merge.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/QC_merge.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module load HTSlib/1.19-GCC-11.3.0
module load BCFtools/1.19-GCC-11.3.0
module load R/4.3.2-foss-2023a

#custom scripts
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
vcf2R=$R/vcf2Rinput.R
gds2plink=$R/gds2plink.R

#tab delimited population file
POP=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_site_info.tsv'

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET'_qc'
TMP=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/tmp

#merge vcf files in tmp dir
bcftools concat -Oz -o $VCF.vcf.gz $( ls -v $TMP/*'_tmp5_qc.vcf.gz' ) --threads 10
bgzip --reindex $VCF.vcf.gz
tabix -p vcf $VCF.vcf.gz

#create gds
Rscript $vcf2R --gzvcf $VCF.vcf.gz	\
		--snprelate_out $VCF

#create plink bed file
Rscript $gds2plink --gds_file $VCF.gds	\
			--out $VCF		\
			--pop_file $POP

#clean up - remove the # if you're sure you dont need the temporary output
#rm -r $TMP
