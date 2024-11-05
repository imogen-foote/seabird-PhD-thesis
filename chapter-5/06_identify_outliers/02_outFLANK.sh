#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:0:00
#SBATCH --job-name=outflank
#SBATCH -o /nesi/nobackup/vuw03922/stdout/outflank.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/outflank.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module purge
module load HTSlib/1.19-GCC-11.3.0
module load R/4.3.1-gimkl-2022a


#Rscript paths
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
outflank=$R/OutFlank_bash.R
vcf2R=$R/vcf2Rinput.R

###resources - outFLANK requires more than 2 pops to run
pop_file=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_site_info.tsv'

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET
OUTFLANK=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/outlier_analyses/outflank

### create directories
#mkdir $TMP
mkdir -p $OUTFLANK

##################################### identify outlier loci  #############################
QVAL=0.05
mkdir -p $OUTFLANK/q${QVAL}
Rscript $outflank	--vcf_file	$VCF'_qc.vcf.gz'	\
			--pop_file	$pop_file		\
			--fst_file	$OUTFLANK/q${QVAL}/$SET	\
			--out		$OUTFLANK/q${QVAL}/$SET	\
			--plots_out	$OUTFLANK/q${QVAL}/$SET	\
			--qval		$QVAL			\
			--filter_method	"FST"			\
			--slw		10000			\
			--min_snp	2			\
			--by_chr	TRUE			\
			--compress	TRUE			\
			--fst_plot	TRUE			\
			--qval_plot	TRUE			\
			--PCA_plot	TRUE			\
			--Manhattan_plot TRUE		



#Cluster: mahuika
#Job ID: 44337309
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   72.5%  00:43:31 of 01:00:00 time limit
#CPU Efficiency:  99.9%  00:43:28 of 00:43:31 core-walltime
#Mem Efficiency:  70.5%  7.05 GB of 10.00 GB

