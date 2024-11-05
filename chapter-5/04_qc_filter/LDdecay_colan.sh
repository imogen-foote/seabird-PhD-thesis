#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=0:30:00
#SBATCH --job-name=LDdecay
#SBATCH -o /nesi/nobackup/vuw03922/stdout/LDdecay.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/LDdecay.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load HTSlib/1.19-GCC-11.3.0
module load PLINK/1.09b6.16
module load R/4.3.2-foss-2023a

###run input
PROJECT=$1
SET=$PROJECT'_'$2

#Rscript paths
LDdecay=./LDdecay_shell.R

###set paths
VCF=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/$SET 
OUT=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/LDdecay/$SET 

### create directories
mkdir -p $OUT

# calc ld with plink
echo "starting plink"

for CHR in {1..65}; do
	echo "Processing chromosome ${CHR}"

	plink	--bfile $VCF'_qc'	--allow-extra-chr	\
	--double-id	--set-missing-var-ids @:#	\
	--maf 0.01	--geno 0.1			\
	--mind 0.5	--chr-set 65			\
	--thin 0.1	-r2 gz	--ld-window 100		\
	--ld-window-kb 1000	--ld-window-r2 0	\
	--chr ${CHR}					\
	--out $OUT/${SET}_chr${CHR}_qc_plink
	
	python ./ld_decay_calc.py -i $OUT/${SET}_chr${CHR}_qc_plink.ld.gz -o $OUT/${SET}_chr${CHR}_maf0.01

	##remove these files each time whether run finishes or not or you'll have issues
	rm $OUT/*nosex $OUT/*bim $OUT/*bed $OUT/*fam $OUT/*gz $OUT/*log
done

#bgzip -cd $OUT/$SET'_qc_plink.ld.gz' | sed 's/[[:blank:]]/\t/g' | sed 's/^\t//g' | sed 's/\t$//g' > $OUT/$SET'_qc_plink.ld'



#visualise with R
#echo "creating LD plots..."
#Rscript ./LDdecay_colan.R $OUT/$SET
#Rscript $LDdecay $OUT/$SET'_qc'
