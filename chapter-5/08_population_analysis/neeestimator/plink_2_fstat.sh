#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:05:00
#SBATCH --job-name=plink2fstat
#SBATCH -o /nesi/nobackup/vuw03922/stdout/plink2fstat.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/plink2fstat.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module purge
module load PLINK/1.09b6.16

#define paths
input=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/${SET}/${SET}_neutral
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/fstat
mkdir -p $out_dir

# Convert PLINK format (.bed, .bim, .fam) to text-based PED format
plink --bfile $input --recode --out $out_dir/${SET}_neutral --chr-set 65

# Convert PED format to Fstat format
awk '{print $1, $2, $5, $6}' $out_dir/${SET}_neutral.ped > $out_dir/${SET}_neutral_temp_fstat.txt

# Convert allele codes to 1 (presence) or 0 (absence)
sed -i 's/2/1/g'  $out_dir/${SET}_neutral_temp_fstat.txt
sed -i 's/1/1/g'  $out_dir/${SET}_neutral_temp_fstat.txt
sed -i 's/0/0/g'  $out_dir/${SET}_neutral_temp_fstat.txt
sed -i 's/ 2/ 0/g'  $out_dir/${SET}_neutral_temp_fstat.txt

# Add locus IDs at the beginning of each line
awk '{print NR, $0}'  $out_dir/${SET}_neutral_temp_fstat.txt >  $out_dir/${SET}_neutral_fstat

# Clean up temporary files
rm $out_dir/${SET}_neutral_temp_fstat.txt
