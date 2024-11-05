#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:10:00
#SBATCH --job-name=create_fasta
#SBATCH -o /nesi/nobackup/vuw03922/stdout/create_fasta.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/create_fasta.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load GATK/4.5.0.0-gimkl-2022a
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1 
module load BCFtools/1.13-GCC-9.2.0

#set environment
tax=$1
ref=$2

#define paths
read_dir=/nesi/nobackup/vuw03922/projects/BycatchID/data/vcf
reference=/nesi/nobackup/vuw03922/projects/BycatchID/resources/mitogenomes/${ref}_mitogenome.fasta
out_dir=/nesi/nobackup/vuw03922/projects/BycatchID/data/fasta/${tax}
mkdir -p $out_dir
sample_list=/nesi/nobackup/vuw03922/projects/BycatchID/resources/sample_info/${tax}.txt

while IFS= read -r sample
do

#select individual VCF
echo $sample > $out_dir/$sample.list
vcftools --gzvcf $read_dir/${tax}_mito.QC.vcf.gz	\
	--keep $out_dir/$sample.list			\
	--recode --recode-INFO-all			\
	--out $out_dir/$sample

#Change missing sites to N and set this as alternative allele
awk '{if ($10 ~ "\\.:" ) {print $1"\t"$2"\t"$3"\t"$4"\t""N""\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} else {print $0}}' $out_dir/$sample.recode.vcf | sed 's/GT:AD:DP:GQ:PL\t\./GT:AD:DP:GQ:PL\t1/g' > $out_dir/$sample.tmp.vcf 

#output alternative and missing alleles only
vcftools --vcf $out_dir/$sample.tmp.vcf	\
	--non-ref-ac 1			\
	--recode --recode-INFO-all	\
	--out $out_dir/$sample.alternative
#index
gatk IndexFeatureFile -I $out_dir/$sample.alternative.recode.vcf

#create fasta for sample
gatk FastaAlternateReferenceMaker	\
	-R $reference			\
	-O $out_dir/$sample.fasta	\
	-V $out_dir/$sample.alternative.recode.vcf

#update sample name in fasta
sed -i "s/>.*/>${sample}/" $out_dir/$sample.fasta

#clean up
rm $out_dir/$sample.list
rm $out_dir/$sample.dict
rm $out_dir/$sample.fasta.fai
rm $out_dir/$sample.recode.vcf
rm $out_dir/$sample.tmp.vcf
rm $out_dir/$sample.alternative.recode.vcf
rm $out_dir/$sample.alternative.recode.vcf.idx

done < $sample_list

