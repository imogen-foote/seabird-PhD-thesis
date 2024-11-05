#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:5:00
#SBATCH --job-name=genotype_gatk
#SBATCH -o /nesi/nobackup/vuw03922/stdout/genotype_gatk.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/genotype_gatk.%j.err
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
#name of mtgenome contig - contig_10 for antipodean, d36_mitogenome for gibsons
interval=$3

#define paths
reference=/nesi/nobackup/vuw03922/projects/BycatchID/resources/mitogenomes/${ref}_mitogenome.fasta
out_dir=/nesi/nobackup/vuw03922/projects/BycatchID/data/vcf
mkdir -p $out_dir
gvcf_dir=/nesi/nobackup/vuw03922/projects/BycatchID/data/gvcf/${tax}
gvcf_list=/nesi/nobackup/vuw03922/projects/BycatchID/resources/sample_info/${tax}_gvcf.txt
database=/nesi/nobackup/vuw03922/projects/BycatchID/data/gatkDB
mkdir -p $database

#create GenomicsDB workspace from the list of gvcf files
gatk GenomicsDBImport						\
	--genomicsdb-workspace-path $database/${tax}_mito_DB/	\
	--create-output-variant-index				\
	--batch-size 0						\
	--sample-name-map $gvcf_list				\
	-L $interval						\
	--reader-threads 2

#perform joint genotyping on samples in the GenomicsDB workspace
gatk GenotypeGVCFs				\
	-R $reference				\
	-V gendb://$database/${tax}_mito_DB	\
	-G StandardAnnotation			\
	-ploidy 1				\
	-new-qual 				\
	-O $out_dir/${tax}_mito.vcf.gz

#unzip vcf
#gunzip -c $out_dir/${tax}_mito.vcf.gz > $out_dir/${tax}_mito.vcf

#filter based on Fisher Strand bias, Strand Odds Ratio, Mapping Quality, and QualByDepth
bcftools filter -i 'FS<60.0 && SOR<4 && MQ>30.0 && QD > 2.0'	\
	-Oz -o $out_dir/${tax}_mito.tmp.vcf.gz			\
	$out_dir/${tax}_mito.vcf.gz

#filter based on min genotype quality and min depth
vcftools	--gzvcf $out_dir/${tax}_mito.tmp.vcf.gz	\
		--minGQ 15	--minDP 3		\
		--recode --recode-INFO-all 		\
		--stdout | bgzip > $out_dir/${tax}_mito.QC.vcf.gz


#remove tmp vcf
rm $out_dir/${tax}_mito.tmp.vcf.gz



#Cluster: mahuika
#Job ID: 44475283
#State: FAILED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    3.4%  00:00:31 of 00:15:00 time limit
#CPU Efficiency: 106.5%  00:00:33 of 00:00:31 core-walltime
#Mem Efficiency:   0.1%  1.09 MB of 2.00 GB

