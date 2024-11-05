#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0-0:15:00
#SBATCH --job-name=haplotype_caller
#SBATCH -o /nesi/nobackup/vuw03922/stdout/haplotype_caller.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/haplotype_caller.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load GATK/4.5.0.0-gimkl-2022a


#set environment
tax=$1
ref=$2

#define paths
read_dir=/nesi/nobackup/vuw03922/projects/BycatchID/data/paleomix_bam
reference=/nesi/nobackup/vuw03922/projects/BycatchID/resources/mitogenomes/${ref}_mitogenome
out_dir=/nesi/nobackup/vuw03922/projects/BycatchID/data/gvcf/${tax}
mkdir -p $out_dir
sample_list=/nesi/nobackup/vuw03922/projects/BycatchID/resources/sample_info/${tax}.txt

#check that reference mitogenome has fasta dict file required for gatk
if [ -f $reference.dict ];
then
	echo "Reference dictionary file found for $reference"
else
	#create dictionary file
	echo "Creating dictionary file for $reference..."
	gatk CreateSequenceDictionary -R $reference.fasta
	echo "Reference dictionary file for $reference created."
fi


#run gatk haplotype caller
while IFS= read -r sample
do

if [ -f "$out_dir/${sample}.mtgenome.g.vcf.gz" ];
then
	echo "$sample g.vcf.gz already present" 

else 

gatk HaplotypeCaller					\
	-R $reference.fasta				\
	-I $read_dir/${sample}_1.0.${ref}_mito.bam	\
	-ERC GVCF					\
	-ploidy 1					\
	-O $out_dir/${sample}.mtgenome.g.vcf.gz
	
fi 

done < $sample_list

#for 43 samples
#Cluster: mahuika
#Job ID: 44441437
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    7.7%  00:09:15 of 02:00:00 time limit
#CPU Efficiency: 173.7%  00:16:04 of 00:09:15 core-walltime
#Mem Efficiency:   7.8%  399.27 MB of 5.00 GB
