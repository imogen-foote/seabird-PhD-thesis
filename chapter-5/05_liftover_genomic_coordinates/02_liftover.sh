#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=0-0:05:00
#SBATCH --job-name=liftover
#SBATCH -o /nesi/nobackup/vuw03922/stdout/liftover.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/liftover.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load BCFtools/1.19-GCC-11.3.0

#define path to picard.jar
picard=/scale_wlg_persistent/filesets/opt_nesi/CS400_centos7_bdw/picard/2.26.10-Java-11.0.4/picard.jar

vcf=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter4/snp/GibsonsAlbatross_full
chain=/nesi/nobackup/vuw03922/projects/AllAlbatross/data/Chapter4/liftover/liftover_antipodean_gibsons/chainnet/liftover.chain
reference=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked
out_dir=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter4/liftover_snp
mkdir -p $out_dir

#check for presence of .dict file for genome, create if doesn't exitst
if [ -e $reference/blue-45G_softmasked_renamed.fasta.dict ]; then
	echo "blue-45G_softmasked_renamed.fasta.dict found."
else
	java -jar $picard CreateSequenceDictionary	\
		R=$reference/blue-45G_softmasked_renamed.fasta	\
		O=$reference/blue-45G_softmasked_renamed.fasta.dict
fi

#make tmp copy of vcf file and unzip 
mkdir $out_dir/tmp
cp $vcf/GibsonsAlbatross_full_qc.vcf.gz $out_dir/tmp/tmp_GibsonsAlbatross_full_qc.vcf.gz
gzip -d $out_dir/tmp/tmp_GibsonsAlbatross_full_qc.vcf.gz

#run liftover
java -jar $picard LiftoverVcf \
	I=$out_dir/tmp/tmp_GibsonsAlbatross_full_qc.vcf		\
	O=$out_dir/GibsonsAlbatross_liftover_full_qc.vcf	\
	CHAIN=$chain		\
	REJECT=$out_dir/GibsonsAlbatross_liftover_qc_reject.vcf \
	R=$reference/blue-45G_softmasked_renamed.fasta

#compress vcf file and index
bgzip -k -i $out_dir/GibsonsAlbatross_liftover_full_qc.vcf

#index vcf file
#bcftools index $out_dir/GibsonsAlbatross_liftover_full_qc.vcf.gz

#remove tmp directory
rm -r $out_dir/tmp


#Cluster: mahuika
#Job ID: 44062197
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   76.3%  00:03:49 of 00:05:00 time limit
#CPU Efficiency: 120.5%  00:04:36 of 00:03:49 core-walltime
#Mem Efficiency:  71.5%  14.30 GB of 20.00 GB

