#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=0:05:00
#SBATCH --job-name=minimap
#SBATCH -o /nesi/nobackup/vuw03922/stdout/minimap.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/minimap.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##Load modules
module purge
module load minimap2/2.24-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

genome=$1

###Set variables
antipodean=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked
gibsons=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/d36_q8_h50_masked
out_dir=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/output/Chapter3/synteny/minimap/gibsons_antipodean_q8_h50_softmasked_$genome
mkdir -p $out_dir
cd $out_dir


minimap2 -x asm5 -t 8 $gibsons/d36_softmasked_$genome.fasta $antipodean/blue-45G_softmasked_$genome.fasta > alignment.paf

#### Convert space into tab delimited and add header for table ####
	sed -e 's/ /\t/g' alignment.paf > alignment.paf.table
	sed -i '1 i\query\tlength\t0_based_start\tq_end\tstrand\tref\tlength\tstart\tr_end\tidentity' alignment.paf.table


#generate genome index files
samtools faidx $antipodean/blue-45G_softmasked_$genome.fasta
cp $antipodean/blue-45G_softmasked_$genome.fasta.fai $out_dir/antipodean.fasta.fai

samtools faidx $gibsons/d36_softmasked_$genome.fasta
cp $gibsons/d36_softmasked_$genome.fasta.fai $out_dir/gibsons.fasta.fai

#add a header
cat <(echo -e "chromosome\tsize\tnot_needed1\tnot_needed2\tnot_needed3") $out_dir/antipodean.fasta.fai > $out_dir/antipodean_circa.tsv
cat <(echo -e "chromosome\tsize\tnot_needed1\tnot_needed2\tnot_needed3") $out_dir/gibsons.fasta.fai > $out_dir/gibsons_circa.tsv




#Cluster: mahuika
#Job ID: 43432566
#State: COMPLETED
#Cores: 4
#Tasks: 1
#Nodes: 1
#Job Wall-time:    1.6%  00:01:54 of 02:00:00 time limit
#CPU Efficiency: 120.6%  00:09:10 of 00:07:36 core-walltime
#Mem Efficiency:  40.3%  12.09 GB of 30.00 GB

