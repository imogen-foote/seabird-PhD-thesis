#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=0:05:00
#SBATCH --job-name=circa
#SBATCH -o /nesi/nobackup/vuw03922/stdout/circa.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/circa.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0

#define variables and paths
antipodean=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked
gibsons=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/d36_q8_h50_masked
out_dir=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/output/Chapter3/synteny/mashmap/gibsons_antipodean_q8_h50_softmasked
#out_dir=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/output/Chapter3/synteny/mashmap/antipodean_gibsons_q8_h50_softmasked

#generate genome index files
#samtools faidx $antipodean/blue-45G_nuclear_sorted_renamed.fasta 
#samtools faidx $antipodean/blue-45G_softmasked_renamed.fasta
#cp $antipodean/blue-45G_nuclear_sorted_renamed.fasta.fai $out_dir/antipodean.fasta.fai
cp $antipodean/blue-45G_softmasked_renamed.fasta.fai $out_dir/antipodean.fasta.fai
#samtools faidx $gibsons/d36_nuclear_sorted_renamed.fasta
#samtools faidx $gibsons/d36_softmasked_renamed.fasta
#cp $gibsons/d36_nuclear_sorted_renamed.fasta.fai $out_dir/gibsons.fasta.fai
cp $gibsons/d36_softmasked_renamed.fasta.fai $out_dir/gibsons.fasta.fai

#add a header
cat <(echo -e "chromosome\tsize\tnot_needed1\tnot_needed2\tnot_needed3") $out_dir/antipodean.fasta.fai > $out_dir/antipodean_circa.tsv
cat <(echo -e "chromosome\tsize\tnot_needed1\tnot_needed2\tnot_needed3") $out_dir/gibsons.fasta.fai > $out_dir/gibsons_circa.tsv

##need to rename contig_X, scaffold_X to gibsons_X and antipodean_X so unique identifiers for each genome
#in .tsv files
#awk 'BEGIN{FS=OFS="\t"} {sub(/^contig_/, "gibsons_", $1)}1' $out_dir/gibsons_circa.tsv | \
#awk 'BEGIN{FS=OFS="\t"} {sub(/^scaffold_/, "gibsons_", $1)}1' > $out_dir/gibsons_circa.tsv

#awk 'BEGIN{FS=OFS="\t"} {sub(/^contig_/, "antipodean_", $1)}1' $out_dir/antipodean_circa.tsv | \
#awk 'BEGIN{FS=OFS="\t"} {sub(/^scaffold_/, "antipodean_", $1)}1' > $out_dir/antipodean_circa.tsv

#in mashmap table
#awk 'BEGIN{FS=OFS="\t"} {sub(/^contig_/, "gibsons_", $1)}1' $out_dir/mashmap.out.table | \
#awk 'BEGIN{FS=OFS="\t"} {sub(/^scaffold_/, "gibsons_", $1)}1' | \
#awk 'BEGIN{FS=OFS="\t"} {sub(/^contig_/, "antipodean_", $6)}1' | \
#awk 'BEGIN{FS=OFS="\t"} {sub(/^scaffold_/, "antipodean_", $1)}1' > $out_dir/mashmap.out.table
