#!/bin/bash -e
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --job-name=mito_rm
#SBATCH -o /nesi/nobackup/vuw03922/stdout/mito_rm.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/mitorm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load BLAST/2.13.0-GCC-11.3.0
#module load SeqKit/2.4.0

#set environment
PROJECT=$1
SAMPLE=$2
FILTER=$3

#set paths
#antipodean
#mitogenome=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/mitogenome/$SAMPLE'_'$FILTER/mitogenome_contig_only.fasta
#gibsons
mitogenome=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/mitogenome/$SAMPLE'_'$FILTER/animal_mt.complete.graph1.1.path_sequence.fasta
genome=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/3_purge_haplotigs/$SAMPLE'_'$FILTER'_purged'/assembly.fasta
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/4_nuclear_only/$SAMPLE'_'$FILTER'_nuclear'
mkdir -p $out_dir

export BLASTDB_LMDB_MAP_SIZE=100000000

makeblastdb -in $genome -dbtype nucl -out genome_reads

echo "blastdb created successfully"

blastn -query $mitogenome -db genome_reads -max_target_seqs 5 -out $out_dir/BLAST_results -outfmt 6

echo "blast results written to $out_dir/BLAST_results"

rm genome_reads.*

echo "blastdb removed"

#find the read with 99-100% match
#then remove on command line
#awk -v RS='>' '!/seq00000092/ {print ">"$0}' $genome > $out_dir/${SAMPLE}_nuclear.fasta 
#creates blank line at start for some reason - check with
head -c 100 $out_dir/${SAMPLE}_nuclear.fasta
#if present remove with
#sed -i '1d' $out_dir/${SAMPLE}_nuclear.fasta
