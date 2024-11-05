#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500M
#SBATCH --time=0:05:00
#SBATCH --job-name=extract_regions
#SBATCH -o /nesi/nobackup/vuw03922/stdout/extract_regions.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/extract_regions.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-01

###define paths
REFERENCE=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked/blue-45G_softmasked_renamed.fasta
outliers=/nesi/nobackup/vuw03922/projects/AllAlbatross/data/Chapter4/outliers
mkdir -p $outliers/fasta
blast=/nesi/nobackup/vuw03922/projects/AllAlbatross/output/Chapter4/outlier_analyses/BLAST
mkdir -p $blast


#Define a function to extract the region around each SNP
extract_region() {
    local chromosome="$1"
    local position="$2"
    local start="$((position - 1000000))"
    local end="$((position + 1000000))"
    local output_file="${outliers}/fasta/${chromosome}_${position}.fasta"
    samtools faidx $REFERENCE "${chromosome}:${start}-${end}" > "$output_file"
}

#Extract regions
while IFS=$'\t' read -r chromosome position; do
    # Extract the region around the SNP
    extract_region "$chromosome" "$position"
done < $outliers/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers.tsv

#concatenate all fasta to one file
cat $outliers/fasta/antipodean*.fasta > $outliers/fasta/all_outliers.fasta

#blast settings
FORMAT="6 qseqid stitle sseqid sallseqid sgi sallgi sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sblastnames sscinames"

#run blast
blastx -query $outliers/fasta/all_outliers.fasta \
       -db nr \
       -out $blast/outliers_blastx.txt \
       -evalue 1e-3 \
       -max_hsps 1 \
       -max_target_seqs 5 \
       -num_threads 64 \
       -outfmt "$FORMAT"
