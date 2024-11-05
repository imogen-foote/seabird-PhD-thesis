#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=00:5:00
#SBATCH --job-name=identify_regions
#SBATCH -o /nesi/nobackup/vuw03922/stdout/identify_regions.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/identify_regions.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###load packages
module purge
module load BEDTools/2.30.0-GCC-11.3.0

#define paths
genome_path=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye
braker=$genome_path/6_gene_annotation/braker_extra_protein
genome=$genome_path/5_repeat_masking/blue-45G_q8_h50_masked/blue-45G_softmasked_renamed.fasta
#blast=$genome_path/6_gene_annotation/BLAST/blastx.txt
regions=/nesi/nobackup/vuw03922/projects/AllAlbatross/data/Chapter4/outliers
out=/nesi/nobackup/vuw03922/projects/AllAlbatross/output/Chapter4/outlier_analyses/genome_annotation
mkdir -p $out

###Define regions
###APPROACH 1
###this approach takes XX distance (specified using bedtools slop command) either side of EVERY SNP###
#create bed file from snp coordinates file
#awk '{print $1"\t"$2"\t"$2+1}' $regions/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers.tsv > $regions/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers.bed
#expand snp coordinates to include a window of 1,000,000 bp
#bedtools slop -b 1000000 -i $regions/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers.bed -g $genome.fai > $regions/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers_expanded.bed

###APPROACH 2
###this approach groups snps that are within XX distance (specified by distance_threshold) of each other and groups them as one region, then expands this region by YY distance (specified by bedtools slop) 
#set distance threshold
distance_threshold=200000

#define function to group snps within distance threshold and print coordinates of region
process_snps() {
    local input_file=$1
    local current_chromosome=""
    local start_pos=""
    local end_pos=""

    while IFS=$'\t' read -r chromosome pos; do
        if [[ -z "$current_chromosome" || "$chromosome" != "$current_chromosome" ]]; then
            if [[ -n "$current_chromosome" ]]; then
                echo -e "$current_chromosome\t$start_pos\t$end_pos"
            fi
            current_chromosome="$chromosome"
            start_pos="$pos"
            end_pos="$pos"
        else
            if (( pos - end_pos <= distance_threshold )); then
                end_pos="$pos"
            else
                echo -e "$current_chromosome\t$start_pos\t$end_pos"
                start_pos="$pos"
                end_pos="$pos"
            fi
        fi
    done < "$input_file"

    if [[ -n "$current_chromosome" ]]; then
        echo -e "$current_chromosome\t$start_pos\t$end_pos"
    fi
}

#run function
process_snps $regions/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers.tsv > $regions/outlier_regions.bed

#expand snp coordinates to include a window of distance specified by extra variable
extra=1000000
bedtools slop -b $extra -i $regions/outlier_regions.bed -g $genome.fai > $regions/outlier_regions_expanded_$extra.bed


###Extract regions
#extract lines from .gff3 file that contain the snps
bedtools intersect -a $braker/braker.gff3 -b $regions/outlier_regions_expanded_$extra.bed -wa -wb > $out/outlier_regions_$extra.gff3

#extract the lines from the blast file where ID matches lines from gff3 extracted lines - can't do for braker_extra_protein because blast was performed on annotation without extra protein info
#grep -o 'ID=[^;]*' $out/outlier_regions_$extra.gff3 | cut -d= -f2 | while read id; do grep -wF "$id" $blast; done > $out/outlier_regions_${extra}_blastx.txt


###extract sequences of genes identified
# Extract unique IDs from outliers.gff3
awk -F '[\t;]' '{for (i=1; i<=NF; i++) if ($i ~ /^ID=/) {sub(/^ID=/, "", $i); if ($i ~ /\.t1$/) print $i}}' $out/outlier_regions_$extra.gff3 | sort | uniq > $out/unique_regions_sequence_ids_$extra.txt 


#extract sequences from braker.codingseq
# Iterate over each sequence ID from previous file
while IFS= read -r sequence_id; do
    # Search for the sequence ID in the FASTA file
    awk -v id=">$sequence_id" '$0 ~ id {print; getline; while (!/^>/) {print; getline; if ($0 ~ /^>/) exit;}}' $braker/braker.codingseq
done < $out/unique_regions_sequence_ids_$extra.txt > $out/outlier_regions_codingseqs_$extra.fasta
