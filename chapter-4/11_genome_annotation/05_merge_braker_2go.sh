#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:05:00
#SBATCH --job-name=merge_braker_b2go
#SBATCH -o /nesi/nobackup/vuw03922/stdout/merge_braker_b2go.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/merge_braker_b2go.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge

#define variables
PROJECT=$1

#define paths
annotation=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation
braker=$annotation/braker_extra_protein/braker.gff3
b2go=$annotation/blast2go/blast2go_table.txt

# Parse BLAST2GO output and store relevant information in a dictionary
declare -A blast_info
declare -A go_info

while IFS=$'\t' read -r tags seq_name description length hits e_value sim_mean go_count go_ids go_names enzyme_codes enzyme_names interpro_ids interpro_go_ids interpro_go_names; do
    if [ "$tags" == "TRUE" ]; then
        blast_info["$seq_name"]="$hits|$e_value"  # Store BLAST hits and e-value
        go_info["$seq_name"]="$go_ids|$go_names" # Store GO IDs and names
    fi
done < $b2go

# Modify GFF3 file to include BLAST and GO annotations
while IFS=$'\t' read -r contig feature type start end dot strand dot attributes; do
    if [[ "$type" == "mRNA" ]]; then
        seq_name=$(echo "$attributes" | grep -oP 'ID=\K[^;]+')
        blast_annotation="${blast_info[$seq_name]}"
        go_annotation="${go_info[$seq_name]}"
        # Append BLAST and GO annotations to the attributes
        attributes="$attributes;blast_annotation=$blast_annotation;go_annotation=$go_annotation"
    fi
    echo -e "$contig\t$feature\t$type\t$start\t$end\t$dot\t$strand\t$dot\t$attributes"
done < $braker > $annotation/blast2go/blast2go.gff3
