#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:20:00
#SBATCH --job-name=map_chicken
#SBATCH -o /nesi/nobackup/vuw03922/stdout/map_chicken.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/map_chicken.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.30.0-GCC-11.3.0
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-01

###define paths
REFERENCE=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked/blue-45G_softmasked_renamed.fasta
regions=/nesi/nobackup/vuw03922/projects/AllAlbatross/data/Chapter4/outliers
out=/nesi/nobackup/vuw03922/projects/AllAlbatross/output/Chapter4/outlier_analyses/birdDB
mkdir -p $out
fasta=$out/fasta
mkdir -p $fasta
ref=/nesi/nobackup/vuw03922/resources/reference_genomes/bird_CDS/all_cds.fna



###EXTRACT SEQUENCE 1,000,000 BP either side of SNP

#Define a function to extract the region around each SNP
#extract_region() {
#    local chromosome="$1"
#    local position="$2"
#    local start="$((position - 1000000))"
#    local end="$((position + 1000000))"
#    local output_file="${regions}/fasta/${chromosome}_${position}.fasta"
#    samtools faidx $REFERENCE "${chromosome}:${start}-${end}" > "$output_file"
#}

#Extract regions
#while IFS=$'\t' read -r chromosome position; do
#    # Extract the region around the SNP
#    extract_region "$chromosome" "$position"
#done < $regions/AllAlbatross_full_qc_q0.05_K1_81_overlapping_outliers.tsv

###use regions identified in previous script
#extract fasta
extra=1000000
bedtools getfasta -fi $REFERENCE -bed $regions/outlier_regions_expanded_$extra.bed -fo $fasta/all_outlier_regions_$extra.fasta


#concatenate all fasta to one file
#cat $regions/fasta/antipodean*.fasta > $regions/fasta/all_outliers.fasta


###BLAST
export BLASTDB_LMDB_MAP_SIZE=100000000

#create blast db with bird protein data
makeblastdb -in $ref -dbtype nucl -out bird_cds_db
echo "blastdb created successfully"

#run blast
blastn -query $fasta/all_outlier_regions_$extra.fasta -db bird_cds_db -max_target_seqs 5 -out $out/bird_cds_BLAST_results_$extra -outfmt 6
echo "blast results written to $out/bird_cds_BLAST_results_$extra"

#remove blast db
rm bird_cds_db.*
echo "blastdb removed"

##extract unique values to a .csv
# Create an associative array to store unique values of the second column
declare -A unique_values

# Read the file line by line
while IFS=$'\t' read -r line; do
    # Extract the value of the second column
    second_column_value=$(echo "$line" | cut -f2)
    
    # Check if the value already exists in the associative array
    if [[ ! -v unique_values["$second_column_value"] ]]; then
        # If not, store the whole line in the associative array
        unique_values["$second_column_value"]="$line"
    fi
done < $out/bird_cds_BLAST_results_$extra

# Output the unique lines
for value in "${!unique_values[@]}"; do
    echo "${unique_values["$value"]}"
done > bird_cds_unique_$extra.csv
