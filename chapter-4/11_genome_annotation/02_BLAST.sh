#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --gres=ssd
#SBATCH --partition=milan
#SBATCH --time=7-0:00:00
#SBATCH --job-name=blastx
#SBATCH -o /nesi/nobackup/vuw03922/stdout/blastx.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/blastx.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#64Gb, 64cpu 1.5day
##NOTE: memory requested must be sufficient to allow the blast db to be stored to tmp_dir (157GB for refseq_protein as at 2023) unless copying to ssd on milan node

#load modules
module purge
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-01

#define variables
PROJECT=$1

#define paths
FASTA=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation/braker/braker.codingseq
BLAST_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/6_gene_annotation/BLAST
mkdir -p $BLAST_DIR

#blast settings
FORMAT="6 qseqid stitle sseqid sallseqid sgi sallgi sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sblastnames sscinames"

#copy the database to local SSD for tmp access during job
cp $BLASTDB/{nr,taxdb}.* $TMPDIR/
export BLASTDB=$TMPDIR

# run BLASTX (nucleotide query vs protein database)
blastx -query $FASTA \
       -db nr \
       -out $BLAST_DIR/blastx.txt \
       -evalue 1e-3 \
       -max_hsps 1 \
       -max_target_seqs 10 \
       -num_threads 64 \
       -outfmt "$FORMAT"

#echo "blastx for ${SEQ_ID} completed."

#evalue sets expectation value - number of hits you can expect to see by chance - threshold to 1e-3
#max_hsps sets maximum number of 'high-scoring segment pairs (HSPs)' - a local alignment with no gaps that achieves one of the highest alignment scores in a given search - to 1
#max_target_seqs sets maximum number of target sequences to report to 10

