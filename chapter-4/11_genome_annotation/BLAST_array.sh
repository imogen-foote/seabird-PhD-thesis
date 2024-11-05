#!/bin/bash -e
#SBATCH --account=vuw03922
#SBATCH --array=1-30
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH --gres=ssd
#SBATCH --partition=milan
#SBATCH --time=4-00:0:00
#SBATCH --job-name=blast_array
#SBATCH -o /nesi/nobackup/vuw03922/stdout/blast_array.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/blast_array.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##NOTE: memory requested must be sufficient to allow the blast db to be stored to tmp_dir (157GB for refseq_protein as at 2023)

#load modules
module purge
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-01

#define variables
PROJECT=$1

#define paths
READ_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/ragtag/gene_annotation/braker/subset_codingseq
BLAST_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/ragtag/gene_annotation/BLAST
mkdir -p $BLAST_DIR

#on command line split big fasta file into 30 smaller, each with 1007 reads (30205/30)
#awk "BEGIN {
#  n_seq=0;
#}
#/^>/ {
#  n_seq++;
#  if (n_seq % 1007 == 1) {
#   file=sprintf(\"subset_codingseq/codingseq_%d.fa\", (n_seq-1)/1007 + 1);
#   close(file);
#  }
#  print >> file;
#  next;
#}
#{
#  print >> file;
#}" < braker.codingseq

#blast settings
FORMAT="6 qseqid stitle sseqid sallseqid sgi sallgi sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sblastnames sscinames"

#copy the database to local SSD for tmp access during job
cp $BLASTDB/{nr,taxdb}.* $TMPDIR/
export BLASTDB=$TMPDIR

# run BLASTX (nucleotide query vs protein database)
blastx -query $READ_DIR/codingseq_${SLURM_ARRAY_TASK_ID}.fa \
       -db nr \
       -out $BLAST_DIR/blastx_results_${SLURM_ARRAY_TASK_ID}.txt \
       -evalue 1e-3 \
       -max_hsps 1 \
       -max_target_seqs 10 \
       -num_threads 20 \
       -outfmt "$FORMAT"
	   
#merge all blastx output
cat blastx_results_{1..30}.txt > all_blastx_results.txt

#evalue sets expectation value - number of hits you can expect to see by chance - threshold to 1e-3
#max_hsps sets maximum number of 'high-scoring segment pairs (HSPs)' - a local alignment with no gaps that achieves one of the highest alignment scores in a given search - to 1
#max_target_seqs sets maximum number of target sequences to report to 10
