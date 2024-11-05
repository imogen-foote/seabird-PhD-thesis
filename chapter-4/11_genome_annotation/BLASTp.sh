#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=64
#SBATCH --mem=1G
#SBATCH --gres=ssd
#SBATCH --partition=milan
#SBATCH --time=0-1:0:00
#SBATCH --job-name=blastp
#SBATCH -o /nesi/nobackup/vuw03922/stdout/blastp.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/blastp.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load BLAST/2.13.0-GCC-11.3.0
module load BLASTDB/2024-01

#define variables
PROJECT=$1

#define paths
FASTA=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/ragtag/gene_annotation/braker/subset_test10.aa.fasta
BLAST_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/ragtag/gene_annotation/BLAST
mkdir -p $BLAST_DIR

#blast settings
FORMAT="6 qseqid stitle sseqid sallseqid sgi sallgi sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore sblastnames sscinames"

#copy the database to local SSD for tmp access during job
cp $BLASTDB/{nr,taxdb}.* $TMPDIR/
export BLASTDB=$TMPDIR

# retrieve the sequence ID based on the SLURM array task ID
#SEQ_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FASTA | awk '/^>/ {print $1}')

# run BLASTX (nucleotide query vs protein database)
blastp -query $FASTA \
       -db nr \
       -out $BLAST_DIR/subset_test.blastp.txt \
       -evalue 1e-3 \
       -max_hsps 1 \
       -max_target_seqs 10 \
       -num_threads 20 \
       -outfmt "$FORMAT"

#retrieve the sequence ID based on the SLURM array task ID
#SEQ_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $FASTA | awk '/^>/ {print substr($1, 2)}') 


#run BLASTX (nucleotide query vs protein database)
#blastx 	-query <(echo ">${SEQ_ID}" && sed -n -e "/>${SEQ_ID}/,/^>/ p" $FASTA | sed '1d') \
#	-db nr \
#	-out $BLAST_DIR/${SEQ_ID}.txt \
#	-evalue 1e-3 \
#	-max_hsps 1 \
#	-max_target_seqs 10 \
#	-num_threads 20 \
#	-outfmt $FORMAT  

echo "blastx for ${SEQ_ID} completed."

#evalue sets expectation value - number of hits you can expect to see by chance - threshold to 1e-3
#max_hsps sets maximum number of 'high-scoring segment pairs (HSPs)' - a local alignment with no gaps that achieves one of the highest alignment scores in a given search - to 1
#max_target_seqs sets maximum number of target sequences to report to 10

