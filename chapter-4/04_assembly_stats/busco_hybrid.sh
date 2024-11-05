#!/bin/bash
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=bigmem
#SBATCH --time=0-12:00
#SBATCH --job-name=BUSCO
#SBATCH -o /nfs/scratch/footeim/stdout/BUSCO.%j.out
#SBATCH -e /nfs/scratch/footeim/stdout/BUSCO.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

PROJECT=$1

module load busco/3.0.2
module load augustus/3.3.1
module load ncbi/blast+/2.7.1
module load hmmer/3.2.1
module load R/CRAN/3.6

BUSCO_DIR=$SCRATCH/scripts/Chapter3/04_assembly_stats/busco
BUSCO=$BUSCO_DIR/run_BUSCO.py
lineage=$BUSCO_DIR/aves_odb9

assembly_name=$2
sample_name=$3
assembly=$SCRATCH/projects/$PROJECT/data/Chapter3/flye/$assembly_name/assembly.fasta
out=$SCRATCH/projects/$PROJECT/output/Chapter3/flye_assembly/$assembly_name/$sample_name.busco

mkdir -p $SCRATCH/projects/$PROJECT/output/Chapter3/flye_assembly/$assembly_name/$sample_name.busco

export BUSCO_CONFIG_FILE="/nfs/scratch/footeim/scripts/Chapter3/04_assembly_stats/busco/config.ini"
export PATH="/home/software/apps/augustus/3.3.1/bin:$PATH"
export PATH="/home/software/apps/augustus/3.3.1/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/nfs/scratch/footeim/scripts/Chapter3/04_assembly_stats/busco/augustus/config"

#manual: https://busco.ezlab.org/

#OTHER FOLDER FOR BUSCO IN APP??:
#/home/software/apps/busco/3.0.2/lib/python3.6/site-packages/busco/
#GenomeAnalysis.py (Genome Assembly assessment)
#GeneSetAnalysis.py (Gene set (proteins) assessment)
#TranscriptomeAnalysis.py (Transcriptome assessment)
#BuscoAnalysis.py (all of the above?)

#BUSCO was used to locate the presence or absence of the Actinopterygii- specific set of 4584 single-copy orthologs (OrthoDB v9). 
#To assess the completeness of the assembly we used BUSCO v3 (16) with the Actinopterygii ortholog dataset. 

#The Bench- marking Universal Single-Copy Orthologs v.3.0.2 (BUSCO, RRID: SCR 015008)[29] programwith default settings 
#(e-value 0.01)was used to screen the Renilla genome assemblies for 978 orthologs from the Metazoan data set as a method 
#to evaluate the com- pleteness of each assembly. 
#BUSCO used BLAST v.2.2.31 [23]and HMMER v.3.1.b2 (HMMER, RRID:SCR 005305)[30] in its pipeline.

time python $BUSCO -i $assembly -o $SAMPLE_NAME'_busco' -l $lineage -m genome --cpu 10

mv 'run_'$SAMPLE_NAME'_busco' $out/ 
