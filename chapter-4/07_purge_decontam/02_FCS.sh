#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=48
#SBATCH --mem=500G
#SBATCH --time=05:00:00
#SBATCH --job-name=fcs
#SBATCH -o /nesi/nobackup/vuw03922/stdout/fcs.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/fcs.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#48cpu 500Gb
#load modules
module purge 
module load Singularity/3.11.3
module load Python/3.11.3-gimkl-2022a

export FCS_DEFAULT_IMAGE=fcs-gx.sif

PROJECT=$1
SAMPLE=$2
FILTER=$3
ASSEMBLER=$4

#specify path to database
GXDB_LOC=/nesi/nobackup/vuw03922/resources/FCS/gxdb

#specify path to genome
genome=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/$ASSEMBLER/3_purge_haplotigs/$SAMPLE'_'$FILTER'_purged'/assembly.fasta 
#specify output
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/$ASSEMBLER/4_FCS/$SAMPLE'_'$FILTER'_FCS'
mkdir $out_dir

#test run with test dataset
#python3 ./fcs.py --image fcs-gx.sif screen genome --fasta ./fcsgx_test.fa.gz --out-dir $out_dir --gx-db "$GXDB_LOC/test-only" --tax-id 6973

#full run with tax-id for 'aves'
python3 ./fcs.py --env-file env.txt --image fcs-gx.sif screen genome --fasta $genome --out-dir $out_dir --gx-db "$GXDB_LOC" --tax-id 8782

cp /nesi/nobackup/vuw03922/stdout/fcs.$SLURM_JOB_ID.err $out_dir



#Cluster: mahuika
#Job ID: 42180872
#State: FAILED
#Cores: 24
#Tasks: 1
#Nodes: 1
#Job Wall-time:   24.6%  02:57:07 of 12:00:00 time limit
#CPU Efficiency:   5.5%  03:54:05 of 2-22:50:48 core-walltime
#Mem Efficiency:   0.9%  4.27 GB of 500.00 GB ##but seems to require all 500 to cache the db to memory

