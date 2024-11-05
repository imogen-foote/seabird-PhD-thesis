#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --time=0-12:00:00
#SBATCH --job-name=BUSCO
#SBATCH -o /nesi/nobackup/vuw03922/stdout/BUSCO.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/BUSCO.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load BUSCO/5.4.7-gimkl-2022a
module load AUGUSTUS/3.5.0-gimkl-2022a
module load BLAST/2.13.0-GCC-11.3.0
module load HMMER/3.3.2-GCC-11.3.0
module load R/4.3.1-gimkl-2022a

#set variables
PROJECT=$1
assembly_name=$2
assembler=$3

#define directories
#BUSCO_DIR=/nesi/nobackup/vuw03922/scripts/Chapter3/04_assembly_stats/busco
#BUSCO=$BUSCO_DIR/run_BUSCO.py
#lineage=$BUSCO_DIR/aves_odb9
assembly=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/$assembler/$assembly_name/*ssembly.fasta
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter3/${assembler}_assembly/$assembly_name/$assembly_name.busco

mkdir -p $out_dir

#BUSCO will try to use the dependencies available in your environment. 
#This can be overridden by passing a config file specifying the paths to the dependencies. 
#In the config/ subfolder of the cloned repository, a config.ini template is provided. 
#In this file, you may declare the paths to all third party components matching what is on your machine. 
#To activate this config file, set the environment variable BUSCO_CONFIG_FILE with the path to the file, as follows:
#export BUSCO_CONFIG_FILE="/nesi/nobackup/vuw03922/scripts/Chapter3/04_assembly_stats/busco/config.ini"

#Augustus requires environment variables to be declared as follows:
export PATH="/opt/nesi/CS400_centos7_bdw/modules/all/AUGUSTUS/3.5.0-gimkl-2022a/bin:$PATH"
export PATH="/opt/nesi/CS400_centos7_bdw/modules/all/AUGUSTUS/3.5.0-gimkl-2022a/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="/nesi/nobackup/vuw03922/scripts/Chapter3/04_assembly_stats/busco/augustus/config"

busco -i $assembly -o $PROJECT'_busco' -l aves_odb10 -m genome --cpu 10 --restart

#manual: https://busco.ezlab.org/

#time python $BUSCO -i $assembly -o $out_dir/$PROJECT'_busco' -l $lineage -m genome --cpu 10

mv 'run_'$PROJECT'_busco' $out/ 


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
