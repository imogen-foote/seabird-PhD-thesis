#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=32
#SBATCH --mem=30G
#SBATCH --time=2-0:0:00
#SBATCH --job-name=repeat_modeler
#SBATCH -o /nesi/nobackup/vuw03922/stdout/repeat_modeler.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/repeat_modeler.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load Apptainer/1.2.2

#export PATH to libraries created
export LIBDIR=/nesi/nobackup/vuw03922/resources/repeat_masker/Libraries

#parameters
PROJECT=$1
SAMPLE=$2
FILTER=$3

#PATHS
ASSEMBLY_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/4_nuclear_only/$SAMPLE'_'$FILTER'_nuclear' 
REPEAT_DIR=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/5_repeat_masking/$SAMPLE'_'$FILTER'_masked' 
mkdir -p $REPEAT_DIR

#1. Create a database for RepeatModeler - build database of repeats found in your assembly
#mkdir $REPEAT_DIR/libraries
#singularity run dfam-tetools-latest.sif BuildDatabase -name $REPEAT_DIR/libraries/repeat_modeler_db \
#                                                      -engine ncbi                                  \
#                                                      $ASSEMBLY_DIR/$SAMPLE'_nuclear.fasta'


#2. Run RepeatModeler - identify TEs
singularity run dfam-tetools-latest.sif RepeatModeler -threads 32                                 \
                                                      -LTRStruct                                  \
                                                      -database $REPEAT_DIR/libraries/repeat_modeler_db \
						-recoverDir ./RM_125449.ThuJan181143212024


#3. check if taxonomic group is available and extract repeats
####DO THIS BIT ON COMMAND LINE ####
##enter singularity
#module load Singularity
#singularity shell dfam-tetools-latest.sif
#famdb.py  ##for info
##to print taxonomy nodes that include 'mus', and the corresponding IDs.
##    The IDs and names are stored in the FamDB file, and are based
##    on the NCBI taxonomy database (https://www.ncbi.nlm.nih.gov/taxonomy).
#famdb.py -i /nesi/nobackup/vuw03922/resources/repeat_masker/Libraries/RepeatMaskerLib.h5 names Aves | head
##to print a taxonomic tree including the given clade and optionally ancestors
##    and/or descendants, with the number of repeats indicated at each level of
##    the hierarchy. With the 'totals' format, prints the number of matching
##    ancestral and lineage-specific entries.
#famdb.py -i /nesi/nobackup/vuw03922/resources/repeat_masker/Libraries/RepeatMaskerLib.h5 lineage -ad 'Aves'
####END####

#extract repeats
singularity exec dfam-tetools-latest.sif famdb.py -i /nesi/nobackup/vuw03922/resources/repeat_masker/Libraries/RepeatMaskerLib.h5 families -ad --add-reverse-complement Procellariiformes > $REPEAT_DIR/libraries/Procellariiformes_library.fa


#4. merge dfam and repeat modeler database
cat $REPEAT_DIR/libraries/repeat_modeler_db-families.fa $REPEAT_DIR/libraries/Procellariiformes_library.fa > $REPEAT_DIR/libraries/combined_library.fa


#5. run repeat masker - -xsmall changes repeat regions to lower case (i.e. softmasking)
singularity exec dfam-tetools-latest.sif RepeatMasker -pa 32                   \
                                                      -dir $REPEAT_DIR         \
                                                      -xsmall                  \
                                                      -lib $REPEAT_DIR/libraries/combined_library.fa $ASSEMBLY_DIR/$SAMPLE'_nuclear.fasta'
mv $REPEAT_DIR/$SAMPLE'_nuclear.fasta.masked' $REPEAT_DIR/assembly.softmasked.fasta

#6. convert softmasked to hardmasked
sed 's/(acgt)/N/g' $REPEAT_DIR/assembly.softmasked.fasta > $REPEAT_DIR/assembly.hardmasked.fasta

#steps 1&2
#Cluster: mahuika
#Job ID: 42250136
#State: FAILED #only last step failed, same for below
#Cores: 16
#Tasks: 1
#Nodes: 1
#Job Wall-time:   62.3%  2-11:46:07 of 4-00:00:00 time limit
#CPU Efficiency: 139.3%  55-12:16:49 of 39-20:17:52 core-walltime
#Mem Efficiency:   5.6%  11.18 GB of 200.00 GB

#steps 4-6
#Cluster: mahuika
#Job ID: 42578445
#State: FAILED
#Cores: 8
#Tasks: 1
#Nodes: 1
#Job Wall-time:   33.2%  01:39:31 of 05:00:00 time limit
#CPU Efficiency: 143.2%  19:00:26 of 13:16:08 core-walltime
#Mem Efficiency:  10.4%  5.19 GB of 50.00 GB

#ALL
#Cluster: mahuika
#Job ID: 42583130
#State: FAILED 
#Cores: 16
#Tasks: 1
#Nodes: 1
#Job Wall-time:   61.0%  2-10:34:34 of 4-00:00:00 time limit
#CPU Efficiency: 141.0%  55-01:47:23 of 39-01:13:04 core-walltime
#Mem Efficiency:  82.0%  16.41 GB of 20.00 GB

