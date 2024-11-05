#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=0:10:00
#SBATCH --job-name=mashmap
#SBATCH -o /nesi/nobackup/vuw03922/stdout/mashmap.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/mashmap.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

# Synteny analysis using mashmap
##Manual:https://github.com/marbl/MashMap

##Load modules
module purge
module load MashMap/3.0.4-Miniconda3


###Set variables
antipodean=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/mitogenome/blue-45G_q8_h50/mt_graph.gfa
gibsons=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter3/mitogenome/d36_q8_h50/animal_mt.complete.graph1.selected_graph.gfa
out_dir=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/output/Chapter3/synteny/mitogenome/q8_h50
mkdir -p $out_dir
cd $out_dir


#### Compute table ####
##### You'll need to adjust seq_length and perc_identity to your sepecific case - more strict as evolutionary relatedness increases. The trade-off 
##### will be reduing off-target hits/noise while ensuring information is not lost.

	mashmap -r $antipodean \
	-q $gibsons \
	--segLength 50000 \
	--perc_identity 95 \
	--kmer 16


#### Convert space into tab delimited and add header for table ####
	sed -e 's/ /\t/g' mashmap.out > mashmap.out.table
	sed -i '1 i\query\tlength\t0_based_start\tend\tstrand\tref\tlength\tstart\tend\tidentity' mashmap.out.table


#### Plot genome-wide comparison ####

## I run this script in the command line; replace with the path to the generateDotPlot.pl file
	#perl /nfs/home/finnca/scripts/generateDotPlot.pl png large mashmap.out
