#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=8:0:00
#SBATCH --job-name=mummer
#SBATCH -o /nesi/nobackup/vuw03922/stdout/mummer.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/mummer.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##Load modules
module purge
module load MUMmer/4.0.0rc1-GCCcore-11.3.0


###Set variables
antipodean=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/blue-45G_q8_h50_masked/blue-45G_softmasked_renamed.fasta
gibsons=/nesi/nobackup/vuw03922/projects/GibsonsAlbatross/data/Chapter3/assemblies/flye/5_repeat_masking/d36_q8_h50_masked/d36_softmasked_renamed.fasta
out_dir=/nesi/nobackup/vuw03922/projects/AntipodeanAlbatross/output/Chapter3/synteny/MUMmer/antipodean_gibsons_q8_h50_softmasked
mkdir -p $out_dir
cd $out_dir

#set the output prefix for MUMmer result
out_prefix=antipodean_gibson

#run NUCmer to perform the genome comparison
nucmer --maxmatch --prefix=$out_prefix $antipodean $gibsons -t 32

#use show-coords to display the coordinates of the alignments
show-coords -rcl $out_prefix.delta > $out_prefix.coords

#generate a graphical representation of the alignment
mummerplot --png --layout --prefix=$out_prefix $out_prefix.delta

echo "MUMmer comparison completed."
