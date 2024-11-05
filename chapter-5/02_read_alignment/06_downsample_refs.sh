#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=1:30:00
#SBATCH --job-name=downsample
#SBATCH -o /nesi/nobackup/vuw03922/stdout/downsample.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/downsample.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz


#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
 
PROJECT=$1
IND=$2 #ANT_SA52 or GIB_AA08

bam=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/paleomix/$IND


#calculate initial coverage
echo "calculating initial coverage stats..."
samtools depth $bam/*'_nuclear_clipped.bam' | awk '{sum += $3} END {print "Average Coverage: ", sum/NR}'

#downsample 
echo "downsampling..."
samtools view -s 0.5 -b -o $bam/downsampled_nuclear_clipped.bam $bam/*'_nuclear_clipped.bam' 

#index new bam file
echo "indexing downsampled bam..."
samtools index $bam/downsampled_nuclear_clipped.bam

#verify coverage
echo "calculating post-downsampling coverage stats..."
samtools depth $bam/downsampled_nuclear_clipped.bam | awk '{sum += $3} END {print "Average Coverage: ", sum/NR}'

#rename original bam file
#mv $bam/*'_nuclear_clipped.bam' $bam/*'_nuclear_clipped_FULL.bam'
#mv $bam/downsampled.bam $bam/*'_nuclear_clipped.bam'


#Cluster: mahuika
#Job ID: 43669588
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   65.2%  00:58:43 of 01:30:00 time limit
#CPU Efficiency: 131.7%  01:17:20 of 00:58:43 core-walltime
#Mem Efficiency:   1.0%  9.87 MB of 1.00 GB

