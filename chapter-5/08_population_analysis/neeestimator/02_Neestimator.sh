#!/bin/bash
#SWATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=5G
#SBATCH --time=0-1:0:00
#SBATCH --job-name=neestimator
#SBATCH -o /nesi/nobackup/vuw03922/stdout/neestimator.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/neestimator.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

###run input
PROJECT=$1
SET=$PROJECT'_'$2

###load packages
neestimator=/nesi/project/vuw03922/NeEstimator2/NeEstimator2x1.jar

###set paths
GPOP=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/genepop/$SET
out=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/Neestimator/$SET
mkdir -p $out

java -jar $neestimator 
