#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=0:10:00
#SBATCH --job-name=GONE
#SBATCH -o /nesi/nobackup/vuw03922/stdout/GONE.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/GONE.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##Copy files from github##
#cd /nesi/project/vuw03922
#git clone https://github.com/esrud/GONE.git
GONE=./script_GONE.sh

#set variables
PROJECT=$1
SET=$PROJECT'_'$2
FILTER=$3
pop=$4

#define paths
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/${SET}_${FILTER}_${pop}
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/GONE/${SET}_${FILTER}/${pop}
mkdir -p $out_dir

#copy files to this directory for running
cp ${read_dir}.map ${read_dir}.ped .

#run GONE
bash $GONE ${SET}_${FILTER}_${pop}

#remove map and ped files
rm ${SET}_${FILTER}_${pop}.map ${SET}_${FILTER}_${pop}.ped

#move all output files to out_dir
mv OUTPUT* $out_dir
mv Output* $out_dir
mv TEMPORARY_FILES/ $out_dir
mv outfileHWD $out_dir
mv timefile $out_dir
mv seedfile $out_dir 



#Cluster: mahuika
#Job ID: 45635625
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:    0.2%  00:06:50 of 3-00:00:00 time limit
#CPU Efficiency: 166.7%  00:22:47 of 00:13:40 core-walltime
#Mem Efficiency:   0.3%  53.82 MB of 20.00 GB

