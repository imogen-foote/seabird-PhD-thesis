#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH -a 1-43
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=2-0:00:00
#SBATCH --job-name=paleomix_bwa_mem
#SBATCH -o /nesi/nobackup/vuw03922/stdout/paleomix.%A_%a.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/paleomix.%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#activate conda env
module purge && module load Miniconda3 
source $(conda info --base)/etc/profile.d/conda.sh 
export PYTHONNOUSERSITE=1
conda activate /nesi/project/vuw03922/paleomix
 
PROJECT=$1
N=${SLURM_ARRAY_TASK_ID}


sample_list=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/paleomix/samples.txt
sample=$( head -n $N $sample_list | tail -n 1 | cut -f 1 )

read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/paleomix/${sample}

#bam pipeline
echo "running PALEOMIX"
cd $read_dir
paleomix bam run *.yaml

#cpt=2, mem=25G, total time: 1h 19m. Max mem used=18.57G

conda deactivate



#Cluster: mahuika
#Job ID: 43223568
#Array Job ID: 43223568_2
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:   35.7%  17:09:06 of 2-00:00:00 time limit
#CPU Efficiency: 178.6%  1-06:38:05 of 17:09:06 core-walltime
#Mem Efficiency:  29.6%  8.88 GB of 30.00 GB

