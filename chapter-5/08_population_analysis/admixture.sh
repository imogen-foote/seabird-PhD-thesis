#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=4
#SBATCH --mem=500M
#SBATCH --time=0-2:0:00
#SBATCH --job-name=admixture
#SBATCH -o /nesi/nobackup/vuw03922/stdout/admixture.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/admixture.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##load modules
module purge
module load R/4.3.1-gimkl-2022a
module load Miniconda3
source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#activate conda env
conda activate /nesi/project/vuw03922/admixture


###run input
PROJECT=$1
SET=$PROJECT'_'$2
FILTER=$3
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/${SET}_${FILTER}
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/admixture/${SET}_${FILTER}
mkdir -p $out_dir


######
#first run across a range of K to determine the best value of K (lowest cross-validation error value)
######
#mkdir $out_dir/CV_check

#for K in 1 2 3 4 5; do
#admixture --cv ${read_dir}.bed $K | tee log${K}.out > $out_dir/CV_check/${SET}.K${K}.out

#mv ${SET}*.P $out_dir/CV_check/
#mv ${SET}*.Q $out_dir/CV_check/
#mv log* $out_dir/CV_check/

#done

#awk '/CV/ {print $3,$4}' $out_dir/CV_check/*out | cut -c 4,7-20 > $out_dir/CV_check/${SET}.cv.error



######
##then run with appropriate value of K
######

K=2
#admixture -C 10000 ${read_dir}.bed $K -j4 > $out_dir/${SET}.K${K}.out

#mv ${SET}*.P $out_dir/
#mv ${SET}*.Q $out_dir/
#mv log* $out_dir/

#name pops correctly
#paste <(cut -f 1-2 ${read_dir}.fam) $out_dir/${SET}_neutral.$K.Q > $out_dir/${SET}_neutral_named.$K.Q

#conda deactivate

#Rscript
R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
admixture_plot=$R/admixture_plot_K$K.R


Rscript $admixture_plot --Q $out_dir/${SET}_neutral_named.$K.Q \
		--out $out_dir



#Cluster: mahuika
#Job ID: 44128906
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   81.1%  00:24:19 of 00:30:00 time limit
#CPU Efficiency: 198.7%  01:36:37 of 00:48:38 core-walltime
#Mem Efficiency:   0.6%  6.18 MB of 1.00 GB

