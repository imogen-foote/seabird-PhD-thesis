#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --time=0-0:05:00
#SBATCH --job-name=fastStructure
#SBATCH -o /nesi/nobackup/vuw03922/stdout/fastStructure.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/fastStructure.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

##load modules
module purge
module load fastStructure/1.0-gimkl-2020a-Python-2.7.18


#activate conda env
#module purge && module load Miniconda3
#source $(conda info --base)/etc/profile.d/conda.sh
#export PYTHONNOUSERSITE=1
#conda activate /nesi/project/vuw03922/faststructure


###run input
PROJECT=$1
SET=$PROJECT'_'$2
FILTER=$3
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/snp/$SET/${SET}_${FILTER}
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter4/fastStructure/${SET}_${FILTER}
mkdir -p $out_dir


######
##determining model complexity
######
mkdir -p ${out_dir}/model_check

for k in 1 2 3 4 5; do
echo "Analysis for K=$k..."	
structure.py K $k --input=${read_dir} --output=${out_dir}/model_check/${SET}_check --cv=5 --format=bed
done

#parse through output from above to determine reasonable value of K
chooseK.py --input=${out_dir}/model_check/${SET}_check


######
##then run with appropriate value of K
######

K=3
structure.py -K $K --input=${read_dir} --output=${out_dir}/${SET}_fastStructure --format=bed

#conda deactivate


#plot
#module purge
#module load  R/4.3.1-gimkl-2022a

#Rscript
#R=/nesi/nobackup/vuw03922/scripts/Chapter4/R_scripts
#admixture_plot=$R/admixture_plot_K$K.R

#name pops correctly
#paste <(cut -f 1-2 ${read_dir}.fam) $out_dir/${SET}_fastStructure.$K.meanQ > $out_dir/${SET}_fastStructure_named.$K.meanQ

#Rscript $admixture_plot --Q $out_dir/${SET}_fastStructure_named.$K.meanQ \
                --out $out_dir


 
distruct.py -K $K					\
	--input=${out_dir}/${SET}_fastStructure		\
	--output=${out_dir}/${SET}_K${K}_distruct.svg

