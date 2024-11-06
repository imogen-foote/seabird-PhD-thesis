#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0-0:10:00
#SBATCH --job-name=circle_plot
#SBATCH -o /nesi/nobackup/vuw03922/stdout/circle_plot.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/circle_plot.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

module load Singularity/3.11.3

PROJECT=$1
SAMPLE=$2
FILTER=$3

ASSEMBLY_FILE=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter3/assemblies/flye/4_nuclear_only/$SAMPLE'_'$FILTER'_nuclear'/$SAMPLE'_nuclear.fasta'

out=/nesi/nobackup/vuw03922/projects/$PROJECT/output/Chapter3/flye_assembly/4_nuclear_only/$SAMPLE'_'$FILTER'_nuclear'/$SAMPLE'_'$FILTER'_nuclear.circle_plot'
mkdir -p $out

perl=/nesi/nobackup/vuw03922/scripts/Chapter3/04_assembly_stats/circle_plot/pl/asm2stats.pl
#perl=$SCRATCH/scripts/Chapter3/04_assembly_stats/circle_plot/pl/asm2stats.minmaxgc.pl

#copy template folder to output
cp -r /nesi/nobackup/vuw03922/scripts/Chapter3/04_assembly_stats/circle_plot/ $out

#create json file
echo "var ${PROJECT} = " > $out/circle_plot/json/${PROJECT}.json
perl $perl $ASSEMBLY_FILE >> $out/circle_plot/json/${PROJECT}.json
echo "localStorage.setItem('${PROJECT}',JSON.stringify(${PROJECT}))" >> $out/circle_plot/json/${PROJECT}.json

#add json to html
sed -i s"%<!--add_jsons_here-->%  <!--add_jsons_here-->\n  <script type=\"text/javascript\" src=\"json/${PROJECT}.json\"></script>%"g $out/circle_plot/assembly-stats.html



#Cluster: mahuika
#Job ID: 40893545
#State: COMPLETED
#Cores: 1
#Tasks: 1
#Nodes: 1
#Job Wall-time:    5.0%  00:03:01 of 01:00:00 time limit
#CPU Efficiency: 100.0%  00:03:01 of 00:03:01 core-walltime
#Mem Efficiency:  40.4%  4.04 GB of 10.00 GB
