#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=1-0:0:00
#SBATCH --job-name=map_check
#SBATCH -o /nesi/nobackup/vuw03922/stdout/map_check.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/map_check.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load gnuplot/5.4.2-GCC-7.4.0
 
PROJECT=$1
bam_list=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_bam.list'
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/samtools_mapcheck
mkdir -p $out_dir

# Declare arrays to store sample names and corresponding files
sample_names=()
mapping_stats_files=()
proper_pair_stats_files=()

# Loop over all samples in the list
while read -r current_bam; do
    # Extract the sample name from the BAM file name
    sample_name=$(basename "$(dirname "$current_bam")")
    sample_names+=($sample_name)  # Store sample names in the array

    # Define output files
    mapping_stats_file=$out_dir/${sample_name}_mapping_stats.txt
    proper_pair_stats_file=$out_dir/${sample_name}_proper_pair_stats.txt

    # Use samtools to get mapping and proper pair rates
    samtools stats $current_bam | grep '^SN' > $mapping_stats_file
    samtools view -f 2 -c $current_bam > $proper_pair_stats_file

    mapping_stats_files+=($mapping_stats_file)
    proper_pair_stats_files+=($proper_pair_stats_file)

    # Print summary to a summary file
    echo -e "Sample: $sample_name" >> $out_dir/summary.txt
    echo -e "Mapping Rates:" >> $out_dir/summary.txt
    cat $mapping_stats_file >> $out_dir/summary.txt
    echo -e "\nProper Pair Count:" >> $out_dir/summary.txt
    cat $proper_pair_stats_file >> $out_dir/summary.txt
    echo -e "\n" >> $out_dir/summary.txt

done < "$bam_list"

# Plotting using gnuplot for mapping rates
#gnuplot << EOF
#set terminal png size 800,400
#set output "$out_dir/all_samples_mapping_plot.png"
#set title "Mapping Rates - All Samples"
#set xlabel "Category"
#set ylabel "Percentage"
#set yrange [0:100]
#set xtic rotate

# Plot all samples on one graph
#plot for [i=1:${#sample_names[*]}] "${mapping_stats_files[i-1]}" using 3:xtic(2) with bars title "${sample_names[i-1]}" linecolor rgb "blue"
#EOF

# Plotting using gnuplot for proper pair counts
#gnuplot << EOF
#set terminal png size 800,400
#set output "$out_dir/all_samples_proper_pair_plot.png"
#set title "Proper Pair Counts - All Samples"
#set xlabel "Category"
#set ylabel "Count"
#set yrange [0:*]
#set xtic rotate

# Plot all samples on one graph
#plot for [i=1:${#sample_names[*]}] "${proper_pair_stats_files[i-1]}" using 1:xtic(2) with boxes title "${sample_names[i-1]}" linecolor rgb "green"
#EOF
