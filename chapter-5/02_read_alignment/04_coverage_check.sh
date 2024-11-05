#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=500MB
#SBATCH --time=1-0:0:00
#SBATCH --job-name=coverage_check
#SBATCH -o /nesi/nobackup/vuw03922/stdout/coverage_check.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/coverage_check.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

# Load modules
module purge
module load SAMtools/1.16.1-GCC-11.3.0
module load gnuplot/5.4.2-GCC-7.4.0

# Define variables and paths
PROJECT=$1

bam_list=/nesi/nobackup/vuw03922/projects/$PROJECT/resources/sample_info/$PROJECT'_bam.list'
out_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/data/Chapter4/samtools_coverage
mkdir -p "$out_dir"
out_file="$out_dir/coverage.txt"
out_plot="$out_dir/coverage.png"
dev_file="$out_dir/t_test.txt"

# Create an array to store individual coverage files
declare -a coverage_files


# Loop through each BAM file path in the text file
while IFS= read -r bam_path
do
    # Extract sample name from the path
    sample=$(basename "$(dirname "$bam_path")")

    # Define output coverage file for the individual
    coverage_file="$out_dir/${sample}_coverage.txt"

    # Calculate coverage using samtools depth
    samtools depth "$bam_path" | awk '{sum += $3} END {print "Average Coverage " $sample ": ", sum/NR; print "Stdev = ", sqrt(sumsq/NR - (sum/NR)**2)}' > "$coverage_file"

    # Add the coverage file to the array 
    coverage_files+=("$coverage_file")

done < "$bam_list"

# Concatenate all individual coverage files into a single file
cat "${coverage_files[@]}" > "$out_file"

# Plot the coverage using gnuplot
#gnuplot <<- GNU_END
#    set terminal png
#    set output "$out_plot"
#    plot "$out_file" using 2 with lines title "Coverage"
#GNU_END

# Clean up individual coverage files
rm -f "${coverage_files[@]}"

echo "Coverage calculation and plotting completed!"

## Check deviation from mean for any samples
# Calculate the mean coverage
#mean_coverage=$(awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' "$out_file")

# Perform t-test using awk
#awk -v mean_coverage="$mean_coverage" '
#    { sum += ($2 - mean_coverage)^2; count++ }
#    END {
#        variance = sum / count;
#        t_statistic = (mean_coverage - $2) / sqrt(variance / count);
#        p_value = 2 * atan2(sqrt(variance / count), abs(t_statistic));
#        print "Mean Coverage:", mean_coverage;
#        print "T-Statistic:", t_statistic;
#        print "P-Value:", p_value;
#    }' "$out_file" > "$dev_file"

#echo "T-test completed."
