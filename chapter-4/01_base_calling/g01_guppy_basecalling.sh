#!/bin/bash
#SBATCH --gpus-per-node=A100:1
#SBATCH --partition=gpu,hgx
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=0-8:00:00
#SBATCH --job-name=guppy_basecaller
#SBATCH -o /nesi/nobackup/vuw03922/stdout/guppy_basecaller.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/guppy_basecaller.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

#load modules
module purge
module load ont-guppy-gpu/6.5.7
module load SAMtools/1.16.1-GCC-11.3.0 

#set up
PROJECT=$1
LIB=$2
MODEL=$3 #either dna_r10.4_e8.1_sup.cfg or dna_r9.4.1_e8.1_sup.cfg
read_dir=/nesi/nobackup/vuw03922/projects/$PROJECT/raw_data/ont/$LIB

mkdir $read_dir/fastq/TMP_fast5
cp $read_dir/fast5/fast5_*/*.fast5 $read_dir/fastq/TMP_fast5

#run guppy
guppy_basecaller -i $read_dir/fastq/TMP_fast5    \
                 -r                         \
                 -s $read_dir/fastq    \
                 --device auto              \
                 --config /opt/nesi/CS400_centos7_bdw/ont-guppy-gpu/6.5.7/data/$MODEL  \
                 --trim_adapters            \
                 --disable_qscore_filtering

#concatenate reads
cat $read_dir/fastq/*fastq > $read_dir/fastq/$LIB.fastq
bgzip $read_dir/fastq/$LIB.fastq

#then remove individual fastq files and tmp directory
rm $read_dir/fastq/*.fastq
rm -r $read_dir/fastq/TMP_fast5


#Cluster: mahuika
#Job ID: 39221455
#State: COMPLETED
#Cores: 2
#Tasks: 1
#Nodes: 1
#Job Wall-time:   24.6%  02:27:20 of 10:00:00 time limit
#CPU Efficiency: 118.1%  05:48:08 of 04:54:40 core-walltime **4CPUs requested**
#Mem Efficiency:  39.4%  2.36 GB of 6.00 GB


### Info parameters
#-i [ --input_path ]	Path to input files.
#-r [ --recursive ]		Search for input file recursively.
#-s [ --save_path ]		Path to save output files.
#-x [ --device ]		Specify GPU device: 'auto', or 'cuda:<device_id>'.
#-c [ --config ]		Configuration file for application.
#--compress_fastq		Compress fastq output files with gzip.
#--num_callers			Number of parallel basecallers to create.
#--duplex_pairing_mode	Read pairing mode to use for duplex basecalling, may be either 'from_1d_summary' or 'from_pair_list'
#--duplex_pairing_file	Input filename to use for duplex pairing

### Flow cell and kit info
#R10.4 (FLO-MIN112, FLO-PRO112) E8.1 (Kit 12) 5mC, 5hmC_5mC

### Chosing the right config file
#I use dna_r10.4_e8.1_sup.cfg
#<strand_type>_<pore_type>_<enzyme_type>_[modbases_specifier]_<model_type>_[instrument_type].cfg
#strand_type		: This will be either the string "dna" or "rna", depending on the type of sequencing being performed.
#pore_type			: The pore the basecalling model was trained for, indicated by the letter "r" followed by a version number. For example:"r9.4.1" or "r10.4".
#enzyme_type		: The enzyme motor the model was trained for. This will either be the letter "e" followed by a version number, or a number indicating the enzyme speed, followed by "bps". For example: "e8.1" or "450bps".
#modbase_specifier	: Optional. If specified, indicates that modified base detection will be performed. 
#This will be the string "modbases_" followed by an indicator of the modification supported, such as "5mc_cg" or "5hmc_5mc_cg".
#model_type			: The type of basecalling model to use, depending on whether you want optimal basecalling speed or accuracy. See below.
#instrument_type	: Optional. If this is not specified, then the configuration is targeted to a GridION device or a PC. 
#The strings "mk1c" or "prom" are used to indicate that the configuration parameters and model are optimised for the MinION Mk1C or PromethION devices, respectively. 
#Note that if the kit and flow cell are specified on the command-line instead of a specific config file, then the config file chosen will be one without an instrument type specified.
#The model types are:
	#sup: Super-accurate basecalling.
	#hac: High accuracy basecalling. These are the configurations that will be selected when a kit and flow cell are specified on the
	#command-line instead of a specific config file.
	#fast: Fast basecalling.
	#sketch: Sketch basecalling. This is primarily for use with adaptive sampling on the MinION Mk1C device to minimise latency.
#For example, to basecall data generated with the R10.4 pore and the E8.1 enzyme, using the Fast CRF model: guppy_basecaller -c dna_r10.4_e8.1_fast.cfg [..]

### Chosing gpu device
#"cuda:0" Use the first GPU in the system, no memory limit
#"cuda:0,1" Use the first two GPUs in the system, no memory limit
#"cuda:0 cuda:1" Same as cuda:0,1
#"cuda:all:100%" Use all GPUs in the system, no memory limit
#"cuda:1,2:50%" Use the second and third GPU in the system, and use only up to half of the GPU memory of each GPU
#"cuda:0 cuda:1,2:8G" Use the first three GPUs in the system. Use a maximum of 8 GiB on each of GPUs 1 and 2.
#"auto" Same as cuda:0
