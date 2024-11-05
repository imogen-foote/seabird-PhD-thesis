#!/bin/bash
#SBATCH --account=vuw03922
#SBATCH --cpus-per-task=2
#SBATCH --mem=1G
#SBATCH --time=0:10:00
#SBATCH --job-name=dgenies
#SBATCH -o /nesi/nobackup/vuw03922/stdout/dgenies.%j.out
#SBATCH -e /nesi/nobackup/vuw03922/stdout/dgenies.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=imogen.foote@vuw.ac.nz

