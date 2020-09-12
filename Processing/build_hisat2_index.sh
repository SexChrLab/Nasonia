#!/bin/bash
#SBATCH --job-name=hisat2_index # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 96:00:00

source activate nasonia_environment

hisat2-build -f /mnt/storage/SAYRES/NASONIA/Heini/Reference/GCF_009193385.1_Nvit_psr_1_genomic.fna /data/storage/SAYRES/NASONIA/Heini/Reference/HISAT2_index/Vitripennis_HISAT2_index
