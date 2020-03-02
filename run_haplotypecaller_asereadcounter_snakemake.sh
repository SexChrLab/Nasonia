#!/bin/bash
#SBATCH --job-name=GATK_HaplotypeCaller_ASE # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=hnatri@asu.edu # send-to address
#SBATCH -t 96:00:00
#SBATCH -n 1

newgrp combinedlab

source activate nasonia_environment

module load java/8u92

export PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/perl5/site_perl

snakemake --snakefile gatk_haplotypecaller_asereadcounter.snakefile -j 20 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 8 -t 96:00:00"