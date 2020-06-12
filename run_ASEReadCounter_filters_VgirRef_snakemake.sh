#!/bin/bash
#SBATCH --job-name=filter_ASE_snakemake # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=kcolney@asu.edu # send-to address
#SBATCH -t 09:00:00
#SBATCH -n 1

newgrp combinedlab

source activate PopInf

module load java/8u92

export PERL5LIB=/packages/6x/vcftools/0.1.12b/lib/perl5/site_perl

#snakemake --snakefile ASEReadCounter_filters_FD_bed.snakefile -j 20 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 8 -t 96:00:00"
#snakemake --snakefile Snakefile -j 40  --latency-wait 15 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 8 -t 96:00:00"
#ASE_VgirRef_minDepth100_step1.snakefile
snakemake --snakefile ASE_VgirRef_minDepth30_step2.snakefile -j 30 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 8 -t 04:00:00"

#snakemake --snakefile ASE_VgirRef_minDepth30.snakefile -j 40  --latency-wait 15 --rerun-incomplete --cluster "sbatch -n 1 --nodes 1 -c 8 -t 96:00:00"

