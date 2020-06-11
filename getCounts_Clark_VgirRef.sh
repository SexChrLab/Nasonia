#!/bin/bash
#SBATCH --job-name=getCounts_snakemake # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=kcolney@asu.edu # send-to address
#SBATCH -t 96:00:00
#SBATCH -n 1

newgrp combinedlab

source activate nasonia_environment

cd /Processing/Clark/RNA/featureCounts
cut -f7 SRR1566022_featurecounts_VgirRef.tsv  | sed 1d > tmp1
cut -f7 SRR1566023_featurecounts_VgirRef.tsv  | sed 1d > tmp2
cut -f7 SRR1566024_featurecounts_VgirRef.tsv  | sed 1d > tmp3
cut -f7 SRR1566028_featurecounts_VgirRef.tsv  | sed 1d > tmp4
cut -f7 SRR1566029_featurecounts_VgirRef.tsv  | sed 1d > tmp5
cut -f7 SRR1566030_featurecounts_VgirRef.tsv  | sed 1d > tmp6
cut -f7 SRR2773794_featurecounts_VgirRef.tsv  | sed 1d > tmp7
cut -f7 SRR2773795_featurecounts_VgirRef.tsv  | sed 1d > tmp8
cut -f7 SRR2773796_featurecounts_VgirRef.tsv  | sed 1d > tmp9
cut -f7 SRR2773797_featurecounts_VgirRef.tsv  | sed 1d > tmp10
cut -f7 SRR2773798_featurecounts_VgirRef.tsv  | sed 1d > tmp11
cut -f7 SRR2773799_featurecounts_VgirRef.tsv  | sed 1d > tmp12

paste tmp1 tmp2 tmp3 tmp4 tmp5 tmp6 tmp7 tmp8 tmp9 tmp10 tmp11 tmp12 > clark_counts_VgirRef.txt
