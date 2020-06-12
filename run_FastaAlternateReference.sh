#!/bin/bash
#SBATCH -n 6
#SBATCH --job-name=FastaAlternateReferenceMaker
#SBATCH -t 1-0:0
#SBATCH -A mwilsons
#SBATCH -o slurm.FastaAlternateReferenceMaker.out
#SBATCH -e slurm.FastaAlternateReferenceMaker.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=kcolney@asu.edu
#--------------------------
module load java/latest
module load samtools/1.9

source activate XY_RNA-Seq

cd Reference/

java -Xmx16g -jar /mnt/storage/SAYRES/NASONIA/01_tools/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R GCF_009193385.1_Nvit_psr_1_genomic.fa -o GCF_009193385.1_pseudoNgir_psr_1_genomic.fa -V /mnt/storage/SAYRES/NASONIA/Heini/FixedDifferences/Sayres_FD

samtools faidx GCF_009193385.1_pseudoNgir_psr_1_genomic.fa

java -jar picard.jar CreateSequenceDictionary REFERENCE=GCF_009193385.1_pseudoNgir_psr_1_genomic.fa OUTPUT=GCF_009193385.1_pseudoNgir_psr_1_genomic.dict

cut -f1 GCF_009193385.1_pseudoNgir_psr_1_genomic.fa.fai > pseduoNgir_chrIDs.txt
cut -f1 GCF_009193385.1_Nvit_psr_1_genomic.fa.fai > Nvit_chrIDs.txt

paste Nvit_chrIDs.txt pseduoNgir_chrIDs.txt > alias_chrIDs.txt

