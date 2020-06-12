# Nasonia
Lack of parent-of-origin effects in Nasonia jewel wasp - a replication and reproduction study

## Background
Parent and species of origin allele specific and differential expression analysis workflow. A replication and reproduction study using RNA and DNA from Nasonia jewel wasp a haplodiploidy species. In diploid cells, the paternal and maternal alleles are on average equally expressed. There are exceptions from this is which a small number of genes express exclusively the maternal or paternal allele copy, known as genomic imprinting. Clark et al. 2016 found no parent-of-origin effects in the hybrids of closely related Nasonia vitripennis and N. giuralti jewel wasp, suggesting a lack of epigenetic reprogramming during embryogenesis in these species. Here, we have reproduced and replicated these findings using the previously published RNA & DNA sequence data from 11 samples as well as a newly generated RNA data from 12 samples of the same species hybrids. Our results from both datasets demonstrated a species-of-origin effect.

<p align="center">
  <img src="https://github.com/SexChrLab/Nasonia/blob/master/Plotting/Figure1.png" width="500"/>
</p>

Figure 1. A) A schematic illustration of the reciprocal F1 crosses. B) Overview of the data processing and analysis workflow.


## Preprint
Please see our preprint for more information:

*Lack of parent-of-origin effects in Nasonia jewel wasp - a replication and reproduction study* 2020. Olney KC*; Natri HM*; Underwood A; Gibson J; Gadau J; Wilson MA.  

\* indicates equal contribution.

If you use Nasonia in a classroom or discuss/correct this work in a manuscript, please cite this preprint.

## Contents:
1. Process data - Download and process files and identify sites that are fixed and different between the inbreed lines of N. Vit and N. Gir 

2. Differential expression - Differential expression for all pairwise comparisons 

3. Allele specific expression - Identify genes with species or parent-of-origin effects

### 1. Process data

The Wilson and Clark data were processed on a HPC cluster using SLURM and snakemake. To create a conda environment with the tools used in data processing, install conda (e.g. Miniconda) and run "conda env create -f nasonia_environment.yaml".

###### QC and read trimming with FastQC, MultiQC, and Trimmomatic

The quality of the FASTQ files was assessed using FastQC and MultiQC. Reads were trimmed for quality and post-trimming quality was assesses again. To run these steps, use the script run_fastqc_trimmomatic.sh and the snakefile fastqc_trimmomatic.snakefile.

###### Read mapping and gene expression quantification

The Nasonia vitripennis reference genome and gene annotation file were downloaded from NCBI (https://www.ncbi.nlm.nih.gov/genome/449?genome_assembly_id=716919). FASTA index file and a sequence dictionary file must be present in the reference directory. These files can be produced with samtools and Picard:

```
samtools faidx GCF_009193385.2_Nvit_psr_1_genomic.fa
```

```
picard CreateSequenceDictionary R=GCF_009193385.1_Nvit_psr_1_genomic.fa O=GCF_009193385.2_Nvit_psr_1_genomic.dict
```

Trimmed reads were mapped to the reference with HISAT2. The HISAT2 index was created using the script build_hisat_index.sh. Gene level counts were quantified using featureCounts. To start the snakemake pipeline for read mapping and gene expression quantification, run run_mapping_quantification_snakemake.sh.

###### BAM processing, variant calling, and ASE quantification

BAM files were processed for variant calling using picard. To run the processing, use the script run_bam_processing.sh and the snakefile process_bam.snakefile. Variants were called using GATK HaplotypeCaller. Allele-specific expression levels were quantified using GATK ASEReadCounter. To run these steps, use the scripts run_gatk_haplotypecaller_asereadcounter.sh and gatk_haplotypecaller_asereadcounter.snakefile.

Read mapping and all downstream processing steps were repeated with the custom generated N. giraulti reference genome.


### 2. Fixed differences, pseudo N.giraulti reference genome

We created a custom Nasonia giraulti reference genome by substituting sites with fixed allele differences between the two inbred species with the N. giraulti allele using GATK’s FastaAlternateReferenceMaker (McKenna et al. 2010). A site was considered to be fixed and different if it was homozygous for the N. vitripennis reference allele among all three of the biological VV samples and homozygous alternate among all three of the biological GG samples. All samples were aligned to both the N.vit and the pseudo N.gir reference genomes as described in step 1. 

To make a pseudo N.giraulti reference genome you first need to have called variants when samples were aligned to the N.vit reference genome from step 1. Then you can run the script run_fixedDifferences.sh which will call a set of python scripts for identifying sits that are fixed and different between VV and GG. After this step has finished, there will be a text file in the FixedDifferences folder that contian the fixed sites infomation, Sayres_Clark_intersect.txt. Now run the script run_FastaAlternateReference.sh which will then make the pseudo N.giraulti reference genome with fixed differences information contained in the Sayres_Clark_intersect.txt file. Re-process the samples as described in step 1 using the new pseudo N.giraulti reference genome for alignment. 

### 3. Differential expression analysis

Gene quanitification from  featureCounts for each dataset was put into a matrix and then read into into R for calling differential expression between species. Differential expression was preformed using limma/voom. 

To run these steps, use the scripts getCounts_Clark.sh and getCounts_Wilson.sh for when the samples were aligned to the N.vit references genome and run the scripts getCounts_Clark_VgirRef.sh and getCounts_Wilson_VgirRef.sh for when the samples were aligned to the pseudo N.gir references genome. This will create the count matrix needed to run differential expression. 

To run differntial expression, use the scripts: Wilson_AvgREFs_limmaVoom.r and Clark_AvgREFs_limmaVoom.r  

For an R markdown of the differntial expression results visit: 
http://rpubs.com/olneykimberly/626815 or /Plotting/Nasonia_Wilson_data.pdf and /Plotting/Nasonia_Clark_data.pdf


### 4. Allele specific expression analysis

For each sample, allele-specific expression (ASE) analysis was done on the reads aligned to the N. vitripennis reference, and on the reads aligned to the pseudo N. giraulti reference for comparison. Counts of reads covering biallelic heterozygous SNP sites were obtained with ASEReadCounter in GATK version 3.8 (Castel et al. 2015). Parameters ensuring adequate coverage with a minimum mapping quality of 10, minimum base quality of 2, and a minimum depth of 30 were used. Only sites with a fixed difference between inbred VV and GG for both Clark and Wilson datasets were used for downstream analysis of allele-specific expression.

To get ASE run first 

ASE_AvgREFs.r

<p align="center">
  <img src="https://github.com/SexChrLab/Nasonia/blob/master/Plotting/Figure4.png" width="900"/>
</p>

Figure 4. Scatterplots of the expression of the N. vitripennis alleles in the two reciprocal hybrids, VG (x-axis) and GV (y-axis), in the Clark (A) and Wilson (B) datasets. Genes with at least two informative SNPs with a min depth of 30 were used (Clark = 6,377, Wilson  = 7,164). Genes exhibiting a significant difference in allelic bias between the hybrids (Fisher’s exact test, FDR-adj. p<0.01) are highlighted in red. Paternally imprinted genes are expected to appear in the upper left corner (light blue box), and maternally imprinted genes in the lower right corner (light pink box). Histograms of the N. vitripennis allele expression are shown for VG (blue) and GV (pink). 


### Publicly available tools used in this analysis

These tools are publicly available and we ask that if you use this workflow to cite the tools used listed in the table below.

Tool | usage | citation
--- | --- |  ---
Trimmomatic | Trim RNA-sequences for quality | Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30: 2114–2120.
HISAT | RNAseq read aligner | Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12: 357–360.
bamtools | analyzing and processing BAM files | Barnett DW, Garrison EK, Quinlan AR, Strömberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011;27: 1691–1692.
FeatureCounts | obtain raw transcriptome counts| Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30: 923–930.
Limma/voom | differenital expression analysis | Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014;15: R29.
GATK | variant calling and reference alternate maker | McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297‐1303. doi:10.1101/gr.107524.110
Picard toolkit | add read groups and mark duplicates | https://github.com/broadinstitute/picard


## Group Members
Name | email | github ID
--- | --- |  ---
Kimberly Olney | olneykimberly@gmail.com | @olneykimberly
Heini M. Natri | heini.natri@gmail.com |@heinin
Melissa A. Wilson | melissa.wilsonsayres@asu.edu | @mwilsonsayres
