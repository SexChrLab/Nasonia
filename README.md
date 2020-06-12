# Nasonia
Lack of parent-of-origin effects in Nasonia jewel wasp - a replication and reproduction study

## Background
Parent and species of origin allele specific and differential expression analysis workflow. A replication and reproduction study using RNA and DNA from Nasonia jewel wasp a haplodiploidy species. In diploid cells, the paternal and maternal alleles are on average equally expressed. There are exceptions from this is which a small number of genes express exclusively the maternal or paternal allele copy, known as genomic imprinting. Clark et al. 2016 found no parent-of-origin effects in the hybrids of closely related Nasonia vitripennis and N. giuralti jewel wasp, suggesting a lack of epigenetic reprogramming during embryogenesis in these species. Here, we have reproduced and replicated these findings using the previously published RNA & DNA sequence data from 11 samples as well as a newly generated RNA data from 12 samples of the same species hybrids. Our results from both datasets demonstrated a species-of-origin effect.

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
samtools faidx GCF_009193385.1_Nvit_psr_1_genomic.fa
```

```
picard CreateSequenceDictionary R=GCF_009193385.1_Nvit_psr_1_genomic.fa O=GCF_009193385.1_Nvit_psr_1_genomic.dict
```

Trimmed reads were mapped to the reference with HISAT2. The HISAT2 index was created using the script build_hisat_index.sh. Gene level counts were quantified using featureCounts. To start the snakemake pipeline for read mapping and gene expression quantification, run run_mapping_quantification_snakemake.sh.

###### BAM processing, variant calling, and ASE quantification

BAM files were processed for variant calling using picard. To run the processing, use the script run_bam_processing.sh and the snakefile process_bam.snakefile. Variants were called using GATK HaplotypeCaller. Allele-specific expression levels were quantified using GATK ASEReadCounter. To run these steps, use the scripts run_gatk_haplotypecaller_asereadcounter.sh and gatk_haplotypecaller_asereadcounter.snakefile.

Read mapping and all downstream processing steps were repeated with the custom generated N. giraulti reference genome.


### 2. Fixed differences, pseudo N.giraulti reference genome


### 3. Differential expression analysis

Gene quanitification from  featureCounts for each dataset was put into a matrix and then read into into R for calling differential expression between species. Differential expression was preformed using limma/voom. 

To run these steps, use the scripts getCounts_Clark.sh and getCounts_Wilson.sh for when the samples were aligned to the N.vit references genome and run the scripts getCounts_Clark_VgirRef.sh and getCounts_Wilson_VgirRef.sh for when the samples were aligned to the pseudo N.gir references genome. This will create the count matrix needed to run differential expression. 

To run differntial expression, use the scripts 

For an R markdown of the differntial expression results visit: 
http://rpubs.com/olneykimberly/626815 



### 4. Allele specific expression analysis




These tools are publicly available and we ask that if you use this workflow to cite the tools used listed in the table below.

### Publicly available tools used in this analysis
Tool | usage | citation
--- | --- |  ---
Trimmomatic | Trim RNA-sequences for quality | Bolger AM, Lohse M, Usadel B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics. 2014;30: 2114–2120.
HISAT | RNAseq read aligner | Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12: 357–360.
bamtools | analyzing and processing BAM files | Barnett DW, Garrison EK, Quinlan AR, Strömberg MP, Marth GT. BamTools: a C++ API and toolkit for analyzing and managing BAM files. Bioinformatics. 2011;27: 1691–1692.
FeatureCounts | obtain raw transcriptome counts| Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2014;30: 923–930.
Limma/voom | differenital expression analysis | Law CW, Chen Y, Shi W, Smyth GK. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol. 2014;15: R29.
GATK | | 
Picard toolkit | | 


## Group Members
Name | email | github ID
--- | --- |  ---
Kimberly Olney | olneykimberly@gmail.com | @olneykimberly
Heini M. Natri | heini.natri@gmail.com |@heinin
Melissa A. Wilson | melissa.wilsonsayres@asu.edu | @mwilsonsayres
