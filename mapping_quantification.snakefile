# Workflow for aligning RNA sequence data to a reference with HISAT2,
# sorting and indexing BAM files with Samtools, and quantifying
# read counts with Subread featureCounts.

configfile: "Nasonia.config.json"

# Tools
HISAT2 = "hisat2"
SAMTOOLS = "samtools"
FEATURECOUNTS = "featureCounts"

# Reference genome files
VITRIPENNIS_GFF = config["Vitripennis_GFF"]
VITRIPENNIS_HISAT2_INDEX = config["Vitripennis_HISAT2_index"]

# Directories
Sayres_RNA_trimmed_FASTQ = config["Sayres_RNA_trimmed_FASTQ"]
Clark_RNA_trimmed_FASTQ = config["Clark_RNA_trimmed_FASTQ"]
Clark_exome_trimmed_FASTQ = config["Clark_exome_trimmed_FASTQ"]

Sayres_RNA_SAM_AL_DIR = config["Sayres_RNA_SAM"]
Sayres_RNA_BAM_AL_DIR = config["Sayres_RNA_BAM"]
Sayres_RNA_SORTED_BAM_AL_DIR = config["Sayres_RNA_processed_BAM"]
Sayres_RNA_FEATURECOUNTS_DIR = config["Sayres_RNA_featureCounts"]

Clark_RNA_SAM_AL_DIR = config["Clark_RNA_SAM"]
Clark_RNA_BAM_AL_DIR = config["Clark_RNA_BAM"]
Clark_RNA_SORTED_BAM_AL_DIR = config["Clark_RNA_processed_BAM"]
Clark_RNA_FEATURECOUNTS_DIR = config["Clark_RNA_featureCounts"]

Clark_exome_SAM_AL_DIR = config["Clark_exome_SAM"]
Clark_exome_BAM_AL_DIR = config["Clark_exome_BAM"]
Clark_exome_SORTED_BAM_AL_DIR = config["Clark_exome_processed_BAM"]

# Samples
SAYRES_SAMPLES = config["SAYRES_SAMPLES"]
CLARK_RNA_SAMPLES= config["CLARK_RNA_SAMPLES"]
CLARK_EXOME_SAMPLES = config["CLARK_EXOME_SAMPLES"]

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(Sayres_RNA_SAM_AL_DIR + "{sample}.sam", Sayres_RNA_SAM_AL_DIR=Sayres_RNA_SAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_BAM_AL_DIR + "{sample}.bam", Sayres_RNA_BAM_AL_DIR=Sayres_RNA_BAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam", Sayres_RNA_SORTED_BAM_AL_DIR=Sayres_RNA_SORTED_BAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam.bai", Sayres_RNA_SORTED_BAM_AL_DIR=Sayres_RNA_SORTED_BAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_FEATURECOUNTS_DIR + "{sample}_featurecounts.tsv", Sayres_RNA_FEATURECOUNTS_DIR=Sayres_RNA_FEATURECOUNTS_DIR, sample=SAYRES_SAMPLES),

		expand(Clark_RNA_SAM_AL_DIR + "{sample}.sam", Clark_RNA_SAM_AL_DIR=Clark_RNA_SAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_BAM_AL_DIR + "{sample}.bam", Clark_RNA_BAM_AL_DIR=Clark_RNA_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam", Clark_RNA_SORTED_BAM_AL_DIR=Clark_RNA_SORTED_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam.bai", Clark_RNA_SORTED_BAM_AL_DIR=Clark_RNA_SORTED_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_FEATURECOUNTS_DIR + "{sample}_featurecounts.tsv", Clark_RNA_FEATURECOUNTS_DIR=Clark_RNA_FEATURECOUNTS_DIR, sample=CLARK_RNA_SAMPLES),

		expand(Clark_exome_SAM_AL_DIR + "{sample}.sam", Clark_exome_SAM_AL_DIR=Clark_exome_SAM_AL_DIR, sample=CLARK_EXOME_SAMPLES),
		expand(Clark_exome_BAM_AL_DIR + "{sample}.bam", Clark_exome_BAM_AL_DIR=Clark_exome_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES),
		expand(Clark_exome_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam", Clark_exome_SORTED_BAM_AL_DIR=Clark_exome_SORTED_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES),
		expand(Clark_exome_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam.bai", Clark_exome_SORTED_BAM_AL_DIR=Clark_exome_SORTED_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES)


rule hisat2_sayres_align_reads:
	input:
		R1 = Sayres_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		R2 = Sayres_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed_paired_2.fastq",
	output:
		SAM = Sayres_RNA_SAM_AL_DIR + "{sample}.sam"
	params:
		hisat2_index = VITRIPENNIS_HISAT2_INDEX,
		threads = 8
	message: "Mapping reads to {params.hisat2_index} with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.hisat2_index} -s no -1 {input.R1} -2 {input.R2} -S {output.SAM}
		"""

rule hisat2_clark_rna_align_reads:
	input:
		R1 = Clark_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed.fastq",
	output:
		SAM = Clark_RNA_SAM_AL_DIR + "{sample}.sam"
	params:
		hisat2_index = VITRIPENNIS_HISAT2_INDEX,
		threads = 8
	message: "Mapping reads to {params.hisat2_index} with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.hisat2_index} -s no -U {input.R1} -S {output.SAM}
		"""

rule hisat2_clark_exome_align_reads:
	input:
		R1 = Clark_exome_trimmed_FASTQ + "{sample}_trimmomatic_trimmed.fastq",
	output:
		SAM = Clark_exome_SAM_AL_DIR + "{sample}.sam"
	params:
		hisat2_index = VITRIPENNIS_HISAT2_INDEX,
		threads = 8
	message: "Mapping reads to {params.hisat2_index} with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.hisat2_index} -s no -U {input.R1} -S {output.SAM}
		"""

rule sayres_sam_to_bam:
	input:
		SAM = Sayres_RNA_SAM_AL_DIR + "{sample}.sam"
	output:
		BAM = Sayres_RNA_BAM_AL_DIR + "{sample}.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		samtools view -b -F 4 {input.SAM} > {output.BAM}
		"""

rule clark_rna_sam_to_bam:
	input:
		SAM = Clark_RNA_SAM_AL_DIR + "{sample}.sam"
	output:
		BAM = Clark_RNA_BAM_AL_DIR + "{sample}.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		samtools view -b -F 4 {input.SAM} > {output.BAM}
		"""


rule clark_exome_sam_to_bam:
	input:
		SAM = Clark_exome_SAM_AL_DIR + "{sample}.sam"
	output:
		BAM = Clark_exome_BAM_AL_DIR + "{sample}.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		samtools view -b -F 4 {input.SAM} > {output.BAM}
		"""


rule sayres_sort_bam:
	input:
		BAM = Sayres_RNA_BAM_AL_DIR + "{sample}.bam"
	output:
		SORTED_BAM = Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		samtools sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule clark_rna_sort_bam:
	input:
		BAM = Clark_RNA_BAM_AL_DIR + "{sample}.bam"
	output:
		SORTED_BAM = Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		samtools sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule clark_exome_sort_bam:
	input:
		BAM = Clark_exome_BAM_AL_DIR + "{sample}.bam"
	output:
		SORTED_BAM = Clark_exome_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		samtools sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule sayres_index_bam:
	input:
		BAM = Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam"
	output: Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		samtools index {input.BAM}
		"""

rule clark_rna_index_bam:
	input:
		BAM = Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam"
	output: Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		samtools index {input.BAM}
		"""

rule clark_exome_index_bam:
	input:
		BAM = Clark_exome_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam"
	output: Clark_exome_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		samtools index {input.BAM}
		"""

rule sayres_featurecounts:
	input:
		BAM = Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam",
		GFF = VITRIPENNIS_GFF
	output:
		COUNTS = Sayres_RNA_FEATURECOUNTS_DIR + "{sample}_featurecounts.tsv"
	params:
		THREADS = 5
	message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
	shell:
		"""
		featureCounts -T {params.THREADS} -O --primary -p -s 0 -t gene -g ID -F GFF -a {input.GFF} -o {output.COUNTS} {input.BAM}
		"""

rule clark_featurecounts:
	input:
		BAM = Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord.bam",
		GFF = VITRIPENNIS_GFF
	output:
		COUNTS = Clark_RNA_FEATURECOUNTS_DIR + "{sample}_featurecounts.tsv"
	params:
		THREADS = 5
	message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
	shell:
		"""
		featureCounts -T {params.THREADS} -O --primary -p -s 0 -t gene -g ID -F GFF -a {input.GFF} -o {output.COUNTS} {input.BAM}
		"""
