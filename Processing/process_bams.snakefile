# Workflow for preprocessing BAM files for variant calling
from os.path import join

configfile: "Nasonia.config.json"

# Directories
Sayres_RNA_trimmed_FASTQ = config["Sayres_RNA_trimmed_FASTQ"]
Clark_RNA_trimmed_FASTQ = config["Clark_RNA_trimmed_FASTQ"]
Clark_exome_trimmed_FASTQ = config["Clark_exome_trimmed_FASTQ"]

Sayres_RNA_SAM_AL_DIR = config["Sayres_RNA_SAM"]
Sayres_RNA_BAM_AL_DIR = config["Sayres_RNA_BAM"]
Sayres_RNA_SORTED_BAM_AL_DIR = config["Sayres_RNA_processed_BAM"]
Sayres_RNA_FEATURECOUNTS_DIR = config["Sayres_RNA_featureCounts"]

Clark_RNA_SAM_AL_DIR = config["Clark_RNA_processed_BAM"]
Clark_RNA_BAM_AL_DIR = config["Clark_RNA_processed_BAM"]
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
		expand(Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs.bam", Sayres_RNA_SORTED_BAM_AL_DIR=Sayres_RNA_SORTED_BAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs.bam.bai", Sayres_RNA_SORTED_BAM_AL_DIR=Sayres_RNA_SORTED_BAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs_mkDups.bam", Sayres_RNA_SORTED_BAM_AL_DIR=Sayres_RNA_SORTED_BAM_AL_DIR, sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs_mkDups.bam.bai", Sayres_RNA_SORTED_BAM_AL_DIR=Sayres_RNA_SORTED_BAM_AL_DIR, sample=SAYRES_SAMPLES),

		expand(Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs.bam", Clark_RNA_SORTED_BAM_AL_DIR=Clark_RNA_SORTED_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs.bam.bai", Clark_RNA_SORTED_BAM_AL_DIR=Clark_RNA_SORTED_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs_mkDups.bam", Clark_RNA_SORTED_BAM_AL_DIR=Clark_RNA_SORTED_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),
		expand(Clark_RNA_SORTED_BAM_AL_DIR + "{sample}_sortedbycoord_wRGs_mkDups.bam.bai", Clark_RNA_SORTED_BAM_AL_DIR=Clark_RNA_SORTED_BAM_AL_DIR, sample=CLARK_RNA_SAMPLES),

		expand(Clark_exome_SORTED_BAM_AL_DIR+ "{sample}_sortedbycoord_wRGs.bam", Clark_exome_SORTED_BAM_AL_DIR=Clark_exome_SORTED_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES),
		expand(Clark_exome_SORTED_BAM_AL_DIR+ "{sample}_sortedbycoord_wRGs.bam.bai", Clark_exome_SORTED_BAM_AL_DIR=Clark_exome_SORTED_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES),
		expand(Clark_exome_SORTED_BAM_AL_DIR+ "{sample}_sortedbycoord_wRGs_mkDups.bam", Clark_exome_SORTED_BAM_AL_DIR=Clark_exome_SORTED_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES),
		expand(Clark_exome_SORTED_BAM_AL_DIR+ "{sample}_sortedbycoord_wRGs_mkDups.bam.bai", Clark_exome_SORTED_BAM_AL_DIR=Clark_exome_SORTED_BAM_AL_DIR, sample=CLARK_EXOME_SAMPLES)


rule sayres_picard_add_readgroups:
	input:
		BAM = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord.bam"),
	output:
		BAM = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam"),
	params:
		SAMPLENAME = "{sample}",
	message: "Adding readgroups to {input.BAM}."
	shell:
		"picard AddOrReplaceReadGroups I={input.BAM} O={output.BAM} RGID=0 RGLB=lib1 RGPL=Illumina RGPU=unit1 RGSM={params.SAMPLENAME}"

rule clark_rna_picard_add_readgroups:
	input:
		BAM = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord.bam"),
	output:
		BAM = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam"),
	params:
		SAMPLENAME = "{sample}",
	message: "Adding readgroups to {input.BAM}."
	shell:
		"picard AddOrReplaceReadGroups I={input.BAM} O={output.BAM} RGID=0 RGLB=lib1 RGPL=Illumina RGPU=unit1 RGSM={params.SAMPLENAME}"

rule clark_exome_picard_add_readgroups:
	input:
		BAM = os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord.bam"),
	output:
		BAM = os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam"),
	params:
		SAMPLENAME = "{sample}",
	message: "Adding readgroups to {input.BAM}."
	shell:
		"picard AddOrReplaceReadGroups I={input.BAM} O={output.BAM} RGID=0 RGLB=lib1 RGPL=Illumina RGPU=unit1 RGSM={params.SAMPLENAME}"

rule sayres_index_rg_bam:
	input:
		BAM = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam"),
	output:
		os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam.bai"), # Missing from rule all
	message:
		"Indexing BAM file {input.BAM} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")

rule clark_rna_index_rg_bam:
	input:
		BAM = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam"),
	output:
		os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam.bai"), # Missing from rule all
	message:
		"Indexing BAM file {input.BAM} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")

rule clark_exome_index_rg_bam:
	input:
		BAM = os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam"),
	output:
		os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam.bai"), # Missing from rule all
	message:
		"Indexing BAM file {input.BAM} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")


rule sayres_picard_mark_duplicates:
	input:
		BAM = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam")
	output:
		BAM = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam"),
		METRICS = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR,  "{sample}_marked_dup_metrics.txt") # Missing from rule all
	params:
		SAMPLENAME = "{sample}",
	message: "Removing duplicates from {input.BAM}."
	shell:
		"picard MarkDuplicates I={input.BAM} O={output.BAM} M={output.METRICS} REMOVE_SEQUENCING_DUPLICATES=true"

rule clark_rna_picard_mark_duplicates:
	input:
		BAM = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam")
	output:
		BAM = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam"),
		METRICS = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR,  "{sample}_marked_dup_metrics.txt") # Missing from rule all
	params:
		SAMPLENAME = "{sample}",
	message: "Removing duplicates from {input.BAM}."
	shell:
		"picard MarkDuplicates I={input.BAM} O={output.BAM} M={output.METRICS} REMOVE_SEQUENCING_DUPLICATES=true"

rule clark_exome__picard_mark_duplicates:
	input:
		BAM = os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs.bam")
	output:
		BAM = os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam"),
		METRICS = os.path.join(Clark_exome_SORTED_BAM_AL_DIR,  "{sample}_marked_dup_metrics.txt") # Missing from rule all
	params:
		SAMPLENAME = "{sample}",
	message: "Removing duplicates from {input.BAM}."
	shell:
		"picard MarkDuplicates I={input.BAM} O={output.BAM} M={output.METRICS} REMOVE_SEQUENCING_DUPLICATES=true"


rule sayres_index_picard_bam:
	input:
		BAM = os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam")
	output:
		os.path.join(Sayres_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam.bai")
	message:
		"Indexing BAM file {input.BAM} and with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")

rule clark_rna_index_picard_bam:
	input:
		BAM = os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam")
	output:
		os.path.join(Clark_RNA_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam.bai")
	message:
		"Indexing BAM file {input.BAM} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")

rule clark_exome_index_picard_bam:
	input:
		BAM = os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam")
	output:
		os.path.join(Clark_exome_SORTED_BAM_AL_DIR, "{sample}_sortedbycoord_wRGs_mkDups.bam.bai")
	message:
		"Indexing BAM file {input.BAM} with Samtools."
	params:
	run:
		for x in input:
			shell("samtools index {x}")
