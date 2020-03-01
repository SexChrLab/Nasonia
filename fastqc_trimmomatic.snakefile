# Workflow for FastQC, MultiQC, and adapter trimming using Trimmomatic.
from os.path import join

configfile: "Nasonia.config.json"

# Directory variables
Sayres_RNA_raw_FASTQ = config["Sayres_RNA_raw_FASTQ"]
Sayres_RNA_raw_FASTQ_FastQC = config["Sayres_RNA_raw_FASTQ_FastQC"]
Sayres_RNA_trimmed_FASTQ = config["Sayres_RNA_trimmed_FASTQ"]
Sayres_RNA_trimmed_FASTQ_FastQC = config["Sayres_RNA_trimmed_FASTQ_FastQC"]
Clark_RNA_raw_FASTQ = config["Clark_RNA_raw_FASTQ"]
Clark_RNA_raw_FASTQ_FastQC = config["Clark_RNA_raw_FASTQ_FastQC"]
Clark_RNA_trimmed_FASTQ = config["Clark_RNA_trimmed_FASTQ"]
Clark_RNA_trimmed_FASTQ_FastQC = config["Clark_RNA_trimmed_FASTQ_FastQC"]
Clark_exome_raw_FASTQ = config["Clark_exome_raw_FASTQ"]
Clark_exome_raw_FASTQ_FastQC = config["Clark_exome_raw_FASTQ_FastQC"]
Clark_exome_trimmed_FASTQ = config["Clark_exome_trimmed_FASTQ"]
Clark_exome_trimmed_FASTQ_FastQC = config["Clark_exome_trimmed_FASTQ_FastQC"]

# Tools
fastqc_path = "fastqc"
multiqc_path = "multiqc"
trimmomatic_path = "trimmomatic"

# Samples
SAYRES_SAMPLES = config["SAYRES_SAMPLES"]
CLARK_RNA_SAMPLES= config["CLARK_RNA_SAMPLES"]
CLARK_EXOME_SAMPLES = config["CLARK_EXOME_SAMPLES"]

#ruleorder: fastqc_analysis > multiqc > trimmomatic > fastqc_analysis_trimmomatic_trimmed_paired > fastqc_analysis_trimmomatic_trimmed_paired > multiqc_trimmed_paired

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		expand(Sayres_RNA_raw_FASTQ_FastQC + "{sample}_fq1_fastqc.html", sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_raw_FASTQ_FastQC + "{sample}_fq2_fastqc.html", sample=SAYRES_SAMPLES),
		expand(Clark_RNA_raw_FASTQ_FastQC + "{sample}_fastqc.html", sample=CLARK_RNA_SAMPLES),
		expand(Clark_exome_raw_FASTQ_FastQC + "{sample}_fastqc.html", sample=CLARK_EXOME_SAMPLES),

		expand(Sayres_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed_paired_1.fastq", sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed_unpaired_1.fastq", sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed_paired_2.fastq", sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed_unpaired_2.fastq", sample=SAYRES_SAMPLES),

		expand(Sayres_RNA_trimmed_FASTQ_FastQC + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html", sample=SAYRES_SAMPLES),
		expand(Sayres_RNA_trimmed_FASTQ_FastQC + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html", sample=SAYRES_SAMPLES),

		expand(Clark_RNA_trimmed_FASTQ + "{sample}_trimmomatic_trimmed.fastq", sample=CLARK_RNA_SAMPLES),
		expand(Clark_exome_trimmed_FASTQ + "{sample}_trimmomatic_trimmed.fastq", sample=CLARK_EXOME_SAMPLES),

		expand(Clark_RNA_trimmed_FASTQ_FastQC + "{sample}_trimmomatic_trimmed_fastqc.html", sample=CLARK_RNA_SAMPLES),
		expand(Clark_exome_trimmed_FASTQ_FastQC + "{sample}_trimmomatic_trimmed_fastqc.html", sample=CLARK_EXOME_SAMPLES),

rule sayres_fastqc_analysis:
	input:
		fq1 = lambda wildcards: Sayres_RNA_raw_FASTQ + config["SAYRES_RAW_FASTQS"][wildcards.sample][0] + ".fastq",
		fq2 = lambda wildcards: Sayres_RNA_raw_FASTQ + config["SAYRES_RAW_FASTQS"][wildcards.sample][1] + ".fastq"
	output:
		fq1_zip =  os.path.join(Sayres_RNA_raw_FASTQ_FastQC, "{sample}_fq1_fastqc.zip"),
		fq1_html = os.path.join(Sayres_RNA_raw_FASTQ_FastQC, "{sample}_fq1_fastqc.html"),
		fq2_zip =  os.path.join(Sayres_RNA_raw_FASTQ_FastQC, "{sample}_fq2_fastqc.zip"),
		fq2_html = os.path.join(Sayres_RNA_raw_FASTQ_FastQC, "{sample}_fq2_fastqc.html")
	params:
		fastqc = fastqc_path,
		fastqc_dir = Sayres_RNA_raw_FASTQ_FastQC,
		fq1_prefix = lambda wildcards: config["SAYRES_RAW_FASTQS"][wildcards.sample][0],
		fq2_prefix = lambda wildcards: config["SAYRES_RAW_FASTQS"][wildcards.sample][1],
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1};
		{params.fastqc} -o {params.fastqc_dir} {input.fq2};
		mv {params.fastqc_dir}{params.fq1_prefix}_fastqc.html {output.fq1_html};
		mv {params.fastqc_dir}{params.fq1_prefix}_fastqc.zip {output.fq1_zip};
		mv {params.fastqc_dir}{params.fq2_prefix}_fastqc.html {output.fq2_html};
		mv {params.fastqc_dir}{params.fq2_prefix}_fastqc.zip {output.fq2_zip}
		"""

rule clark_rna_fastqc_analysis:
	input:
		fq1 = Clark_RNA_raw_FASTQ + "{sample}.fastq",
	output:
		fq1_zip =  os.path.join(Clark_RNA_raw_FASTQ_FastQC, "{sample}_fastqc.zip"),
		fq1_html = os.path.join(Clark_RNA_raw_FASTQ_FastQC, "{sample}_fastqc.html"),
	params:
		fastqc = fastqc_path,
		fastqc_dir = Clark_RNA_raw_FASTQ_FastQC,
		fq1_prefix = "{sample}",
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1}
		"""

rule clark_exome_fastqc_analysis:
	input:
		fq1 = Clark_exome_raw_FASTQ + "{sample}.fastq"
	output:
		fq1_zip =  os.path.join(Clark_exome_raw_FASTQ_FastQC, "{sample}_fastqc.zip"),
		fq1_html = os.path.join(Clark_exome_raw_FASTQ_FastQC, "{sample}_fastqc.html")
	params:
		fastqc = fastqc_path,
		fastqc_dir = Clark_exome_raw_FASTQ_FastQC,
		fq1_prefix = "{sample}",
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1}
		"""

rule sayres_trimmomatic:
	input:
		fq1 = lambda wildcards: Sayres_RNA_raw_FASTQ + config["SAYRES_RAW_FASTQS"][wildcards.sample][0],
		fq2 = lambda wildcards: Sayres_RNA_raw_FASTQ + config["SAYRES_RAW_FASTQS"][wildcards.sample][1],
		ADAPTER_FASTA = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/adapter_sequences.fa"
	output:
		paired_1 = os.path.join(Sayres_RNA_trimmed_FASTQ,"{sample}_trimmomatic_trimmed_paired_1.fastq"),
		unpaired_1 = os.path.join(Sayres_RNA_trimmed_FASTQ,"{sample}_trimmomatic_trimmed_unpaired_1.fastq"),
		paired_2 = os.path.join(Sayres_RNA_trimmed_FASTQ,"{sample}_trimmomatic_trimmed_paired_2.fastq"),
		unpaired_2 = os.path.join(Sayres_RNA_trimmed_FASTQ,"{sample}_trimmomatic_trimmed_unpaired_2.fastq"),
		logfile = os.path.join(Sayres_RNA_trimmed_FASTQ,"{sample}_trimmomatic.log")
	params:
		trimmomatic = trimmomatic_path,
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 10,
		trailing = 10,
		winsize = 4,
		winqual = 15,
		minlen = 80
	shell:
		"""
		{params.trimmomatic} PE -threads {params.threads} -phred33 -trimlog {output.logfile} \
		{input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} \
		{output.paired_2} {output.unpaired_2} \
		ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
		LEADING:{params.leading} TRAILING:{params.trailing} \
		SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}
		"""

rule clark_rna_trimmomatic:
	input:
		fq1 = os.path.join(Clark_RNA_raw_FASTQ, "{sample}.fastq"),
		ADAPTER_FASTA = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/adapter_sequences.fa"
	output:
		trimmed_1 =   os.path.join(Clark_RNA_trimmed_FASTQ,"{sample}_trimmomatic_trimmed.fastq"),
		logfile = os.path.join(Clark_RNA_trimmed_FASTQ,"{sample}_trimmomatic.log")
	params:
		trimmomatic = trimmomatic_path,
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 8,
		trailing = 8,
		winsize = 4,
		winqual = 15,
		minlen = 40
	shell:
		"""
		{params.trimmomatic} SE -threads {params.threads} -phred33 -trimlog {output.logfile} \
		{input.fq1} {output.trimmed_1} \
		ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
		LEADING:{params.leading} TRAILING:{params.trailing} \
		SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}
		"""

rule clark_exome_trimmomatic:
	input:
		fq1 = os.path.join(Clark_exome_raw_FASTQ, "{sample}.fastq"),
		ADAPTER_FASTA = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/adapter_sequences.fa"
	output:
		trimmed_1 = os.path.join(Clark_exome_trimmed_FASTQ,"{sample}_trimmomatic_trimmed.fastq"),
		logfile = os.path.join(Clark_exome_trimmed_FASTQ,"{sample}_trimmomatic.log")
	params:
		trimmomatic = trimmomatic_path,
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 8,
		trailing = 8,
		winsize = 4,
		winqual = 15,
		minlen = 40
	shell:
		"""
		{params.trimmomatic} SE -threads {params.threads} -phred33 -trimlog {output.logfile} \
		{input.fq1} {output.trimmed_1} \
		ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
		LEADING:{params.leading} TRAILING:{params.trailing} \
		SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}
		"""

rule sayres_fastqc_analysis_trimmomatic_trimmed_paired:
	input:
		fq1 = os.path.join(Sayres_RNA_trimmed_FASTQ , "{sample}_trimmomatic_trimmed_paired_1.fastq"),
		fq2 = os.path.join(Sayres_RNA_trimmed_FASTQ , "{sample}_trimmomatic_trimmed_paired_2.fastq")
	output:
		html1 = os.path.join(Sayres_RNA_trimmed_FASTQ_FastQC, "{sample}_trimmomatic_trimmed_paired_1_fastqc.html"),
		html2 = os.path.join(Sayres_RNA_trimmed_FASTQ_FastQC, "{sample}_trimmomatic_trimmed_paired_2_fastqc.html")
	params:
		fastqc = fastqc_path,
		fastqc_dir = Sayres_RNA_trimmed_FASTQ_FastQC
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1} {input.fq2}
		"""

rule clark_rna_fastqc_analysis_trimmomatic_trimmed_se:
	input:
		fq1 = os.path.join(Clark_RNA_trimmed_FASTQ , "{sample}_trimmomatic_trimmed.fastq"),
	output:
		html1 = os.path.join(Clark_RNA_trimmed_FASTQ_FastQC, "{sample}_trimmomatic_trimmed_fastqc.html")
	params:
		fastqc = fastqc_path,
		fastqc_dir = Clark_RNA_trimmed_FASTQ_FastQC
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1}
		"""

rule clark_exome_fastqc_analysis_trimmomatic_trimmed_se:
	input:
		fq1 = os.path.join(Clark_exome_trimmed_FASTQ , "{sample}_trimmomatic_trimmed.fastq"),
	output:
		html1 = os.path.join(Clark_exome_trimmed_FASTQ_FastQC, "{sample}_trimmomatic_trimmed_fastqc.html"),
	params:
		fastqc = fastqc_path,
		fastqc_dir = Clark_exome_trimmed_FASTQ_FastQC
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1}
		"""

