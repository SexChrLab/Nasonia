# Workflow for calling germline SNP and indel variants using GATK HaplotypeCaller and ASE quantification using GATK ASEReadCounter.

configfile: "Nasonia.config.json"

# Reference files: reference genome sequences in FASTA format abd an index .fa.fai
# files created with Samtools faidx.
REF = config["Vitripennis_ref_index_path"]
REF_FA_INDEX = config["Vitripennis_ref_index_path"] # Reference FASTA index file (.fa.fai)

# Directories
Sayres_RNA_SORTED_BAM_AL_DIR = config["Sayres_RNA_processed_BAM"]
Clark_RNA_SORTED_BAM_AL_DIR = config["Clark_RNA_processed_BAM"]
Clark_exome_SORTED_BAM_AL_DIR = config["Clark_exome_processed_BAM"]

Sayres_RNA_gVCF = config["Sayres_RNA_gVCF"]
Clark_RNA_gVCF = config["Clark_RNA_gVCF"]
Clark_exome_gVCF = config["Clark_exome_gVCF"]

VCF_DIR = "/data/storage/SAYRES/NASONIA/Heini/VCF/"
SAYRES_ASE_DIR = config["Sayres_RNA_ASE"]
CLARK_ASE_DIR = config["Clark_RNA_ASE"]

# Samples
SAYRES_SAMPLES = config["SAYRES_SAMPLES"]
CLARK_RNA_SAMPLES= config["CLARK_RNA_SAMPLES"]
CLARK_EXOME_SAMPLES = config["CLARK_EXOME_SAMPLES"]

# Formatting strings for command line
GATK_GVCF_LIST = []
ALL_GVCFS = config["ALL_gir_gVCF_paths"]
for item in ALL_GVCFS:
	gvcf_list_item = "--variant " + item
	GATK_GVCF_LIST.append(gvcf_list_item)

GATK_GVCF_LIST_str = " ".join(GATK_GVCF_LIST)

# Output files
rule all:
	input:
		expand(Sayres_RNA_gVCF + "{sample}_rawLikehoods.g.vcf", Sayres_RNA_gVCF=Sayres_RNA_gVCF, sample=SAYRES_SAMPLES),
		expand(Clark_RNA_gVCF + "{sample}_rawLikehoods.g.vcf", Clark_RNA_gVCF=Clark_RNA_gVCF, sample=CLARK_RNA_SAMPLES),
		expand(Clark_exome_gVCF + "{sample}_rawLikehoods.g.vcf", Clark_exome_gVCF=Clark_exome_gVCF, sample=CLARK_EXOME_SAMPLES),
		VCF_DIR + "SAYRES_CLARK_all_gVCFs.g.vcf",
		VCF_DIR + "SAYRES_CLARK_all_genotypes.vcf",
		expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth100.csv", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
		expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth100.csv", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),

rule sayres_gatk_haplotypecaller_gvcf:
	input:
		BAM = lambda wildcards: Sayres_RNA_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
		REF = REF
	output:
		VCF = Sayres_RNA_gVCF + "{sample}_rawLikehoods.g.vcf"
	params:
		xmx = "48g",
		xms = "48g",
		threads = "16",
		min_pruning = 2,
		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
		maxReadsInRegionPerSample=10000, # Maximum reads in an active region
		min_base_quality_score=10, # Minimum base quality required to consider a base for calling
		minReadsPerAlignmentStart=10, # Minimum number of reads sharing the same alignment start for each genomic location in an active region
		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
		output_mode="EMIT_ALL_CONFIDENT_SITES" # Which type of calls we should output: EMIT_VARIANTS_ONLY, EMIT_ALL_SITES
	message: "Calculating genotype likelihoods from {input.BAM} using GATK HaplotypeCaller"
	shell:
		"""
		gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
		HaplotypeCaller \
		--native-pair-hmm-threads {params.threads} \
		--emit-ref-confidence GVCF \
		--reference {input.REF} \
		--input {input.BAM} \
		--min-pruning {params.min_pruning} \
		--heterozygosity {params.heterozygosity} \
		--indel-heterozygosity {params.indel_heterozygosity} \
		--max-reads-per-alignment-start {params.maxReadsInRegionPerSample} \
		--min-base-quality-score {params.min_base_quality_score} \
		--sample-ploidy {params.sample_ploidy} \
		--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
		--max-alternate-alleles {params.max_alternate_alleles} \
		--output-mode {params.output_mode} \
		--output {output.VCF}
		"""

rule clark_exome_gatk_haplotypecaller_gvcf:
	input:
		BAM = lambda wildcards: Clark_exome_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
		REF = REF
	output:
		VCF = Clark_exome_gVCF + "{sample}_rawLikehoods.g.vcf"
	params:
		xmx = "48g",
		xms = "48g",
		threads = "16",
		min_pruning = 2,
		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
		maxReadsInRegionPerSample=10000, # Maximum reads in an active region
		min_base_quality_score=10, # Minimum base quality required to consider a base for calling
		minReadsPerAlignmentStart=10, # Minimum number of reads sharing the same alignment start for each genomic location in an active region
		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
		output_mode="EMIT_ALL_CONFIDENT_SITES" # Which type of calls we should output: EMIT_VARIANTS_ONLY, EMIT_ALL_SITES
	message: "Calculating genotype likelihoods from {input.BAM} using GATK HaplotypeCaller"
	shell:
		"""
		gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
		HaplotypeCaller \
		--native-pair-hmm-threads {params.threads} \
		--emit-ref-confidence GVCF \
		--reference {input.REF} \
		--input {input.BAM} \
		--min-pruning {params.min_pruning} \
		--heterozygosity {params.heterozygosity} \
		--indel-heterozygosity {params.indel_heterozygosity} \
		--max-reads-per-alignment-start {params.maxReadsInRegionPerSample} \
		--min-base-quality-score {params.min_base_quality_score} \
		--sample-ploidy {params.sample_ploidy} \
		--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
		--max-alternate-alleles {params.max_alternate_alleles} \
		--output-mode {params.output_mode} \
		--output {output.VCF}
		"""

rule clark_rna_gatk_haplotypecaller_gvcf:
	input:
		BAM = lambda wildcards: Clark_RNA_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
		REF = REF
	output:
		VCF = Clark_RNA_gVCF + "{sample}_rawLikehoods.g.vcf"
	params:
		xmx = "48g",
		xms = "48g",
		threads = "16",
		min_pruning = 2,
		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
		maxReadsInRegionPerSample=10000, # Maximum reads in an active region
		min_base_quality_score=10, # Minimum base quality required to consider a base for calling
		minReadsPerAlignmentStart=10, # Minimum number of reads sharing the same alignment start for each genomic location in an active region
		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
		output_mode="EMIT_ALL_CONFIDENT_SITES" # Which type of calls we should output: EMIT_VARIANTS_ONLY, EMIT_ALL_SITES
	message: "Calculating genotype likelihoods from {input.BAM} using GATK HaplotypeCaller"
	shell:
		"""
		gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
		HaplotypeCaller \
		--native-pair-hmm-threads {params.threads} \
		--emit-ref-confidence GVCF \
		--reference {input.REF} \
		--input {input.BAM} \
		--min-pruning {params.min_pruning} \
		--heterozygosity {params.heterozygosity} \
		--indel-heterozygosity {params.indel_heterozygosity} \
		--max-reads-per-alignment-start {params.maxReadsInRegionPerSample} \
		--min-base-quality-score {params.min_base_quality_score} \
		--sample-ploidy {params.sample_ploidy} \
		--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
		--max-alternate-alleles {params.max_alternate_alleles} \
		--output-mode {params.output_mode} \
		--output {output.VCF}
		"""

rule all_samples_merge_gvcfs:
	input:
		REF = REF
	output:
		VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_gVCFs.g.vcf")
	params:
		GVCF_LIST = GATK_GVCF_LIST_str,
		xmx = "62g",
		xms = "62g",
		threads = "24"
	message: "Merging gVCF records."
	shell:
		"""
		java -Xmx16g -jar /mnt/storage/SAYRES/NASONIA/01_tools/GenomeAnalysisTK.jar \
		-T CombineGVCFs \
		-R {input.REF} \
		{params.GVCF_LIST} \
		-o {output.VCF}
		"""

rule all_samples_gatk_genotype_gvcfs:
	input:
		REF = REF,
		VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_gVCFs.g.vcf")
	output:
		VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_genotypes.vcf")
	params:
		#GVCF_LIST = GATK_GVCF_LIST_str,
		xmx = "62g",
		xms = "62g",
		threads = "24",
		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
	message: "Merging gVCF records."
	shell:
		"""
		gatk-launch --java-options "-Xmx{params.xmx} -Xms{params.xms}" \
		GenotypeGVCFs \
		--reference {input.REF} \
		--variant {input.VCF} \
		--heterozygosity {params.heterozygosity} \
		--indel-heterozygosity {params.indel_heterozygosity} \
		--sample-ploidy {params.sample_ploidy} \
		--standard-min-confidence-threshold-for-calling {params.standard_min_confidence_threshold_for_calling} \
		--max-alternate-alleles {params.max_alternate_alleles} \
		--output {output.VCF}
		"""


rule sayres_asereadcounter_gvcfs:
	input:
		REF = REF,
		BAM = lambda wildcards: Sayres_RNA_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
		VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_genotypes.vcf")
	output:
		ASE = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth100.csv")
	params:
		xmx = "62g",
		xms = "62g",
		threads = "24",
		minDepth=100, 
		minMappingQuality=10,
		minBaseQuality=2
	message: "Running ASEReadCounter."
	shell:
		"""
		java -Xmx16g -jar /mnt/storage/SAYRES/NASONIA/01_tools/GenomeAnalysisTK.jar \
		-T ASEReadCounter \
		-R {input.REF} \
		-I {input.BAM} \
		--sitesVCFFile {input.VCF} \
		-U ALLOW_N_CIGAR_READS \
		-minDepth {params.minDepth} \
		--minMappingQuality {params.minMappingQuality} \
		--minBaseQuality {params.minBaseQuality} \
		-drf DuplicateRead \
		-o {output.ASE}
		"""

rule clark_asereadcounter_gvcfs:
	input:
		REF = REF,
		BAM = lambda wildcards: Clark_RNA_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
		VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_genotypes.vcf")
	output:
		ASE = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth100.csv")
	params:
		xmx = "62g",
		xms = "62g",
		threads = "24",
		minDepth=100, 
		minMappingQuality=10,
		minBaseQuality=2
	message: "Running ASEReadCounter."
	shell:
		"""
		java -Xmx16g -jar /mnt/storage/SAYRES/NASONIA/01_tools/GenomeAnalysisTK.jar \
		-T ASEReadCounter \
		-R {input.REF} \
		-I {input.BAM} \
		--sitesVCFFile {input.VCF} \
		-U ALLOW_N_CIGAR_READS \
		-minDepth {params.minDepth} \
		--minMappingQuality {params.minMappingQuality} \
		--minBaseQuality {params.minBaseQuality} \
		-drf DuplicateRead \
		-o {output.ASE}
		"""
