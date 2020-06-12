# 1) Filter ASEReadCounter .csv for only heterozygous sites
# 2) Filter .csv to only include sites that are fixed and different between VV and GG <- FD sites are determined in before hand
# 3) Intersect .csv with gene annotation file .gff using bedtools to determine geneID with snp information
# 4) Determine if all SNPs within a gene show bias for the same haplotype or not (i.e VV or GG haplotype)

# Heini M Natri, hnatri@asu.edu
# Kimberly Olney, kcolney@asu.edu

configfile: "Nasonia.config.json"

# Reference files: reference genome sequences in FASTA format abd an index .fa.fai
# files created with Samtools faidx.
REF = config["Vitripennis_ref_path"]
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

FD_DIR = "/data/storage/SAYRES/NASONIA/Heini/FixedDifferences/"

# Samples
SAYRES_SAMPLES = config["SAYRES_SAMPLES"]
CLARK_RNA_SAMPLES= config["CLARK_RNA_SAMPLES"]
CLARK_EXOME_SAMPLES = config["CLARK_EXOME_SAMPLES"]

# Formatting strings for command line
GATK_GVCF_LIST = []
ALL_GVCFS = config["ALL_gVCF_paths"]
for item in ALL_GVCFS:
    gvcf_list_item = "--variant " + item
    GATK_GVCF_LIST.append(gvcf_list_item)

GATK_GVCF_LIST_str = " ".join(GATK_GVCF_LIST)

# Output files
rule all:
    input:
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30.csv", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30.csv", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET.csv", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET.csv", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD.csv", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD.csv", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_intersect.csv", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_intersect.csv", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD.bed", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD.bed", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),

        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF.bed", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_sed.bed", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_cut.bed", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_grep.bed", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),
        expand(SAYRES_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_ratio.bed", SAYRES_ASE_DIR=SAYRES_ASE_DIR, sample=SAYRES_SAMPLES),

        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF.bed", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_sed.bed", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_cut.bed", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_grep.bed", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),
        expand(CLARK_ASE_DIR + "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_ratio.bed", CLARK_ASE_DIR=CLARK_ASE_DIR, sample=CLARK_RNA_SAMPLES),


# Sayres 
rule sayres_asereadcounter_gvcfs:
    input:
        REF = REF,
        BAM = lambda wildcards: Sayres_RNA_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
        VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_genotypes.vcf")
    output:
        ASE = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30.csv")
    params:
        xmx = "62g",
        xms = "62g",
        threads = "24",
        minDepth=30, 
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

rule sayres_asereadcounter_gvcfs_HET:
    input:
        ASE = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30.csv")
    output:
        ASE_HET = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET.csv")
    params:
        xmx = "62g",
        xms = "62g",
        threads = "24",
        minDepth=30, 
        minMappingQuality=10,
        minBaseQuality=2
    message: "Filtering ASEReadCounter for HET sites."
    shell:
        """
        python3 getOnlyHeterozygousSites.py {input.ASE} {output.ASE_HET}
        """

rule sayres_asereadcounter_gvcfs_HET_FD:
    input:
        Sayres_FD = os.path.join(FD_DIR, "Sayres_FD_NoHeader_sites.txt"),
        ASE_HET = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET.csv")        
    output:
        ASE_HET_FD = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.csv")
    message: "Filtering ASEReadCounter for FD sites."
    shell:
        """
        awk -F'\t' 'NR==FNR{{c[$1,$2]++;next}};c[$1,$2]>0' {input.Sayres_FD} {input.ASE_HET} > {output.ASE_HET_FD}
        """

rule sayres_asereadcounter_gvcfs_HET_FD_intersect:
    input:
        Intersect_FD = os.path.join(FD_DIR, "Sayres_Clark_intersect.txt"),
        ASE_HET = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET.csv")        
    output:
        ASE_HET_FD_inter = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_intersect.csv")
    message: "Filtering ASEReadCounter for FD sites."
    shell:
        """
        awk -F'\t' 'NR==FNR{{c[$1,$2]++;next}};c[$1,$2]>0' {input.Intersect_FD} {input.ASE_HET} > {output.ASE_HET_FD_inter}
        """
        
rule sayres_ase_to_bed:
    input:
        ASE_HET_FD = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.csv")    
    output:
        ASE_to_BED = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.bed")
    message: "convert ASE to BED."
    shell:
        """
        python3 convert_ase_to_bed_fmt.py --ase_filename {input.ASE_HET_FD} --output_filename {output.ASE_to_BED}
        """
        
rule sayres_ase_to_bed_gff_intersect:
    input:
        ASE_HET_FD_BED = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.bed")
    output:
        ASE_to_BED_GFF = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF.bed")
    message: "ASE BED with GFF gene IDs."
    shell:
        """
        bedtools intersect  -wa -wb -a {input.ASE_HET_FD_BED} -b Reference/GCF_009193385.2_Nvit_psr_1.1_genomic.gff > {output.ASE_to_BED_GFF}
        """

rule sayres_ase_gff_sed:
    input:
        ASE_to_BED_GFF = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF.bed")
    output:
        ASE_GFF_sed = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_sed.bed")
    message: "ASE BED with GFF gene IDs, sed;."
    shell:
        """
        sed 's/;/\t/g' {input.ASE_to_BED_GFF} > {output.ASE_GFF_sed}
        """

rule sayres_ase_gff_cut:
    input:
        ASE_GFF_sed = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_sed.bed")
    output:
        ASE_GFF_cut = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_cut.bed")
    message: "ASE BED with GFF gene IDs, cut."
    shell:
        """
        cut -f1,2,7,8,9,23 {input.ASE_GFF_sed} > {output.ASE_GFF_cut}
        """

rule sayres_ase_gff_grep:
    input:
        ASE_GFF_cut = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_cut.bed")
    output:
        ASE_GFF_grep = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_grep.bed")
    message: "ASE BED with GFF gene IDs, grep."
    shell:
        """
        grep "gene" {input.ASE_GFF_cut} | sed 's/Parent=//g' > {output.ASE_GFF_grep}
        """

rule sayres_ase_gff_ratio:                                                       
    input:                                                                      
        ASE_GFF_grep = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_grep.bed")
    output:                                                                     
        ASE_GFF_ratio = os.path.join(SAYRES_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_ratio.bed")
    message: "ASE BED with GFF gene IDs, ratio of REF to TOTAL."                
    shell:                                                                      
        """                                                                     
        python3 ASE_ratio.py --ase_filename {input.ASE_GFF_grep} --output_filename {output.ASE_GFF_ratio}
        """ 

# Clark 
rule clark_asereadcounter_gvcfs:
    input:
        REF = REF,
        BAM = lambda wildcards: Clark_RNA_SORTED_BAM_AL_DIR + wildcards.sample + "_sortedbycoord_wRGs_mkDups.bam",
        VCF = os.path.join(VCF_DIR, "SAYRES_CLARK_all_genotypes.vcf")
    output:
        ASE = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30.csv")
    params:
        xmx = "62g",
        xms = "62g",
        threads = "24",
        minDepth=30, 
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

rule clark_asereadcounter_gvcfs_HET:
    input:
        ASE = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30.csv")
    output:
        ASE_HET = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET.csv")
    params:
        xmx = "62g",
        xms = "62g",
        threads = "24",
        minDepth=30, 
        minMappingQuality=10,
        minBaseQuality=2
    message: "Filtering ASEReadCounter for HET sites."
    shell:
        """
        python3 getOnlyHeterozygousSites.py {input.ASE} {output.ASE_HET}
        """

rule clark_asereadcounter_gvcfs_HET_FD:
    input:
        Clark_FD = os.path.join(FD_DIR, "Clark_FD_NoHeader_sites.txt"),
        ASE_HET = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET.csv")        
    output:
        ASE_HET_FD = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.csv")
    message: "Filtering ASEReadCounter for FD sites."
    shell:
        """
        awk -F'\t' 'NR==FNR{{c[$1,$2]++;next}};c[$1,$2]>0' {input.Clark_FD} {input.ASE_HET} > {output.ASE_HET_FD}
        """


rule clark_asereadcounter_gvcfs_HET_FD_intersect:
    input:
        Intersect_FD = os.path.join(FD_DIR, "Sayres_Clark_intersect.txt"),
        ASE_HET = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET.csv")        
    output:
        ASE_HET_FD_inter = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_intersect.csv")
    message: "Filtering ASEReadCounter for FD sites."
    shell:
        """
        awk -F'\t' 'NR==FNR{{c[$1,$2]++;next}};c[$1,$2]>0' {input.Intersect_FD} {input.ASE_HET} > {output.ASE_HET_FD_inter}
        """
  
        
rule clark_ase_to_bed:
    input:
        ASE_HET_FD = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.csv")    
    output:
        ASE_to_BED = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.bed")
    message: "convert ASE to BED."
    shell:
        """
        python3 convert_ase_to_bed_fmt.py --ase_filename {input.ASE_HET_FD} --output_filename {output.ASE_to_BED}
        """
        
rule clark_ase_to_bed_gff_intersect:
    input:
        ASE_HET_FD_BED = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD.bed")
    output:
        ASE_to_BED_GFF = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF.bed")
    message: "ASE BED with GFF gene IDs."
    shell:
        """
        bedtools intersect  -wa -wb -a {input.ASE_HET_FD_BED} -b Reference/GCF_009193385.2_Nvit_psr_1.1_genomic.gff > {output.ASE_to_BED_GFF}
        """

rule clark_ase_gff_sed:
    input:
        ASE_to_BED_GFF = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF.bed")
    output:
        ASE_GFF_sed = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_sed.bed")
    message: "ASE BED with GFF gene IDs, sed;."
    shell:
        """
        sed 's/;/\t/g' {input.ASE_to_BED_GFF} > {output.ASE_GFF_sed}
        """

rule clark_ase_gff_cut:
    input:
        ASE_GFF_sed = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_sed.bed")
    output:
        ASE_GFF_cut = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_cut.bed")
    message: "ASE BED with GFF gene IDs, cut."
    shell:
        """
        cut -f1,2,7,8,9,23 {input.ASE_GFF_sed} > {output.ASE_GFF_cut}
        """

rule clark_ase_gff_grep:
    input:
        ASE_GFF_cut = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_cut.bed")
    output:
        ASE_GFF_grep = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_grep.bed")
    message: "ASE BED with GFF gene IDs, grep."
    shell:
        """
        grep "gene" {input.ASE_GFF_cut} | sed 's/Parent=//g' > {output.ASE_GFF_grep}
        """
        
rule clark_ase_gff_ratio:
    input:
        ASE_GFF_grep = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_grep.bed")
    output:
        ASE_GFF_ratio = os.path.join(CLARK_ASE_DIR, "{sample}_ASEReadCounter_minDepth30_HET_FD_GFF_ratio.bed")
    message: "ASE BED with GFF gene IDs, ratio of REF to TOTAL."
    shell:
        """
        python3 ASE_ratio.py --ase_filename {input.ASE_GFF_grep} --output_filename {output.ASE_GFF_ratio}
        """



