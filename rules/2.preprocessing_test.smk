#!/usr/bin/env python

def get_ploidy(wildcards):
    ploidy = pd.read_csv(config["META"],sep='\t').set_index(["sample"], drop=False)[['sample','ploidy']].drop_duplicates()
    print("wildcard is %s ploidy is %s"%(wildcards.sample,ploidy.loc[(wildcards.sample), "ploidy"]))
    ploidy = ploidy.loc[(wildcards.sample), "ploidy"]
    return {str(ploidy)}

def get_pcr(wildcards):
    pcr = pd.read_csv(config["META"],sep='\t').set_index(["sample"], drop=False)[['sample','pcr_indel_model']].drop_duplicates()
    print("wildcard is %s pcr_indel_model is %s"%(wildcards.sample,pcr.loc[(wildcards.sample), "pcr_indel_model"]))
    pcr = pcr.loc[(wildcards.sample), "pcr_indel_model"]
    return {str(pcr)}

def get_pixeldist(wildcards):
    pixeldist = pd.read_csv(config["META"],sep='\t').set_index(["prefix"], drop=False)[['prefix','flow_cell_type']].drop_duplicates()
    print("wildcard is %s flow_cell_type is %s"%(wildcards.prefix,pixeldist.loc[(wildcards.prefix), "flow_cell_type"]))
    pixeldist = pixeldist.loc[(wildcards.prefix), "flow_cell_type"]
    return {str(pixeldist)}


"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""

#Part 2: separate raw SNPs and INDELs in two separated vcf files since filtration parameters will be different for each (see further)
rule BQSRselectSNPs:
    input:
        "analysis/gvcf/{sample}_raw_variants_1.vcf.gz"
    output:
        "analysis/gvcf/{sample}_raw_snp_1.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type SNP \
        -O {output}
        """


rule BQSRselectINDELs:
    input:
        "analysis/gvcf/{sample}_raw_variants_1.vcf.gz"
    output:
        "analysis/gvcf/{sample}_raw_indel_1.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type INDEL \
        -O {output}
        """

#Part 3: hard-filter raw SNPs and INDELs vcf files with [VariantFiltration] to mark the bad ones with FILTER (many variants in both raw files are not genuine variants)
#SNPs and INDELs matching any of the specified criteria will be considered bad and marked FILTER
#for _raw_SNPs_1

rule BQSRhardfilterSNPs:
    input:
        "analysis/gvcf/{sample}_raw_snp_1.vcf.gz"
    output:
        "analysis/gvcf/{sample}_raw_snp_2.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
        -V {input} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O {output}
        """

#for _raw_ INDELs_1
rule BQSRhardfilterINDELs:
    input:
        "analysis/gvcf/{sample}_raw_indel_1.vcf.gz"
    output:
        "analysis/gvcf/{sample}_raw_indel_2.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
        -V {input} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O {output}
        """

#exclude filtered
rule BQSR2selectSNPs:
    input:
        "analysis/gvcf/{sample}_raw_snp_2.vcf.gz"
    output:
        "analysis/gvcf/{sample}_filtered_snp_2.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        --exclude-filtered \
        -V {input} \
        -O {output}
        """

rule BQSR2selectINDELs:
    input:
        "analysis/gvcf/{sample}_raw_indel_2.vcf.gz"
    output:
        "analysis/gvcf/{sample}_filtered_indel_2.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        --exclude-filtered \
        -V {input} \
        -O {output}
        """


# base recalibration
rule BQSR_Pass1:         
    input:
        bam = "analysis/merged_samples/{sample}_aligned_duplicates_marked_sorted_merged.bam",
        Indels = "analysis/gvcf/{sample}_filtered_indel_2.vcf.gz",
        DbSNP = "analysis/gvcf/{sample}_filtered_snp_2.vcf.gz",
        fasta = REFERENCE
    output:
        Recall =  "analysis/Recal1/{sample}_recal.table"
    resources:
        mem_mb = 50000
    shell:"""
           gatk BaseRecalibrator \
           -I {input.bam} \
           -R {input.fasta} \
           --known-sites  {input.Indels}  \
           --known-sites  {input.DbSNP} \
           -O {output.Recall}
           """

rule ApplyBQSR:
    input:
        bam = "analysis/merged_samples/{sample}_aligned_duplicates_marked_sorted_merged.bam",
        fasta = REFERENCE,
        recal = "analysis/Recal1/{sample}_recal.table"
    output:
        Rbam = "analysis/BQSR1/{sample}_recal.pass1.bam"
    resources:
        mem_mb = 50000
    shell:"""
          gatk ApplyBQSR \
          -I {input.bam}  \
          -R {input.fasta} \
          --bqsr-recal-file {input.recal} \
          -O {output.Rbam}
          """

#Base Recalibration 
rule BQSR_Pass2:         
    input:
        bam = "analysis/BQSR1/{sample}_recal.pass1.bam",
        Indels = "analysis/gvcf/{sample}_filtered_indel_2.vcf.gz",
        DbSNP = "analysis/gvcf/{sample}_filtered_snp_2.vcf.gz",
        fasta = REFERENCE
    output:
        Recall =  "analysis/Recal2/{sample}_recal.table"
    resources:
        mem_mb = 50000
    shell:"""
           gatk BaseRecalibrator \
           -I {input.bam} \
           -R {input.fasta} \
           --known-sites  {input.Indels}  \
           --known-sites  {input.DbSNP} \
           -O {output.Recall}
           """ 

#detects systematic errors made by the sequencer when it estimates the quality score of each base call
rule ApplyBQSR2:
    input:
        bam = "analysis/BQSR1/{sample}_recal.pass1.bam",
        fasta = REFERENCE,
        recal = "analysis/Recal2/{sample}_recal.table"
    output:
        Rbam = "analysis/BQSR_2/{sample}_recal.pass2.bam"
    resources:
        mem_mb = 50000
    shell:"""
          gatk ApplyBQSR \
          -I {input.bam}  \
          -R {input.fasta} \
          --bqsr-recal-file {input.recal} \
          -O {output.Rbam}

          """  
#plots the two models of recalibration to see whether convergence is observed or not
rule AnalyzeCovariates:
    input:
        genome = REFERENCE,
        before = "analysis/Recal1/{sample}_recal.table",
        after = "analysis/Recal2/{sample}_recal.table"
    output:
        out="analysis/analyzeCovariated/{sample}_recalibration_1_2.pdf"
    threads: 8
    shell:
        """
        gatk AnalyzeCovariates \
        -before {input.before} -after {input.after} -plots {output.out}
        """

rule samtools_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    benchmark:
        "benchmarks/processbam/{prefix}.samtools_index.json"
    wrapper:
        "0.27.1/bio/samtools/index"


rule HaplotypeCaller2:
    input:
        genome = REFERENCE,
        bam = "analysis/BQSR_2/{sample}_recal.pass2.bam"
    output:
        out="analysis/gvcf2/{sample}.g.vcf.gz"
    benchmark:
        "benchmarks/HaplotypeCaller2/HaplotypeCaller2_{sample}.json"
    threads: 12
    params:
        ploidy = get_ploidy,
        pcr = get_pcr
    shell:
        """
        gatk --java-options "-Xmx20G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads {threads} \
        --pcr-indel-model {params.pcr} \
        -R {input.genome} \
        -I {input.bam} \
        -L scaffolds.list \
        -O {output} \
        -contamination 0 -ERC GVCF \
        -ploidy {params.ploidy} \
        --tmp-dir ./tmp
        """

