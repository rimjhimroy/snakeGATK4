#!/usr/bin/env python
from datetime import datetime
import os


def get_ploidy(wildcards):
    ploidy = pd.read_csv(config["META"],sep='\t').set_index(["sample"], drop=False)[['sample','ploidy']].drop_duplicates()
    print("wildcard is %s ploidy is %s"%(wildcards.sample,ploidy.loc[(wildcards.sample), "ploidy"]))
    ploidy = ploidy.loc[(wildcards.sample), "ploidy"]
    return {str(ploidy)}

"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""

# base recalibration
rule BQSR_Pass1:         
     input:
        bam = "analysis/merged_samples/{sample}_aligned_duplicates_marked_sorted_merged.bam",
        Indels = "analysis/vcf/all_filtered_indel.vcf.gz",
        DbSNP = "analysis/vcf/all_filtered_snp.vcf.gz",
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
        fasta = REFERENCE
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
        Indels = "analysis/vcf/all_filtered_indel.vcf.gz",
        DbSNP = "analysis/vcf/all_filtered_snp.vcf.gz",
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
rule ApplyBQSR:
     input:
        bam = "analysis/BQSR1/{sample}_recal.pass1.bam",
        fasta = REFERENCE
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

#today = datetime.today().strfmt('%Y-%m-%d')

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
    threads: 8
    params:
        ploidy = get_ploidy
    shell:
        """
        gatk --java-options "-Xmx20G -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        --pair-hmm-implementation AVX_LOGLESS_CACHING_OMP \
        --native-pair-hmm-threads {threads} \
        -R {input.genome} \
        -I {input.bam} \
        -L scaffolds.list \
        -O {output} \
        -contamination 0 -ERC GVCF \
        -ploidy {params.ploidy} \
        --tmp-dir ./tmp
        """

rule GenomicsDBImport2:
    input:
        "sample_map.txt",
        lambda wildcards: expand("analysis/gvcf2/{vcf}.g.vcf.gz", vcf=samples["sample"].drop_duplicates())
    output:
        db = directory("analysis/genomicsDB2/{scaffold}.db"),
        tar = "analysis/genomicsDB2/{scaffold}.db.tar"
    
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output.db} \
        --L {wildcards.scaffold} \
        --sample-name-map sample_map.txt \
        --reader-threads 5 \
        --batch-size 50 \
        --tmp-dir=./tmp

        tar -cf {output.tar} {output.db}
        """

rule GenotypeGVCFs2:
    input:
        "analysis/genomicsDB2/{scaffold}.db"
    output:
        "analysis/vcf2/{scaffold}.vcf.gz"
    benchmark:
        "benchmarks/GenotypeGVCFs2/GenotypeGVCFs_{scaffold}.json"
    shell:
        """
        gatk --java-options "-Xmx8g -Xms8g" \
        GenotypeGVCFs \
        -R {REFERENCE} \
        -O {output} \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        --include-non-variant-sites \
        -V gendb://{input} \
        -L {wildcards.scaffold} \
        --TMP_DIR ./tmp
        """

rule GatherVcfs2:
    input:
        expand("analysis/vcf2/{scaffold}.vcf.gz",scaffold = scaffolds)
    output:
        file="analysis/vcf2/all_raw.vcf.gz",
        index="analysis/vcf2/all_raw.vcf.gz.tbi"
    params:
        " -I ".join("analysis/vcf2/" + s + ".vcf.gz" for s in scaffolds)
    benchmark:
        "benchmarks/GatherVcfs/GatherVcfs.json"
    shell:
        """
        gatk --java-options "-Xmx6g -Xms6g" \
        GatherVcfsCloud \
        --ignore-safety-checks \
        --gather-type BLOCK \
        -I {params} \
        -O {output.file}

        gatk --java-options "-Xmx6g -Xms6g" \
        IndexFeatureFile \
        -F {output.file} \
        -O {output.index}
        """

rule selectSNPs2:
    input:
        "analysis/vcf2/all_raw.vcf.gz"
    output:
        "analysis/vcf2/all_raw_snp.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type SNP \
        -select-type NO_VARIATION \
        -O {output}
        """


rule selectINDELs2:
    input:
        "analysis/vcf2/all_raw.vcf.gz"
    output:
        "analysis/vcf2/all_raw_indel.vcf.gz"
    shell:
        """
        gatk SelectVariants \
        -V {input} \
        -select-type INDEL \
        -O {output}
        """

rule hardfilterSNPs2:
    input:
        "analysis/vcf2/all_raw_snp.vcf.gz"
    output:
        "analysis/vcf2/all_filtered_snp.vcf.gz"
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
        --set-filtered-genotype-to-no-call true \
        -O {output}
        """

rule hardfilterINDELs2:
    input:
        "analysis/vcf2/all_raw_indel.vcf.gz"
    output:
        "analysis/vcf2/all_filtered_indel.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
        -V {input} \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        --set-filtered-genotype-to-no-call true \
        -O {output}
        """

