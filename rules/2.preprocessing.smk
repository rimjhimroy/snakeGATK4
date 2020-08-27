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

rule SamToFastqAndBwaMem:
    input:
        ubam = "data/ubam/{prefix}.unaligned_reads.bam",
        reference=REFERENCE

    output:
        out=protected("analysis/mapping/{prefix}_unmerged.bam")
    benchmark:
        "benchmarks/SamToFastqAndBwaMem/SamToFastqAndBwaMem_{prefix}.json"
    log:
        "logs/{prefix}.bwa.stderr.log"
    threads: config["bwa_threads"]
    shell:
        """
        picard SamToFastq \
        INPUT={input.ubam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        TMP_DIR=./tmp \
        NON_PF=true | \
        bwa mem -K 100000000 -p -v 3 -t {threads} -Y {input.reference} /dev/stdin - 2> >(tee {log} >&2) | \
        samtools view -1 - > {output.out}
        """

bwa_version = subprocess.check_output("bwa 2>&1 | grep -e 'Version'", shell=True).decode("utf-8").rstrip()

rule MergeBamAlignment:
    input:
        unmerged_bam = "analysis/mapping/{prefix}_unmerged.bam",
        unmapped_bam = "data/ubam/{prefix}.unaligned_reads.bam",
        reference = REFERENCE
    output:
        temp(os.path.join(config["tmpdir"],"{prefix}_aligned_unsorted.bam"))
    benchmark:
        "benchmarks/MergeBamAlignment/MergeBamAlignment_{prefix}.json"
    params:
        v_bwa = bwa_version

    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms3000m" \
        MergeBamAlignment \
        --VALIDATION_STRINGENCY SILENT \
        --EXPECTED_ORIENTATIONS FR \
        --ATTRIBUTES_TO_RETAIN X0 \
        --ALIGNED_BAM {input.unmerged_bam}  \
        --UNMAPPED_BAM {input.unmapped_bam} \
        --OUTPUT {output} \
        --REFERENCE_SEQUENCE {input.reference} \
        --PAIRED_RUN true \
        --SORT_ORDER "unsorted" \
        --IS_BISULFITE_SEQUENCE false \
        --ALIGNED_READS_ONLY false \
        --CLIP_ADAPTERS false \
        --MAX_RECORDS_IN_RAM 2000000 \
        --ADD_MATE_CIGAR true \
        --MAX_INSERTIONS_OR_DELETIONS -1 \
        --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
        --PROGRAM_RECORD_ID "bwamem" \
        --PROGRAM_GROUP_VERSION "{params.v_bwa}" \
        --PROGRAM_GROUP_COMMAND_LINE "bwa mem -K 100000000 -p -v 3 -t 2 -Y {input.reference}" \
        --PROGRAM_GROUP_NAME "bwamem" \
        --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
        --ALIGNER_PROPER_PAIR_FLAGS true \
        --UNMAP_CONTAMINANT_READS true \
        --TMP_DIR ./tmp
        """


rule MarkDuplicates:
    input:
        os.path.join(config["tmpdir"],"{prefix}_aligned_unsorted.bam")
    output:
        bam = temp(os.path.join(config["tmpdir"],"{prefix}_aligned_unsorted_duplicates_marked.bam")),
        metrics_filename = "analysis/mapping/{prefix}.duplicate_metrics"
    params:
        pixeldist = get_pixeldist
    benchmark:
        "benchmarks/MarkDuplicates/MarkDuplicates_{prefix}.json"
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m -XX:+UseParallelGC -XX:ParallelGCThreads=2" \
        MarkDuplicates \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --METRICS_FILE {output.metrics_filename} \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE {params.pixeldist} \
        --ASSUME_SORT_ORDER "queryname" \
        --CREATE_MD5_FILE true \
        --TMP_DIR ./tmp
        """

rule SortAndFixTags:
    input:
        bam = os.path.join(config["tmpdir"],"{prefix}_aligned_unsorted_duplicates_marked.bam"),
        reference = REFERENCE
    output:
        out=protected("analysis/mapping/{prefix}_aligned_duplicates_marked_sorted.bam")
    benchmark:
        "benchmarks/SortAndFixTags/SortAndFixTags_{prefix}.json"
    shell:
        """
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms4000m" \
        SortSam \
        --INPUT {input.bam} \
        --OUTPUT /dev/stdout \
        --SORT_ORDER "coordinate" \
        --CREATE_INDEX false \
        --CREATE_MD5_FILE false \
        --TMP_DIR ./tmp \
        | \
        gatk --java-options "-Dsamjdk.compression_level=5 -Xms500m" \
        SetNmMdAndUqTags \
        --INPUT /dev/stdin \
        --OUTPUT {output.out} \
        --CREATE_INDEX true \
        --CREATE_MD5_FILE true \
        --REFERENCE_SEQUENCE {input.reference} \
        --TMP_DIR ./tmp
        """
rule MergeBamLanes:
    input:
        bamin = lambda wildcards: expand("analysis/mapping/{unit}_aligned_duplicates_marked_sorted.bam",unit=dunits[wildcards.sample]["prefix"]),
    output:
        out="analysis/merged_samples/{sample}_aligned_duplicates_marked_sorted_merged.bam"
    benchmark:
        "benchmarks/merge/{sample}.txt"
    params:
        inbam= lambda wildcards, input: " -I ".join(input.bamin)
    shell:
        """
        gatk MergeSamFiles \
        -I {params.inbam} \
        --ASSUME_SORTED true \
        -O {output.out} \
        --CREATE_INDEX true \
        --USE_THREADING
        """
#Part 1: generate a first list of raw SNPs and INDELs (${sample}_raw_variants_1.vcf) from the .bam file with GATK [HaplotypeCaller]
rule HaplotypeCaller:
    input:
        genome = REFERENCE,
        bam = "analysis/merged_samples/{sample}_aligned_duplicates_marked_sorted_merged.bam"
    output:
        out="analysis/gvcf/{sample}_raw_variants_1.vcf.gz"
    benchmark:
        "benchmarks/HaplotypeCaller/HaplotypeCaller_{sample}.json"
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

