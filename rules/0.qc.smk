#!/usr/bin/env python


rule fastqc:
    input:
        unpack(get_fastq)
    output:
        fwd = "analysis/qc/fastqc/{prefix}" + config["suffix"][0] + "_fastqc.zip",
        rev = "analysis/qc/fastqc/{prefix}" + config["suffix"][1] + "_fastqc.zip"
    threads: config["fastqc_threads"]
    shell:
        """
        fastqc --quiet --outdir analysis/qc/fastqc --noextract -f fastq {input} -t {threads}
        """
# this wrapper has a problem with multiple threads, it will use pigz if compressed file provided.Not available on this machine. Not sure about the reason. TODO:remove wrapper.
rule trim_pe:
    input:
        unpack(get_fastq)
    output:
        r1 = "data/trimmed/{prefix}_R1_001.trim.fastq.gz",
        r2 = "data/trimmed/{prefix}_R2_001.trim.fastq.gz",
        r1_unpaired="data/trimmed/{prefix}_R1_001.unpaired.fastq.gz",
        r2_unpaired="data/trimmed/{prefix}_R2_001.unpaired.fastq.gz"
    log:
        "logs/trimmomatic/{prefix}.log"
    benchmark:
        "benchmarks/trimmomatic/trimmomatic_{prefix}.json"
    params:
        trimmer=["LEADING:3","TRAILING:3","SLIDINGWINDOW:4:15","MINLEN:25"],
        compression_level="-9"
    wrapper:
        "0.35.2/bio/trimmomatic/pe"

rule bamstats:
    input:
        "analysis/mapping/{prefix}_aligned_duplicates_marked_sorted.bam" 
    output:
        "analysis/qc/bamtools/{prefix}_bamtools.stats"
    benchmark:
        "benchmarks/bamstats/bamstats_{prefix}.json"
    shell:
        """
        bamtools stats -in {input} | grep -v "*" > {output}
        
        """

rule multiqc:
    input:
        expand("analysis/qc/fastqc/{prefix}{R}_fastqc.zip",prefix=samples.index,R=config["suffix"]),
        expand("analysis/qc/bamtools/{prefix}_bamtools.stats",prefix=samples.index),
        expand("analysis/mapping/{prefix}.duplicate_metrics",prefix=samples.index)
    output:
        "analysis/qc/multiqc.html"
    benchmark:
        "benchmarks/multiqc/multiqc.json"
    wrapper:
        "0.35.1/bio/multiqc"

rule genomeCov:
    input:
        "analysis/mapping/{prefix}_aligned_duplicates_marked_sorted.bam"
    output:
        "analysis/qc/coverage/{prefix}_coverage.txt"
    benchmark:
        "benchmarks/genomeCov/genomeCov_{prefix}.json"
    shell:
        """
        bamcov {input} -o {output}
        """

