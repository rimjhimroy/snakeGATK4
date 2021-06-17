#!/usr/bin/env python

def get_read_group(wildcards):
    flowcell_barcode = samples.loc[(wildcards.prefix),"flowcell_barcode"]
    lane = samples.loc[(wildcards.prefix),"lane"]
    sample_barcode = samples.loc[(wildcards.prefix),"sample_barcode"]
    sampleid = samples.loc[(wildcards.prefix),"sample"]
    return {"flowcell_barcode":flowcell_barcode,"lane":lane,"sample_barcode":sample_barcode,"sampleid":sampleid}


rule fastq2ubam:
    input:
        unpack(get_fastq2)
    output:
        bam="data/ubam/{prefix}.unaligned_reads.bam",
	    bai="data/ubam/{prefix}.unaligned_reads.bam.bai"
    benchmark:
        "benchmarks/fastq2ubam/fastq2ubam_{prefix}.json"
    params:
        flowcell_barcode = lambda wildcards : get_read_group(wildcards)['flowcell_barcode'],
        lane=lambda wildcards : get_read_group(wildcards)['lane'],
        sample_barcode = lambda wildcards : get_read_group(wildcards)['sample_barcode'],
        sampleid = lambda wildcards : get_read_group(wildcards)['sampleid'],
	platform = config['platform']
    shell:
       """ 
        picard FastqToSam \
         FASTQ={input.r1} \
         FASTQ2={input.r2} \
         OUTPUT={output.bam} \
         READ_GROUP_NAME={params.flowcell_barcode}{params.lane} \
         PLATFORM=illumina \
         PLATFORM_UNIT={params.flowcell_barcode}.{params.lane}.{params.sample_barcode} \
         LIBRARY_NAME={params.sampleid}_{params.sample_barcode} \
         SAMPLE_NAME={params.sampleid} \
         PLATFORM_MODEL={params.platform} \
         CREATE_INDEX=True \
         TMP_DIR=./tmp
	samtools index {output.bam}
	
        """

