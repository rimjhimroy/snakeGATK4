#!/usr/bin/env python
from datetime import datetime
import os

"""
This workflow involves pre-processing the raw sequence data (uBAM format) to produce analysis-ready BAM files. Recalibrate Base Quality Scores is not included since our organism does not have know SNPs database.
"""
VCF,=glob_wildcards("analysis/gvcf2/{sample}.g.vcf.gz")

#today = datetime.today().strfmt('%Y-%m-%d')

rule WriteSampleMap:
    input:
        config["META"]
    output:
        "sample_map.txt"
    run:
        import pandas as pd
        sam = pd.read_csv(config["META"],sep='\t').set_index(["sample"], drop=False)
        sam['gvcf']=["analysis/gvcf2/"+y+".g.vcf.gz" for y in sam.index]
        fsam = sam[['sample','gvcf']].drop_duplicates()
        fsam.to_csv(path_or_buf="sample_map.txt",columns=['sample','gvcf'],sep='\t',index=False,header=False)


rule GenomicsDBImport:
    input:
        "sample_map.txt",
        lambda wildcards: expand("analysis/gvcf2/{vcf}.g.vcf.gz", vcf=samples["sample"].drop_duplicates())
    output:
        tar = "analysis/genomicsDB/{scaffold}.db.tar"
    
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output.db} \
        --L {wildcards.scaffold} \
        --sample-name-map sample_map.txt \
        --reader-threads 5 \
        --batch-size 10 \
        --tmp-dir ./tmp

        tar -cf {output.tar} {output.db}
	rm -r {output.db}
        """

rule GenotypeGVCFs:
    input:
        "analysis/genomicsDB/{scaffold}.db.tar"
    output:
        "analysis/vcf/{scaffold}.vcf.gz"
    benchmark:
        "benchmarks/GenotypeGVCFs/GenotypeGVCFs_{scaffold}.json"
    shell:
        """
	tar -xf {scaffold}.db.tar -C {scaffold}.db
        gatk --java-options "-Xmx8g -Xms8g" \
        GenotypeGVCFs \
        -R {REFERENCE} \
        -O {output} \
        --only-output-calls-starting-in-intervals \
        --include-non-variant-sites \
        -V gendb://{input} \
        -L {wildcards.scaffold} \
        --tmp-dir ./tmp

        rm -r {input}
        """

rule GatherVcfs:
    input:
        expand("analysis/vcf/{scaffold}.vcf.gz",scaffold = scaffolds)
    output:
        file="analysis/vcf/all_raw.vcf.gz",
        index="analysis/vcf/all_raw.vcf.gz.tbi"
    params:
        " -I ".join("analysis/vcf/" + s + ".vcf.gz" for s in scaffolds)
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
        -I {output.file} \
        -O {output.index}
        """
    
