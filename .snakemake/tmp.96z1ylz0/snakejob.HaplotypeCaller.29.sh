#!/bin/sh
# properties = {"type": "single", "rule": "HaplotypeCaller", "local": false, "input": ["data/genome/biscutV2.1_new.fa", "analysis/merged_samples/pg42_aligned_duplicates_marked_sorted_merged.bam"], "output": ["analysis/gvcf/pg42_raw_variants_1.vcf.gz"], "wildcards": {"sample": "pg42"}, "log": [], "threads": 12, "resources": {}, "jobid": 29, "cluster": {"account": "rchoudhury", "time": "4-0", "cpu": 12, "partition": "all", "mem": "100G", "name": "HaplotypeCaller.pg42", "output": "output/slurm_out/HaplotypeCaller_pg42_%j.out", "error": "output/slurm_out/HaplotypeCaller_pg42_%j.err"}}
 cd /gpfs/homefs/ips/rchoudhury/snakemake-gatk4-non-model && \
/home/ubelix/ips/rchoudhury/miniconda3/envs/gatk4/bin/python3.6 \
-m snakemake analysis/gvcf/pg42_raw_variants_1.vcf.gz --snakefile /gpfs/homefs/ips/rchoudhury/snakemake-gatk4-non-model/SnakeGATK1 \
--force -j --keep-target-files --keep-remote \
--wait-for-files /gpfs/homefs/ips/rchoudhury/snakemake-gatk4-non-model/.snakemake/tmp.96z1ylz0 data/genome/biscutV2.1_new.fa analysis/merged_samples/pg42_aligned_duplicates_marked_sorted_merged.bam --latency-wait 120 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules HaplotypeCaller --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /gpfs/homefs/ips/rchoudhury/snakemake-gatk4-non-model/.snakemake/tmp.96z1ylz0/29.jobfinished || (touch /gpfs/homefs/ips/rchoudhury/snakemake-gatk4-non-model/.snakemake/tmp.96z1ylz0/29.jobfailed; exit 1)

