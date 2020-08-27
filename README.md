# snakemake-gatk4-non-model
#### germline variant calling
Working with non-model organism means you don't have known SNPs and prepared interval list in [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle). Alternatively I did:
- **Optional Base Quality Score Calibration(BQSR)** (working...)
- Apply hard filters to call sets.
- Split genome into scaffolds as intervals.

#### Note
When *not* using the latest patterned flow cell data (usually illuma sequencers which are older versions than the HiSeq 3000/HiSeq 4000 Systems or the NovaSeq 6000 System [link](https://emea.illumina.com/science/technology/next-generation-sequencing/sequencing-technology/patterned-flow-cells.html), remove `--OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \` line from `MarkDuplicates` rule in `2.preprocessing.smk` which will revert it to default value of `100`.  
When working with PCR-free libraries (e.g., trueseq), add `--pcr-indel-model NONE \` in rule `HaplotypeCaller` in `2.preprocessing.smk`.  

<img src="https://github.com/bakeronit/snakemake-gatk4-non-model/blob/master/gatk-snakemake.png" alt="workflow" width="420" height="494">

**How to run:**

0. `conda env create -f environment.yaml` and install [bamcov](https://github.com/fbreitwieser/bamcov) ...
1. prepare metadata.tsv.
2. modify config.yaml if needed.
3. create `output/slurm_out` and `tmp` folders
4. `snakemake -n` 

**To run on Slurm HPC:**

`nohup snakemake -p --cluster-config slurm_hpc.json --cluster "sbatch -J {cluster.name} -t {cluster.time} -c {cluster.cpu} --mem={cluster.mem} -p {cluster.partition} -o {cluster.output} -e {cluster.error}" --jobs 100 --latency-wait 5 &`

**To run on trimed reads without QC:**
`snakemake -s SnakeGATK1 -p -j 9999 --cluster-config slurm_hpc.json --cluster "sbatch -J {cluster.name} -t {cluster.time} -c {cluster.cpu} --mem={cluster.mem} -p {cluster.partition} -o {cluster.output} -e {cluster.error}" --latency-wait 120`

**To submit everything in one go:**
`snakemake -s SnakeGATK1 -p -j 9999 --immediate-submit --notemp --until HaplotypeCaller --cluster-config slurm_hpc.json --cluster "./parseJobID.py {dependencies}" --latency-wait 120`

**Clean everything:**

`snakemake clean`


Reference:

https://snakemake.readthedocs.io/en/stable/
https://github.com/gatk-workflows/gatk4-germline-snps-indels
https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling
https://zhuanlan.zhihu.com/p/33891718
https://software.broadinstitute.org/gatk/documentation/article?id=11097
https://gatkforums.broadinstitute.org/gatk/discussion/12443/genomicsdbimport-run-slowly-with-multiple-samples
