conda activate gatk4

# copy /gpfs/homefs/ips/rchoudhury/miniconda3/envs/gatk4/share/snpeff-4.5covid19-1/snpEff.config to current folder
Add the following:
```
# Biscutella genome, Version 2.1
biscutella.genome: biscutV2.1_new.fa
```
mkdir -p data/biscutella
# move the gff to data/biscutella/genes.gff
# move the genome do data/biscutella/sequences.fa

snpEff build -gff3 -v biscutella




java -jar -Xmx4G ~/miniconda3/envs/gatk4/share/snpeff-4.5covid19-1/snpEff.jar eff biscutella ../analysis/vcf_final/all_filtered_snp.vcf.gz > all_filtered_snp_snpEff.vcf

cd Eff
snpeff="~/miniconda3/envs/gatk4/share/snpeff-4.5covid19-1"
vcf="all_filtered_snp_snpEff"
scripts="/home/ubelix/ips/rchoudhury/snakemake-gatk4-non-model/Scripts/popgen"

conda activate gatk4
# genic (includes 5', 3' UTRs)
SnpSift filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" ${vcf}.vcf > ${vcf}_gene.vcf
# coding
SnpSift filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" ${vcf}.vcf > ${vcf}_coding.vcf
# non-synonymous
SnpSift filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" ${vcf}.vcf > ${vcf}_nonsyn.vcf
# synonymous
SnpSift filter "(ANN[0].EFFECT has 'synonymous_variant')" ${vcf}.vcf > ${vcf}_syn.vcf
# Four-fold degenrate sites (output file suffix: 4fd)
python $scripts/summary_stats/parse_snpeff_synonymous.py ${vcf}_syn.vcf
