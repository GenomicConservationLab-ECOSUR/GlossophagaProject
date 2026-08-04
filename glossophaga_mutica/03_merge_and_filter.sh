#!/bin/bash
# Concatenate per-scaffold BCFs, then split out G. mutica-only callset (SNVs + indels retained)

OUTPUT_DIR="split"
FINAL_OUT="all.bcf"
STATS_OUT="mutica_plus_stats.txt"

ls -v ${OUTPUT_DIR}/*.bcf > bcf_to_combine.txt
bcftools concat -f bcf_to_combine.txt -Ob -o "$FINAL_OUT"
bcftools index "$FINAL_OUT"
bcftools stats "$FINAL_OUT" > "$STATS_OUT"

# Remove outgroups (G. leachii, G. commissarisi) to keep G. mutica only
bcftools view -s ^Morenoi.bam,Ocosingo.bam "$FINAL_OUT" -Ob -o mutica.bcf
bcftools index mutica.bcf
bcftools stats mutica.bcf > mutica_stats.txt

# High-confidence biallelic SNPs: Q>30, MAC>=1 (relaxed to keep singletons)
bcftools filter -i 'TYPE="snp" && QUAL>30 && MAC>=1' mutica.bcf | \
bcftools view -Oz -o mutica_filtered.vcf.gz
bcftools index mutica_filtered.vcf.gz
bcftools stats mutica_filtered.vcf.gz > mutica_filtered_stats.txt
