#!/bin/bash
# Author: Jesús Rocamontes-Morales
# Date: 01/07/2026
# Description: RDA
# Pipeline: RDA

BASE="Your directory"
VCF="$BASE/mutica_filtered.vcf.gz" # From the filtering step.
PLINK="$BASE/plink" # Use same version of Plink (PLINK v1.9.0-b.7.7)

# Step 1: Extract biallelic SNPs for IBS distance
bcftools view -v snps -m2 -M2 $VCF -Oz -o $BASE/mutica_rda.vcf.gz
bcftools index $BASE/mutica_rda.vcf.gz

# Step 2: Convert to PLINK and calculate IBS distance matrix
$PLINK --vcf $BASE/mutica_rda.vcf.gz \
    --double-id \
    --allow-extra-chr \
    --set-missing-var-ids @:# \
    --make-bed \
    --out $BASE/mutica_rda

$PLINK --bfile $BASE/mutica_rda \
    --allow-extra-chr \
    --distance 1-ibs square \
    --out $BASE/mutica_ibs
