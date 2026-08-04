#!/bin/bash
# PLINK conversion, LD pruning, and PCA on the filtered G. mutica SNP set
# For our analysis we used: PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)

./plink --vcf mutica_filtered.vcf.gz \
  --double-id --allow-extra-chr \
  --biallelic-only strict \
  --set-missing-var-ids @:# \
  --make-bed --out mutica_full

./plink --bfile mutica_full \
  --allow-extra-chr \
  --indep-pairwise 50 10 0.2 \
  --out mutica_pruning

./plink --bfile mutica_full \
  --allow-extra-chr \
  --extract mutica_pruning.prune.in \
  --make-bed --out mutica_final_pruned_popgen

# Fix chromosome IDs for ADMIXTURE compatibility
awk '{$1=1; print $0}' mutica_final_pruned_popgen.bim > tmp && mv tmp mutica_final_pruned_popgen.bim

./plink --bfile mutica_final_pruned_popgen \
  --allow-extra-chr \
  --pca 10 \
  --out mutica_pca_results
