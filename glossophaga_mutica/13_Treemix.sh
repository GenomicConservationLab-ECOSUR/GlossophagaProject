#!/bin/bash
# Prepare filtered, LD-pruned, population-labeled allele frequency input
# for TreeMix, then run TreeMix across m = 0 to 5 migration edges.
# Requires bcftools, PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024) or compatible,
# plink2treemix.py (distributed with TreeMix), and TreeMix v1.13.

# Filter for biallelic SNPs, quality > 30, site depth > 10, MAC >= 1,
# and < 10% missing data per site
bcftools view -v snps -m2 -M2 all.bcf | \
bcftools filter -i 'QUAL>30 && INFO/DP>10 && MAC>=1 && F_MISSING<0.1' -Oz -o mutica_treemix_filtered.vcf.gz
bcftools index mutica_treemix_filtered.vcf.gz

# LD pruning (window 50, step 10, r2 = 0.1)
./plink --vcf mutica_treemix_filtered.vcf.gz \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --indep-pairwise 50 10 0.1 \
  --out treemix_pruning

# Extract the pruned set (outgroups retained)
./plink --vcf mutica_treemix_filtered.vcf.gz \
  --allow-extra-chr \
  --set-missing-var-ids @:# \
  --extract treemix_pruning.prune.in \
  --make-bed --out mutica_tm_pruned

# Population/cluster assignments, matching the regional groupings used
# throughout the manuscript (Pacifico, Golfo, Yucatan, Oaxaca, Nayarit)
cat << CLUST > mutica.clust
Concha.bam Concha.bam Pacific
Poana.bam Poana.bam Gulf
Quilamula.bam Quilamula.bam Gulf
Hobonil.bam Hobonil.bam Yucatan
Ixtapa.bam Ixtapa.bam Pacific
Nizanda.bam Nizanda.bam Oaxaca
Nayarit.bam Nayarit.bam Nayarit
Morenoi.bam Morenoi.bam Leachii
Ocosingo.bam Ocosingo.bam Commissarisi
CLUST

# Allele frequencies by population
./plink --bfile mutica_tm_pruned \
  --allow-extra-chr \
  --within mutica.clust \
  --freq \
  --out mutica_tm_ready

# Convert to TreeMix input format
gzip mutica_tm_ready.frq.strat
python plink2treemix.py mutica_tm_ready.frq.strat.gz mutica_tm.input.gz

# NOTE: the header order below must exactly match the population column
# order plink2treemix.py actually writes -- verify with
# `zcat mutica_tm.input.gz | head -1` before trusting this header, since a
# mismatch will silently mislabel every population in all downstream output.
echo "Commissarisi Gulf Leachii Nayarit Oaxaca Pacific Yucatan" > header.txt
zcat mutica_tm.input.gz > temp_data.txt
cat header.txt temp_data.txt | gzip > mutica_tm_fixed.input.gz

# Verify
zcat mutica_tm_fixed.input.gz | head -n 2

# Run TreeMix for m = 0 to 5 migration edges
for i in {0..5}; do
  echo ">>> Starting TreeMix: m = $i"
  treemix -i mutica_tm_fixed.input.gz \
          -m $i \
          -root Leachii \
          -bootstrap \
          -k 500 \
          -noss \
          -o mutica_tm_m$i > tm_${i}.log
done
