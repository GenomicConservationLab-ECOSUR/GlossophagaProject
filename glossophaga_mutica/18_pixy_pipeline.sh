#!/bin/bash
# Prepare input and run pixy (Dxy, pi) in 200 kb sliding windows (50 kb step)
# for two population designs: (1) five environmentally-defined mainland
# clusters, each vs. all others, and (2) Island vs. Mainland.
# Requires bcftools, tabix, bedtools, pixy.

bcftools view -O z -o mutica_pixy_ready.vcf.gz mutica.bcf # from the SNPs pipeline.
tabix -p vcf mutica_pixy_ready.vcf.gz

# Environmental cluster map (Nayarit/Belice excluded; see PCA/UPGMA
# clustering script for how clusters were assigned)
cat << EOF > pop_map_env.txt
Concha.bam	Cluster_Concha
Poana.bam	Cluster_Poana
Quilamula.bam	Cluster_Quilamula
Hobonil.bam	Cluster_Hobo_Ixta
Ixtapa.bam	Cluster_Hobo_Ixta
Nizanda.bam	Cluster_Nizanda
EOF

# Island vs. Mainland map (Belice excluded)
cat << EOF > pop_map_island.txt
Nayarit.bam	Island
Concha.bam	Mainland
Poana.bam	Mainland
Quilamula.bam	Mainland
Hobonil.bam	Mainland
Ixtapa.bam	Mainland
Nizanda.bam	Mainland
EOF

# Sliding windows (200 kb, 50 kb step) from VCF contig lengths
bcftools view -h mutica_pixy_ready.vcf.gz | grep "^##contig=" | \
  sed 's/.*ID=\([^,]*\),length=\([0-9]*\).*/\1\t\2/' > genome.lengths
bedtools makewindows -g genome.lengths -w 200000 -s 50000 > sliding_windows_200k_50k.bed

# Run pixy: each environmental cluster vs. all others
mkdir -p output_environmental
for TARGET in Cluster_Concha Cluster_Poana Cluster_Quilamula Cluster_Hobo_Ixta Cluster_Nizanda; do
    awk -v tgt="$TARGET" '{if ($2 == tgt) print $1 "\t" tgt; else print $1 "\tOthers"}' \
        pop_map_env.txt > "pop_map_${TARGET}_vs_others.txt"
    pixy --stats dxy pi \
         --vcf mutica_pixy_ready.vcf.gz \
         --populations "pop_map_${TARGET}_vs_others.txt" \
         --bed_file sliding_windows_200k_50k.bed \
         --output_prefix "${TARGET}_vs_others" \
         --output_folder output_environmental \
         --bypass_invariant_check \
         --n_cores 8
done

# Run pixy: Island vs. Mainland
mkdir -p output_island
pixy --stats dxy pi \
     --vcf mutica_pixy_ready.vcf.gz \
     --populations pop_map_island.txt \
     --bed_file sliding_windows_200k_50k.bed \
     --output_prefix "Island_vs_Mainland" \
     --output_folder output_island \
     --bypass_invariant_check \
     --n_cores 8
