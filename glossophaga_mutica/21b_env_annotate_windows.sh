#!/bin/bash
# Convert each cluster's candidate windows to BED and intersect with the
# reference genome's structural annotation (gene features only).

GFF_FILE="GCA_039655065.1_mGloMut1.hap1_genomic.gff"
CLUSTERS=("Concha" "Hobo_Ixta" "Nizanda" "Poana" "Quilamula")

awk '$3 == "gene"' "$GFF_FILE" > mutica_genes_cluster.gff

for cluster in "${CLUSTERS[@]}"; do
    CANDIDATES_CSV="Cluster_${cluster}_Robust_Candidates_Top200.csv"
    if [ -f "$CANDIDATES_CSV" ]; then
        tail -n +2 "$CANDIDATES_CSV" | \
          awk -F, -v cl="$cluster" '{print $1 "\t" $2 "\t" $3 "\t" cl "_" $1 "_" $2 "_" $3}' > "${cluster}_candidates.bed"

        bedtools intersect -a "${cluster}_candidates.bed" -b mutica_genes_cluster.gff -wa -wb > "annotated_${cluster}_windows.txt"

        awk -F'\t' '{print $4 "\t" $NF}' "annotated_${cluster}_windows.txt" | \
          sed -E 's/ID=[^;]+;//; s/Name=[^;]+;//; s/description=//; s/;gbkey.*//' > "clean_${cluster}_annotations.txt"

        rm "${cluster}_candidates.bed"
    fi
done
