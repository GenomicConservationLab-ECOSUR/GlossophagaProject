#!/bin/bash
# Convert candidate windows to BED and intersect with the reference genome's
# structural annotation. Requires bedtools and the GFF3 annotation file
# (GCA_039655065.1_mGloMut1.hap1_genomic.gff -- same reference sequence as
# PRJNA1112473; see Methods for accession correspondence).

tail -n +2 Island_Robust_Candidates_Top200.csv | \
  awk -F, '{print $1 "\t" $2 "\t" $3 "\tIsland_" $1 "_" $2 "_" $3}' > candidates.bed

awk '$3 == "gene"' GCA_039655065.1_mGloMut1.hap1_genomic.gff > mutica_genes_island.gff

bedtools intersect -a candidates.bed -b mutica_genes_island.gff -wa -wb > annotated_island_windows.txt

awk -F'\t' '{print $4 "\t" $NF}' annotated_island_windows.txt | \
  sed -E 's/ID=[^;]+;//; s/Name=[^;]+;//; s/description=//; s/;gbkey.*//' > clean_island_annotations.txt
