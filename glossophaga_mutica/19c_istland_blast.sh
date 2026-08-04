#!/bin/bash
# Extract gene-model sequences overlapping candidate windows and identify
# them via BLASTX against UniProtKB/Swiss-Prot. Requires bedtools, BLAST+,
# and a local uniprot_db (see repo README for database build instructions).

bedtools intersect -a GCA_039655065.1_mGloMut1.hap1_genomic.gff -b candidates.bed -wa | \
  awk '$3 == "gene"' > isolated_genes.gff

# Parse GFF column 9 to extract a clean gene name (Name= preferred,
# locus_tag= fallback for anonymous NCBI automated gene calls)
awk -F'\t' '{
    split($9, attributes, ";");
    name = "Unknown";
    for (i in attributes) {
        if (attributes[i] ~ /^Name=/) {
            name = substr(attributes[i], 6);
            break;
        } else if (attributes[i] ~ /^locus_tag=/) {
            name = substr(attributes[i], 11);
        }
    }
    print $1 "\t" ($4-1) "\t" $5 "\t" name;
}' isolated_genes.gff | sort -u > exact_genes.bed

bedtools getfasta -fi Glossophaga_mutica_autosomes.fna -bed exact_genes.bed -nameOnly > exact_genes.fasta

blastx -query exact_genes.fasta \
       -db uniprot_db \
       -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
       -max_target_seqs 1 \
       -num_threads 8 \
       -out specific_genes_blast_island.txt
