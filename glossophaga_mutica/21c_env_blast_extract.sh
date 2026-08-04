#!/bin/bash
# Extract gene-model sequences overlapping candidate windows and identify
# them via BLASTX against UniProtKB/Swiss-Prot, for each environmental
# cluster. Requires bedtools, BLAST+, and a local uniprot_db (see repo
# README for database build instructions).

GFF_FILE="GCA_039655065.1_mGloMut1.hap1_genomic.gff"
GENOME_FASTA="Glossophaga_mutica_autosomes.fna"
UNIPROT_DB="uniprot_db"
CLUSTERS=("Concha" "Hobo_Ixta" "Nizanda" "Poana" "Quilamula")

for cluster in "${CLUSTERS[@]}"; do
    CANDIDATES_CSV="Cluster_${cluster}_Robust_Candidates_Top200.csv"
    if [ -f "$CANDIDATES_CSV" ]; then
        echo "Processing ${cluster}..."

        # 1. Convert CSV to temporary BED
        tail -n +2 "$CANDIDATES_CSV" | awk -F, '{print $1 "\t" $2 "\t" $3}' > "${cluster}_temp.bed"

        # 2. Intersect windows with GFF to isolate true genes
        bedtools intersect -a "$GFF_FILE" -b "${cluster}_temp.bed" -wa | \
          awk '$3 == "gene"' > "${cluster}_isolated_genes.gff"

        # 3. Parse GFF to extract exact coordinates and standard gene names
        #    (Name= preferred, locus_tag= fallback for anonymous NCBI calls)
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
        }' "${cluster}_isolated_genes.gff" | sort -u > "${cluster}_exact_genes.bed"

        # 4. Extract pure gene sequences for BLAST
        bedtools getfasta -fi "$GENOME_FASTA" -bed "${cluster}_exact_genes.bed" \
          -nameOnly > "${cluster}_exact_genes.fasta"

        # 5. Run targeted BLASTX
        blastx -query "${cluster}_exact_genes.fasta" \
               -db "$UNIPROT_DB" \
               -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
               -max_target_seqs 1 \
               -num_threads 8 \
               -out "specific_genes_blast_${cluster}.txt"

        rm "${cluster}_temp.bed" "${cluster}_isolated_genes.gff" \
           "${cluster}_exact_genes.bed" "${cluster}_exact_genes.fasta"
        echo "  - Saved: specific_genes_blast_${cluster}.txt"
    fi
done
