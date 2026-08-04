#!/bin/bash
# Align each PCG gene folder independently using MAFFT.
# Requires MAFFT. Run from inside the extracted_pcgs folder produced by
# 08a_extract_pcgs.py (one subfolder per gene, e.g. ATP6/, COX1/, ...).

for gene_folder in */; do
    gene_name=${gene_folder%/}
    echo "-----------------------------------"
    echo "Aligning $gene_name..."
    cd "$gene_name" || continue
    cat *.fasta > combined.fasta
    mafft --auto combined.fasta > "aligned_${gene_name}.fasta"
    rm combined.fasta
    cd ..
done
