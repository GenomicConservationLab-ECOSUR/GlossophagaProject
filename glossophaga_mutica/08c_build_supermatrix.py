#!/usr/bin/env python3
# Build the concatenated mitogenomic supermatrix and its gene partitions,
# then write both to a single NEXUS file for BEAUti/BEAST.
# Requires biopython. Run from the directory containing extracted_pcgs/.
#
# Combines three originally separate steps: partition generation,
# gap-padded concatenation, and NEXUS assembly. Gene order must be
# consistent across all steps (enforced here via sorted() gene names).

import os
from Bio import SeqIO

ROOT_DIR = "extracted_pcgs"
FASTA_OUT = "supermatrix.fasta"
PARTITIONS_OUT = "partitions.txt"
NEXUS_OUT = "supermatrix.nex"


def generate_partitions(root_dir):
    gene_lengths = {}
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.startswith("aligned_") and file.endswith(".fasta"):
                gene_name = file.replace("aligned_", "").replace(".fasta", "")
                for record in SeqIO.parse(os.path.join(root, file), "fasta"):
                    gene_lengths[gene_name] = len(record.seq)
                    break

    sorted_genes = sorted(gene_lengths.keys())
    current_pos = 1
    with open(PARTITIONS_OUT, "w") as f:
        for gene in sorted_genes:
            length = gene_lengths[gene]
            end_pos = current_pos + length - 1
            f.write(f"DNA, {gene} = {current_pos}-{end_pos}\n")
            current_pos = end_pos + 1

    print(f"partitions.txt written. Total alignment length: {current_pos - 1} bp")
    return sorted_genes, gene_lengths


def concatenate_with_padding(root_dir, sorted_genes, gene_lengths):
    all_species = set()
    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.startswith("aligned_") and file.endswith(".fasta"):
                for record in SeqIO.parse(os.path.join(root, file), "fasta"):
                    all_species.add(record.id)

    print(f"Found {len(all_species)} species and {len(sorted_genes)} genes.")
    supermatrix = {sp: "" for sp in all_species}

    for gene in sorted_genes:
        print(f"Processing gene: {gene} (Length: {gene_lengths[gene]})")
        gene_file = None
        for root, dirs, files in os.walk(root_dir):
            if f"aligned_{gene}.fasta" in files:
                gene_file = os.path.join(root, f"aligned_{gene}.fasta")

        gene_seqs = {r.id: str(r.seq) for r in SeqIO.parse(gene_file, "fasta")}

        for sp in all_species:
            if sp in gene_seqs:
                supermatrix[sp] += gene_seqs[sp]
            else:
                supermatrix[sp] += "-" * gene_lengths[gene]
                print(f"  [!] Padding missing {gene} for {sp}")

    with open(FASTA_OUT, "w") as out:
        for sp, seq in supermatrix.items():
            out.write(f">{sp}\n{seq}\n")
    print(f"{FASTA_OUT} written (rectangular matrix).")


def create_nexus(fasta_file, partition_file, output_nex):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    n_taxa = len(records)
    n_char = len(records[0].seq)

    with open(output_nex, "w") as f:
        f.write("#NEXUS\n\nBEGIN DATA;\n")
        f.write(f"  DIMENSIONS NTAX={n_taxa} NCHAR={n_char};\n")
        f.write("  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n")
        f.write("  MATRIX\n")
        for record in records:
            f.write(f"  {record.id:<20} {str(record.seq)}\n")
        f.write("  ;\nEND;\n\n")

        f.write("BEGIN SETS;\n")
        with open(partition_file, "r") as p:
            for line in p:
                if "=" in line:
                    parts = line.strip().split(",")[1].split("=")
                    name = parts[0].strip()
                    rng = parts[1].strip()
                    f.write(f"  CHARSET {name} = {rng};\n")
        f.write("END;\n")

    print(f"{output_nex} written for BEAUti.")


if __name__ == "__main__":
    sorted_genes, gene_lengths = generate_partitions(ROOT_DIR)
    concatenate_with_padding(ROOT_DIR, sorted_genes, gene_lengths)
    create_nexus(FASTA_OUT, PARTITIONS_OUT, NEXUS_OUT)
