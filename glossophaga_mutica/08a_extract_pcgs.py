#!/usr/bin/env python3
# Extract the 13 mitochondrial protein-coding genes (PCGs) from MITOS2-annotated
# GenBank (.gb) files. Handles common gene-name spelling/formatting variants.
# Download all the correct gb files from Genbank to produce our results (read Mitogenomic divergence time estimation section)
# Requires biopython. Run from the directory containing the .gb files.

import os
from Bio import SeqIO

gene_patterns = {
    'ATP6': ['atp6', 'atp 6', 'atp-6', 'apt6', 'apt 6', 'apt-6'],
    'ATP8': ['atp8', 'atp 8', 'atp-8', 'apt8', 'apt 8', 'apt-8'],
    'COX1': ['cox1', 'cox 1', 'cox-1', 'coi', 'co i'],
    'COX2': ['cox2', 'cox 2', 'cox-2', 'coii', 'co ii'],
    'COX3': ['cox3', 'cox 3', 'cox-3', 'coiii', 'co iii'],
    'CYTB': ['cob', 'cytb', 'cyt b', 'cyt-b'],
    'ND1': ['nad1', 'nd1', 'nad 1', 'nd 1'],
    'ND2': ['nad2', 'nd2', 'nad 2', 'nd 2'],
    'ND3': ['nad3', 'nd3', 'nad 3', 'nd 3'],
    'ND4': ['nad4', 'nd4', 'nad 4', 'nd 4'],
    'ND4L': ['nad4l', 'nd4l', 'nad 4l', 'nd 4l'],
    'ND5': ['nad5', 'nd5', 'nad 5', 'nd 5'],
    'ND6': ['nad6', 'nd6', 'nad 6', 'nd 6']
}

output_dir = "extracted_pcgs"
os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir("."):
    if filename.endswith(".gb"):
        print(f"Processing {filename}...")
        try:
            record = SeqIO.read(filename, "genbank")
            sample_name = filename.replace(".gb", "")
            for feature in record.features:
                if feature.type == 'CDS':
                    qualifiers = str(feature.qualifiers).lower()
                    for standard_name, patterns in gene_patterns.items():
                        if any(p in qualifiers for p in patterns):
                            gene_folder = os.path.join(output_dir, standard_name)
                            os.makedirs(gene_folder, exist_ok=True)
                            seq_record = feature.extract(record.seq)
                            out_path = os.path.join(gene_folder, f"{sample_name}_{standard_name}.fasta")
                            with open(out_path, "w") as f:
                                f.write(f">{sample_name}\n{str(seq_record)}\n")
                            break
        except Exception as e:
            print(f"Skipping {filename}: {e}")

print("Extraction complete. Check the 'extracted_pcgs' folder.")
