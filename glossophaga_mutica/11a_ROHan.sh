#!/bin/bash
# Author: Jesús Rocamontes-Morales
# Date:
# Description: Run ROHan (rohmu=8e-4, tstv=2.20) for all Glossophaga mutica
# samples plus the G. leachii outgroup. G. commissarisi was excluded from
# this analysis due to its lower sequencing depth (~10x) relative to the
# other samples (~30x). Run one sample at a time, no loops or
# parallelization, to keep full core allocation per sample.
#
# NOTE: set the paths below to match your own local setup before running.
# ROHAN_BIN: path to your local ROHan installation (https://github.com/grenaud/ROHan),
#            not included in this repository and must be installed separately.
# REFERENCE: path to the reference genome FASTA (autosomes only).
# AUTOSOMES: path to the autosomes list file (one scaffold/chromosome ID per line).
# BAM_DIR:   directory containing the sorted, deduplicated BAM files for each sample.
# OUTPUT_DIR: directory where ROHan output will be written (created if it
#             does not already exist).
# Pipeline:

ROHAN_BIN="/path/to/your/ROHan/bin/rohan"
REFERENCE="/path/to/your/reference/Glossophaga_mutica_autosomes.fna"
AUTOSOMES="/path/to/your/reference/autosomes.txt"
BAM_DIR="/path/to/your/bam/files"
OUTPUT_DIR="/path/to/your/output/directory"

mkdir -p "$OUTPUT_DIR"

# Islas Marias (IslNa)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/IslNa_8e-4_tstv" "$REFERENCE" "$BAM_DIR/GNAY_sorted.bam"

# Concha (ConPa)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/ConPa_8e-4_tstv" "$REFERENCE" "$BAM_DIR/G4CC1_sorted.bam"

# Ixtapa (IxtPa)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/IxtPa_8e-4_tstv" "$REFERENCE" "$BAM_DIR/GIX_sorted.bam"

# Hobonil (HobYu)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/HobYu_8e-4_tstv" "$REFERENCE" "$BAM_DIR/G5YH_sorted.bam"

# Poana (PoaGo)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/PoaGo_8e-4_tstv" "$REFERENCE" "$BAM_DIR/G3TP1_sorted.bam"

# Quilamula (QuiGo)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/QuiGo_8e-4_tstv" "$REFERENCE" "$BAM_DIR/G1QM_sorted.bam"

# Nizanda (NizOa)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/NizOa_8e-4_tstv" "$REFERENCE" "$BAM_DIR/G2ON_sorted.bam"

# G. leachii outgroup (MOR)
$ROHAN_BIN --rohmu 8e-4 --tstv 2.20 -t 12 --auto "$AUTOSOMES" -o "$OUTPUT_DIR/Gleachii_8e-4_tstv" "$REFERENCE" "$BAM_DIR/MOR_sorted.bam"
