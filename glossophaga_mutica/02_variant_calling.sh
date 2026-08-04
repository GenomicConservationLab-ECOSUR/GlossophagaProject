#!/bin/bash
# Per-scaffold variant calling with bcftools mpileup/call, parallelized
# locally across available CPU cores instead of a SLURM job array.
# Requires bcftools >= 1.10.2.
#
# Note on resources: each parallel bcftools process needs enough RAM to
# hold pileup data for all BAMs at that scaffold. Set CORES conservatively
# based on available RAM (e.g. total_RAM_GB / ~10-15 GB per process is a
# reasonable starting estimate for ~9 BAMs at 30-50x) rather than just
# using every core, or the system may swap/OOM under full parallelism.

REFERENCE="Glossophaga_mutica_autosomes.fna"
BAMLIST="bamlist.txt"
SCAFFOLD_LIST="scaffold_list.txt" # Usable for a reference genome with properly defined scaffolds/chromosomes names only. For instance G. mutica has 15 chromosomic scaffolds. Consider your species.
OUTPUT_DIR="split"
CORES=4   # adjust to your machine; see note above

mkdir -p "$OUTPUT_DIR"

call_scaffold() {
  REGION="$1"
  bcftools mpileup -r "$REGION" \
    -a FORMAT/DP,FORMAT/AD \
    -q 20 -Q 20 \
    -Ou -f "$REFERENCE" -b "$BAMLIST" --max-depth 500 | \
  bcftools call -mv -Ob -o "${OUTPUT_DIR}/${REGION}.bcf"
}
export -f call_scaffold
export REFERENCE BAMLIST OUTPUT_DIR

cat "$SCAFFOLD_LIST" | xargs -P "$CORES" -I{} bash -c 'call_scaffold "$@"' _ {}
