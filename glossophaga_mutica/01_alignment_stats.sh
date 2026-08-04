#!/bin/bash
# Per-sample raw read pairs, mapping rate, and mean autosomal coverage
# Requires samtools >= 1.10 (for samtools coverage). Run in the directory
# containing bamlist.txt and the indexed BAMs produced by the alignment step.

BAMLIST="bamlist.txt"
OUT="alignment_stats.tsv"

# Create output with header only if it doesn't exist yet (resumable)
if [ ! -f "$OUT" ]; then
  echo -e "sample\traw_read_pairs\tmapped_reads\tmapping_rate_pct\tmean_coverage_x" > "$OUT"
fi

while read -r BAM; do
  SAMPLE=$(basename "$BAM" .bam)

  # Skip samples already processed
  if grep -q "^${SAMPLE}	" "$OUT"; then
    echo "Skipping ${SAMPLE}, already done."
    continue
  fi

  FLAGSTAT=$(samtools flagstat "$BAM")
  TOTAL=$(echo "$FLAGSTAT" | awk 'NR==1{print $1}')
  MAPPED_LINE=$(echo "$FLAGSTAT" | grep " mapped (" | head -1)
  MAPPED=$(echo "$MAPPED_LINE" | awk '{print $1}')
  RATE=$(echo "$MAPPED_LINE" | grep -oP '\(\K[0-9.]+(?=%)')
  RAW_PAIRS=$((TOTAL / 2))

  DEPTH=$(samtools coverage "$BAM" | \
    awk 'NR>1{sum+=$7*($3-$2+1); len+=($3-$2+1)} END{printf "%.1f", sum/len}')

  echo -e "${SAMPLE}\t${RAW_PAIRS}\t${MAPPED}\t${RATE}\t${DEPTH}" >> "$OUT"
  echo "Done: ${SAMPLE}"
done < "$BAMLIST"
