#!/bin/bash
# Align paired-end reads to the G. mutica reference genome with BWA-MEM,
# then sort and index with samtools. Expects input FASTQs named
# <sample>_R1.fq.gz / <sample>_R2.fq.gz in the working directory.
# Requires bwa >= 0.7.17 and samtools >= 1.10.

THREADS=15 # or your preferred number of threads.
REFERENCE="Glossophaga_mutica_autosomes.fna"

# Index the reference only if it hasn't been indexed already
if [ ! -f "${REFERENCE}.bwt" ]; then
  bwa index "$REFERENCE"
fi

> bamlist.txt

# Example of the bamlist.txt content:
#Concha.bam
#Hobonil.bam
#Ixtapa.bam
#Nayarit.bam
#Nizanda.bam
#.. etc


for R1 in *_R1.fq.gz; do
  BASE=$(basename "$R1" _R1.fq.gz)
  R2="${BASE}_R2.fq.gz"
  echo "Processing $BASE..."

  bwa mem -t "$THREADS" "$REFERENCE" "$R1" "$R2" > "${BASE}.sam"
  samtools view -@ "$THREADS" -Sb "${BASE}.sam" > "${BASE}.raw.bam"
  samtools sort -@ "$THREADS" -o "${BASE}.bam" "${BASE}.raw.bam"
  samtools index "${BASE}.bam"
  rm "${BASE}.sam" "${BASE}.raw.bam"

  echo "${BASE}.bam" >> bamlist.txt
done
