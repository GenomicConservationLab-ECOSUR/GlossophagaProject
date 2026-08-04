#!/bin/bash
# Run beta-PSMC (main estimate + 100 bootstrap replicates) for one individual.
# Requires bcftools, samtools, and beta-psmc (github.com/emiliocarrera1/beta-psmc
# or equivalent). Parallelized locally across available CPU cores.
#
# CHANGE THIS: set SAMPLE and BAM to match the individual you are running.
SAMPLE="concha"
BAM="Concha.bam"

REFERENCE="Glossophaga_mutica_autosomes.fna"
CORES=8   # adjust to your machine; each beta-psmc bootstrap process is single-threaded

# Generate diploid consensus and PSMC input
bcftools mpileup -C50 -f "$REFERENCE" "$BAM" | \
    bcftools call -c -V indels | \
    vcfutils.pl vcf2fq -d 13 -D 76 | gzip > "${SAMPLE}.fq.gz"

fq2psmcfa -q20 "${SAMPLE}.fq.gz" > "${SAMPLE}.psmcfa"

# Main (point-estimate) PSMC run, no bootstrap
beta-psmc -N25 -t15 -r5 -p "4+25*2+4+6" -B5 -o "${SAMPLE}.psmc" "${SAMPLE}.psmcfa"

# Split genome for bootstrap resampling
splitfa "${SAMPLE}.psmcfa" > "split_${SAMPLE}.psmcfa"

# 100 bootstrap replicates, run in parallel, resumable (skips any round
# already completed, e.g. after a job restart)
for i in $(seq 1 100); do
    if [ ! -f "round-${i}.psmc" ]; then
        beta-psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -B5 -o "round-${i}.psmc" "split_${SAMPLE}.psmcfa" &
        if (( $(jobs -r | wc -l) >= CORES )); then
            wait -n
        fi
    fi
done
wait

# Concatenate main estimate + all bootstrap replicates into one file,
# ready for plotting (see 16_psmc_plot.R)
cat "${SAMPLE}.psmc" round-*.psmc > "${SAMPLE}_bootstrap.psmc"
