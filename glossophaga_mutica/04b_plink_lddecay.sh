#!/bin/bash
# 04_plink_lddecay.sh (revised)
# Standalone LD decay assessment for G. mutica.

# Requires: mutica_full.bed, mutica_full.bim, mutica_full.fam
# PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)

# Thin to ~1% of SNPs, kept at random and roughly evenly spaced genome-wide.
# Adjust --thin if you want a denser/sparser curve (e.g. 0.05 for more points,
# 0.005 if 0.01 is still too large).
./plink --bfile mutica_full \
  --allow-extra-chr \
  --thin 0.01 \
  --make-bed --out mutica_thinned_for_lddecay

./plink --bfile mutica_thinned_for_lddecay \
  --allow-extra-chr \
  --r2 \
  --ld-window-kb 500 \
  --ld-window 999999 \
  --ld-window-r2 0 \
  --out mutica_ld_decay
