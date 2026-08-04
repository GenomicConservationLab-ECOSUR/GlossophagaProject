#!/bin/bash
# Nuclear phylogenetic inference: VCF filtering, SVDquartets (PAUP*), and
# IQ-TREE maximum likelihood analysis.
# Requires: vcftools, vcf2phylip.py (Ortiz 2019), PAUP* v4.0a build 169,
# IQ-TREE v3.0.1.

BASE="your directory"
BCF="directory/all.bcf" # result of variant calling script (02_variant_calling.sh)
SVDQ_DIR="$BASE/SVDquartets"
IQTREE_DIR="$BASE/IQTree"
PAUP="$BASE/paup4a169_ubuntu64"
VCF2PHY="$BASE/vcf2phylip.py"

# --- Environment note ---
# The libcurl symlink fix and LD_LIBRARY_PATH export below are specific to
# a conda environment named "svdq" on the original machine, and are only
# needed if PAUP* fails to start due to a missing/mismatched libcurl. Adjust
# or remove for your own environment; not required on all systems.
LIBPATH="/home/rocamontes/miniconda3/envs/svdq/lib"
ln -sf "$(conda info --base)/envs/svdq/lib/libcurl.so.4" \
       "$(conda info --base)/envs/svdq/lib/libcurl-gnutls.so.4"
export LD_LIBRARY_PATH=$LIBPATH

# Step 1: Filter VCF (biallelic SNPs, MAF >= 0.05, <=20% missing per site)
cd "$BASE"
vcftools --gzvcf "$BCF" \
  --max-missing 0.8 \
  --maf 0.05 \
  --recode \
  --out mutica_phylo_filtered

# Step 2: Convert to NEXUS/PHYLIP for SVDquartets
cd "$SVDQ_DIR"
python3 "$VCF2PHY" \
  -i "$BASE/mutica_phylo_filtered.recode.vcf" \
  -n \
  -o Ocosingo.bam \
  -m 7 \
  --output-prefix mutica_svdq

# Step 3: SVDquartets - exhaustive quartet evaluation, 1000 bootstrap replicates
head -n 18 mutica_svdq.min7.nexus > mutica_svdq_run.nexus
cat >> mutica_svdq_run.nexus << 'EOF'
BEGIN PAUP;
    log file=mutica_svdq_log.txt start replace;
    svdquartets evalquartets=all bootstrap=standard nreps=1000 treefile=mutica_svdq_tree.tre replace;
    quit;
END;
EOF
LD_LIBRARY_PATH=$LIBPATH "$PAUP" mutica_svdq_run.nexus

# Step 4: 50% majority-rule consensus tree from bootstrap replicates
cat > mutica_consensus.nex << 'EOF'
#NEXUS
BEGIN PAUP;
    log file=mutica_consensus_log.txt start replace;
    set maxtrees=1000;
    gettrees file=mutica_svdq_tree.tre;
    contree all / strict=no majrule=yes percent=50 treefile=mutica_consensus_bootstrap.tre replace;
    quit;
END;
EOF
LD_LIBRARY_PATH=$LIBPATH "$PAUP" mutica_consensus.nex

# Step 5: IQ-TREE ML analysis (GTR+ASC, 1000 ultrafast bootstrap replicates)
cd "$IQTREE_DIR"
python3 "$VCF2PHY" \
  -i "$BASE/mutica_phylo_filtered.recode.vcf" \
  -o Ocosingo.bam \
  -m 7 \
  --output-prefix mutica_iqtree

iqtree \
  -s mutica_iqtree_output.varsites.phy \
  -st DNA \
  -m GTR+ASC \
  -B 1000 \
  -T AUTO \
  --boot-trees \
  --prefix mutica_iqtree_output3
