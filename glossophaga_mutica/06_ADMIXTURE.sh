#!/bin/bash
# Run ADMIXTURE in unsupervised mode for K = 1 to 6, then collect CV errors.
# Requires ADMIXTURE v1.3.0. Uses 2-fold cross-validation (see Methods),
# chosen to avoid numerical instability with N=7 individuals.
# K = 6 is expected to fail to converge (see Results); it is still run here
# so the failure and CV.txt output are reproducible.

CV_FOLDS=2
THREADS=10

for K in {1..6}; do
  admixture -j${THREADS} --cv=${CV_FOLDS} mutica_final_pruned_popgen.bed $K | tee log${K}.out
done

grep -h "CV error" log*.out > CV.txt
