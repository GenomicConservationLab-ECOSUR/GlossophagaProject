#!/bin/bash
# Run ADMIXTURE in unsupervised mode for K = 2 to 6 with replicate runs, then keep the
# highest-likelihood run per K. Requires ADMIXTURE v1.3.0.
#
# - 2-fold cross-validation (--cv=2), given N = 7 individuals (see Methods).
# - ADMIXTURE can converge to different local optima, so K = 2-5 are run from 20 random
#   seeds each; the run with the highest log-likelihood is reported.
# - Runs whose log-likelihood diverges to nan are killed after TIMEOUT and marked FAILED.
# - K = 6 is expected to fail to converge (see Results); it is attempted once so the
#   failure is reproducible.
# - The script is resumable: finished and FAILED runs are skipped.

BED=mutica_final_pruned_popgen
CV_FOLDS=2
THREADS=10
TIMEOUT=1h

declare -A NSEEDS=( [2]=20 [3]=20 [4]=20 [5]=20 [6]=1 )

for K in 2 3 4 5 6; do
  for S in $(seq 1 ${NSEEDS[$K]}); do
    DIR=seeds/K${K}_s${S}
    [ -f $DIR/FAILED ] && continue
    [ -f $DIR/log.out ] && grep -q "^Loglikelihood:" $DIR/log.out && continue
    mkdir -p $DIR
    (cd $DIR && ln -sf ../../${BED}.{bed,bim,fam} . \
      && timeout ${TIMEOUT} admixture -j${THREADS} --cv=${CV_FOLDS} -s $S ${BED}.bed $K > log.out; \
      grep -q "^Loglikelihood:" log.out || touch FAILED)
  done
done

# Summary of all finished replicates: K, seed, log-likelihood, CV error
# (best log-likelihood first within each K; CV = "-nan" where the CV step failed)
for d in seeds/K*_s*; do
  b=$(basename $d); K=${b%%_*}; K=${K#K}; s=${b##*_s}
  ll=$(grep '^Loglikelihood:' $d/log.out 2>/dev/null | awk '{print $2}')
  cv=$(grep 'CV error' $d/log.out 2>/dev/null | sed 's/.*: //')
  [ -n "$ll" ] && echo -e "$K\t$s\t$ll\t$cv"
done | sort -k1,1n -k3,3gr > seed_summary.tsv

# Copy the highest-likelihood .Q file for each converged K to best/ (used by the plotting script)
mkdir -p best
for K in 2 3 4 5; do
  S=$(awk -v K=$K '$1==K {print $2; exit}' seed_summary.tsv)
  cp seeds/K${K}_s${S}/${BED}.${K}.Q best/
done
