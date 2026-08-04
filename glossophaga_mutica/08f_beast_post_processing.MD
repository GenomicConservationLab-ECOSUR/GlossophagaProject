# 08f — BEAST post-processing (LogCombiner / TreeAnnotator)

Steps 08a-08d prepare the input; 08e (Gmutica.xml) is the BEAST run
configuration itself, executed via BEAST v2.7.5. This step covers what
comes after the two independent MCMC runs finish, done via the BEAST
GUI tools rather than a script.

## 1. Combine independent runs (LogCombiner)

Combine the log and tree files from the two independent MCMC runs
(Dated_1 / Dated_2) into a single combined posterior, discarding an
appropriate burn-in from each.

GUI: LogCombiner -> add both .log files -> set burn-in -> output
Combined.log. Repeat separately for the .trees files -> Combined.trees.

Command-line equivalent:
```
logcombiner -log Dated_1.log -log Dated_2.log -burnin 10 -o Combined.log
logcombiner -log Dated_1.trees -log Dated_2.trees -burnin 10 -o Combined.trees
```
(Adjust -burnin to the percentage/generation count appropriate for your
runs; confirm convergence and ESS > 200 in Tracer before combining.)

## 2. Build the maximum clade credibility tree (TreeAnnotator)

From Combined.trees, produce the single summary tree (Mutica.tree) with
mean node heights and 95% HPD intervals.

GUI: TreeAnnotator -> input Combined.trees -> burn-in already removed at
the LogCombiner step, so set to 0 here -> Node heights: Mean heights ->
output Mutica.tree.

Command-line equivalent:
```
treeannotator -burnin 0 -heights mean Combined.trees Mutica.tree
```

## 3. Check convergence (Tracer)

Before treating any date as final, open Combined.log in Tracer and
confirm ESS > 200 for posterior, tree height, and clock rate.
