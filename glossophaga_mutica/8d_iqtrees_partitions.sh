#!/bin/bash
# ModelFinder partition/model search on the mitogenomic supermatrix.
# Produces partitions.txt.best_scheme.nex, used to inform substitution
# model choices for the BEAST XML. Requires IQ-TREE.

iqtree -s supermatrix.fasta -q partitions.txt -m MFP -bb 1000 -alrt 1000
