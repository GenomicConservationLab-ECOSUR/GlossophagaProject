README
======

Scripts for the Glossophaga mutica whole-genome analysis pipeline.
Numbering reflects the order of analysis, not strict execution dependency
across sections (e.g. 20_env_pca.R must be run before 18_pixy_pipeline.sh,
since its output feeds the pixy population map).

00-04   Read alignment, alignment stats, variant calling, and filtering
        to the population-genetics SNP dataset (PLINK-ready).

05-07   Population structure: PCA and ADMIXTURE (K=1-6).

08a-08f Mitogenomic phylogenetics: coding gene extraction, alignment,
        supermatrix assembly, partition search (IQ-TREE), and BEAST
        divergence dating (post-processing notes in 08f).

09-10   Redundancy analysis (RDA) of genetic structure vs. geography
        and environment.

11a-11b Genome-wide heterozygosity and runs of homozygosity (ROHan),
        plus the per-chromosome heterozygosity summary figure.

12      Nuclear phylogenetics: maximum-likelihood (IQ-TREE) and
        coalescent-based (SVDquartets) inference.

13-14   TreeMix: population graph and migration edge estimation.

15-16   Demographic history via beta-PSMC.

17      Species distribution modeling (present, LIG, LGM).

18      Pixy pipeline: sliding-window Dxy/pi for both population designs
        (Island vs. Mainland; environmental clusters). Requires the
        population map produced by 20_env_pca.R.

19a-19g Island vs. Mainland selection scan: candidate window selection,
        gene annotation, BLASTX identification, significance filtering,
        publication table, and Manhattan plot.

20      Environmental cluster assignment via PCA/UPGMA (mainland
        samples only; Islas Marias analyzed separately as Island vs.
        Mainland, see 19a-19g).

21a-21f Environmental cluster selection scan: same pipeline as
        19a-19g, applied to the five environmentally-defined mainland
        clusters.

Each script's header comment documents required input paths and any
tool-specific setup (e.g. local ROHan installation, UniProt database).
