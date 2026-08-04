library(tidyverse)

# Build the publication-ready main-text table: genuine, BLAST-confirmed
# protein-coding candidates only. TE hits, below-significance hits, and
# non-coding (lncRNA) windows are excluded here entirely -- not just
# de-emphasized -- since they are not, by construction, host-gene
# divergence candidates. Full detail for all windows (including excluded
# categories) remains available in Validated_Island_Annotations_Corrected.csv
# for supplementary reporting/reproducibility.

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_island")

df <- read_csv("Validated_Island_Annotations_Corrected.csv", show_col_types = FALSE)

n_total <- nrow(df)
n_te <- sum(df$Is_TE_Hit, na.rm = TRUE)
n_no_hit <- sum(df$No_Blast_Hit, na.rm = TRUE)
n_below_threshold <- sum(!df$Is_Significant & !df$No_Blast_Hit, na.rm = TRUE)

publication_table <- df %>%
  filter(Is_Significant, !Is_TE_Hit, !No_Blast_Hit) %>%
  select(
    Chromosome,
    `Position (Mb)` = Pos_Mb,
    Gene = Gene_Symbol_Final,
    Dxy = Max_Dxy,
    Pi = Avg_Pi,
    Function = Function_Note
  ) %>%
  mutate(CHR_NUM = as.numeric(str_remove(Chromosome, "Chr"))) %>%
  arrange(CHR_NUM, `Position (Mb)`) %>%
  select(-CHR_NUM)

# write_excel_csv() adds a UTF-8 byte-order mark so Excel/LibreCalc
# correctly auto-detect encoding on open (real UniProt Function text can
# contain non-ASCII characters, e.g. Greek letters in compound names).
write_excel_csv(publication_table, "Table_Island_Publication.csv")

message(sprintf("Publication table: %d protein-coding candidates (of %d total candidate windows)",
                nrow(publication_table), n_total))
message(sprintf("Excluded: %d LINE-1/retrotransposon hits, %d below significance threshold, %d non-coding (no BLAST hit)",
                n_te, n_below_threshold, n_no_hit))
message("Full detail for all windows (including excluded categories) remains in Validated_Island_Annotations_Corrected.csv")
