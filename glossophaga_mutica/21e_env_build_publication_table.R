library(tidyverse)

# Build publication-ready tables for each environmental cluster: genuine,
# BLAST-confirmed protein-coding candidates only. TE/repeat-element hits,
# below-significance hits, and non-coding windows are excluded entirely
# (matching the island analysis). Full detail for all windows remains in
# each Validated_<cluster>_Annotations_Corrected.csv for supplementary
# reporting/reproducibility.
#
# Clusters with no gene-overlapping candidates (e.g. Poana, Quilamula --
# confirmed genuine gene deserts, not a pipeline error) are skipped
# automatically if their Corrected file does not exist.

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_environmental")

CLUSTERS <- c("Concha", "Hobo_Ixta", "Nizanda", "Poana", "Quilamula")

chr_map <- c(
  "CM077303.1" = "1",  "CM077304.1" = "2",  "CM077305.1" = "3",
  "CM077306.1" = "4",  "CM077307.1" = "5",  "CM077308.1" = "6",
  "CM077309.1" = "7",  "CM077310.1" = "8",  "CM077311.1" = "9",
  "CM077312.1" = "10", "CM077313.1" = "11", "CM077314.1" = "12",
  "CM077315.1" = "13", "CM077316.1" = "14", "CM077317.1" = "15"
)

for (cluster in CLUSTERS) {

  input_file <- paste0("Validated_", cluster, "_Annotations_Corrected.csv")
  output_file <- paste0("Table_", cluster, "_Publication.csv")

  if (!file.exists(input_file)) {
    cat(sprintf("Skipping %s -- no corrected annotation file (no gene-overlapping candidates).\n", cluster))
    next
  }

  df <- read_csv(input_file, show_col_types = FALSE)

  n_total <- nrow(df)
  n_te <- sum(df$Is_TE_Hit, na.rm = TRUE)
  n_no_hit <- sum(df$No_Blast_Hit, na.rm = TRUE)
  n_below_threshold <- sum(!df$Is_Significant & !df$No_Blast_Hit, na.rm = TRUE)

  publication_table <- df %>%
    filter(Is_Significant, !Is_TE_Hit, !No_Blast_Hit) %>%
    mutate(
      Chromosome = paste0("Chr", recode(chromosome, !!!chr_map)),
      `Position (Mb)` = round(window_pos_1 / 1e6, 3)
    ) %>%
    select(
      Chromosome,
      `Position (Mb)`,
      Gene = Gene_Symbol_Final,
      Dxy = avg_dxy,
      Pi = avg_pi,
      Function = Function_Note
    ) %>%
    arrange(Chromosome, `Position (Mb)`)

  write_excel_csv(publication_table, output_file)

  cat(sprintf("\n%s: %d protein-coding candidates (of %d total candidate windows)\n",
              cluster, nrow(publication_table), n_total))
  cat(sprintf("  Excluded: %d TE/repeat-element hits, %d below significance threshold, %d non-coding (no BLAST hit)\n",
              n_te, n_below_threshold, n_no_hit))
  cat(sprintf("  Saved: %s\n", output_file))
}

cat("\nAll clusters processed.\n")
