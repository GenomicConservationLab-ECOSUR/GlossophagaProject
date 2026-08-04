library(tidyverse)

# Identify candidate selective-sweep windows for each environmental cluster
# (vs. all other clusters): top 1% Dxy (empirical), combined with
# below-average pi in the focal cluster. Requires pixy output from
# 18_pixy_pipeline.sh.

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_environmental")

environmental_clusters <- c("Concha", "Hobo_Ixta", "Nizanda", "Poana", "Quilamula")

for (cluster in environmental_clusters) {

  dxy_raw <- read_tsv(paste0("Cluster_", cluster, "_vs_others_dxy.txt"), show_col_types = FALSE)
  pi_raw  <- read_tsv(paste0("Cluster_", cluster, "_vs_others_pi.txt"), show_col_types = FALSE)

  dxy_clean <- dxy_raw %>%
    filter(!is.na(avg_dxy)) %>%
    select(chromosome, window_pos_1, window_pos_2, avg_dxy)

  pi_clean <- pi_raw %>%
    filter(!is.na(avg_pi)) %>%
    select(chromosome, window_pos_1, window_pos_2, avg_pi)

  merged <- inner_join(dxy_clean, pi_clean, by = c("chromosome", "window_pos_1", "window_pos_2"))

  dxy_threshold <- quantile(merged$avg_dxy, 0.99, na.rm = TRUE)
  pi_mean       <- mean(merged$avg_pi, na.rm = TRUE)

  cluster_candidates <- merged %>%
    filter(avg_dxy >= dxy_threshold) %>%
    filter(avg_pi < pi_mean) %>%
    arrange(desc(avg_dxy)) %>%
    slice_head(n = 200)

  write_csv(cluster_candidates, paste0("Cluster_", cluster, "_Robust_Candidates_Top200.csv"))
}
