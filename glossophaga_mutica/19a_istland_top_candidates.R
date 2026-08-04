library(tidyverse)

# Identify candidate selective-sweep windows for Island vs. Mainland:
# top 1% Dxy (empirical), combined with below-average pi in the Island
# population. Requires pixy output from 18_pixy_pipeline.sh.

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_island")

dxy_data <- read_tsv("Island_vs_Mainland_dxy.txt", show_col_types = FALSE)
pi_data  <- read_tsv("Island_vs_Mainland_pi.txt", show_col_types = FALSE)

dxy_clean <- dxy_data %>%
  filter((pop1 == "Island" & pop2 == "Mainland") | (pop1 == "Mainland" & pop2 == "Island")) %>%
  select(chromosome, window_pos_1, window_pos_2, avg_dxy) %>%
  filter(!is.na(avg_dxy))

pi_clean <- pi_data %>%
  filter(pop == "Island") %>%
  select(chromosome, window_pos_1, window_pos_2, avg_pi) %>%
  filter(!is.na(avg_pi))

merged <- inner_join(dxy_clean, pi_clean, by = c("chromosome", "window_pos_1", "window_pos_2"))

dxy_threshold <- quantile(merged$avg_dxy, 0.99, na.rm = TRUE)
pi_mean       <- mean(merged$avg_pi, na.rm = TRUE)

island_candidates <- merged %>%
  filter(avg_dxy >= dxy_threshold) %>%
  filter(avg_pi < pi_mean) %>%
  arrange(desc(avg_dxy)) %>%
  slice_head(n = 200)

print(island_candidates)
write_csv(island_candidates, "Island_Robust_Candidates_Top200.csv")
