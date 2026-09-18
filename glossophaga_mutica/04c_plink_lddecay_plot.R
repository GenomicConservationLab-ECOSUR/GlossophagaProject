# Plot LD decay from PLINK --r2 output (mutica_ld_decay.ld)
# Produces a supplementary figure showing mean pairwise r2 as a function of
# physical distance between SNPs, with reference lines marking the r2=0.2
# threshold and 50kb window used for LD pruning in 04_plink.sh.

library(dplyr)
library(ggplot2)

ld <- read.table("mutica_ld_decay.ld", header = TRUE, stringsAsFactors = FALSE)

ld <- ld %>%
  mutate(dist_bp = abs(BP_B - BP_A),
         dist_kb = dist_bp / 1000)

bin_width_kb <- 5
ld <- ld %>%
  mutate(dist_bin = floor(dist_kb / bin_width_kb) * bin_width_kb)

ld_binned <- ld %>%
  group_by(dist_bin) %>%
  summarise(mean_r2 = mean(R2, na.rm = TRUE),
            n_pairs = n(),
            .groups = "drop")

# Diagnostic: inspect number of pairs per bin
# Sparse bins (few pairs) produce noisy, unreliable mean r2 estimates and
# can create misleading upticks/artifacts at the edges of the distance range.
# Print a summary so you can judge where the data gets too thin to trust.

cat("\n--- Pairs per distance bin (summary) ---\n")
print(summary(ld_binned$n_pairs))

cat("\n--- Bins with the fewest pairs (potential artifacts) ---\n")
print(ld_binned %>% arrange(n_pairs) %>% head(15))

# Save the full per-bin table (with n_pairs) for manual inspection/reporting
write.csv(ld_binned, "mutica_ld_decay_binned.csv", row.names = FALSE)

# Note on bin pair counts
# Each distance bin here contains on the order of 1e5-1e6 raw SNP pairs
# (see diagnostic above), so sparsity of pairs is NOT the issue. With only
# N=7 individuals (14 haplotypes), the number of *independent* haplotype
# comparisons is capped regardless of how many SNP pairs are computed, which
# is why the curve still becomes unstable at longer distances despite large
# pair counts (see truncation below rather than a pair-count filter).

# Truncate to the distance range over which the curve is stable
# Beyond ~450kb, the curve becomes unstable (e.g., an upward tail) despite
# large raw pair counts, consistent with the limited number of independent
# haplotypes available (N=7 individuals) rather than a binning artifact.
# Estimates are reported/plotted only up to max_dist_kb for this reason.

max_dist_kb <- 400
n_dropped <- sum(ld_binned$dist_bin > max_dist_kb)
cat(sprintf("\nTruncating plot to <= %d kb: dropping %d of %d bins beyond this distance.\n",
            max_dist_kb, n_dropped, nrow(ld_binned)))

ld_binned_truncated <- ld_binned %>% filter(dist_bin <= max_dist_kb)

# Plot

p <- ggplot(ld_binned_truncated, aes(x = dist_bin, y = mean_r2)) +
  geom_point(aes(size = n_pairs), alpha = 0.5, color = "grey30") +
  scale_size_continuous(name = "N pairs", range = c(0.5, 3)) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "firebrick") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "darkgreen") +
  labs(x = "Distance between SNP pairs (kb)",
       y = expression(paste("Mean pairwise ", r^2)),
       title = "",) +
  theme_minimal(base_size = 13)

print(p)
