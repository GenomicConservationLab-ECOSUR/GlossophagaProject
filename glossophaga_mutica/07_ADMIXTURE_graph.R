setwd("your_directory")

library(ggplot2)
library(dplyr)
library(tidyr)

# Map raw BAM-derived sample names to publication codes
code_map <- c(
  "Concha"    = "ConPa",
  "Ixtapa"    = "IxtPa",
  "Poana"     = "PoaGo",
  "Quilamula" = "QuiGo",
  "Hobonil"   = "HobYu",
  "Nizanda"   = "NizOa",
  "Nayarit"   = "IslNa"
)

# Geographic order for plotting
geo_order <- c("IslNa", "PoaGo", "QuiGo", "ConPa", "IxtPa", "NizOa", "HobYu")

fam <- read.table("mutica_final_pruned_popgen.fam")

# Consistent colours across K: each cluster takes the colour of the first reference
# individual (in the order below) whose dominant cluster it is. Colours are from Set3.
reference_cols <- c(IslNa = "#FB8072",   # red
                    ConPa = "#F5E663",   # yellow (central group)
                    HobYu = "#80B1D3",   # blue
                    NizOa = "#8DD3C7",   # teal
                    QuiGo = "#FDB462")   # orange
fallback_col <- "grey70"                 # a cluster none of the references dominate
colour_levels <- c(unname(reference_cols), fallback_col)   # fixed stacking order in every panel

save_single_panels <- TRUE   # also write one PDF per K (same size), as before
fig_w <- 7                   # single-panel size (inches)
fig_h <- 3

# --- Read every K, assign colours, and collect the long-format data ---
all_long <- list()

for (K_value in 2:5) {
  
  Q <- read.table(paste0("mutica_final_pruned_popgen.", K_value, ".Q"))
  
  df <- cbind(fam, Q)
  colnames(df) <- c("FamilyID", "IndID", "PatID", "MatID", "Sex", "Phenotype",
                    paste0("Cluster", 1:K_value))
  
  df$IndID <- gsub("\\.bam$", "", df$IndID)
  df$IndID <- code_map[df$IndID]
  df$IndID <- factor(df$IndID, levels = geo_order)
  
  Qmat <- as.matrix(df[, paste0("Cluster", 1:K_value)])
  rownames(Qmat) <- as.character(df$IndID)
  
  cluster_cols <- rep(NA_character_, K_value)
  for (ind in names(reference_cols)) {
    j <- which.max(Qmat[ind, ])
    if (is.na(cluster_cols[j])) cluster_cols[j] <- reference_cols[[ind]]
  }
  cluster_cols[is.na(cluster_cols)] <- fallback_col
  names(cluster_cols) <- paste0("Cluster", 1:K_value)
  
  df_long <- df %>%
    pivot_longer(cols = starts_with("Cluster"),
                 names_to = "Cluster",
                 values_to = "Ancestry")
  df_long$Colour <- factor(unname(cluster_cols[df_long$Cluster]), levels = colour_levels)
  df_long$K <- K_value
  all_long[[as.character(K_value)]] <- df_long
  
  write.csv(df, paste0("admixture_K", K_value, "_ancestry_proportions.csv"), row.names = FALSE)
  
  if (save_single_panels) {
    single_plot <- ggplot(df_long, aes(x = IndID, y = Ancestry, fill = Colour)) +
      geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.4) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "italic"),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(size = 16, face = "bold"),
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      labs(title = paste0("K = ", K_value), x = NULL, y = "Ancestry proportion") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_identity()
    ggsave(paste0("admixture_K", K_value, ".pdf"), single_plot, width = fig_w, height = fig_h)
  }
}

admix_all <- bind_rows(all_long)

# --- Multipanel figure: one row per K, shared x axis and y title ---
multi_plot <- ggplot(admix_all, aes(x = IndID, y = Ancestry, fill = Colour)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.4) +
  facet_grid(K ~ ., labeller = label_bquote(rows = italic(K) ~ "=" ~ .(K))) +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 12, face = "italic"),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 16),
        strip.text.y = element_text(angle = 0, size = 16),
        strip.background = element_blank(),
        panel.spacing.y = unit(1.2, "lines"),
        plot.margin = margin(t = 10, r = 5, b = 5, l = 5),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

multi_plot
ggsave("admixture_K2-5_multipanel.pdf", multi_plot, width = 7, height = 9.5)
ggsave("admixture_K2-5_multipanel.png", multi_plot, width = 7, height = 9.5, dpi = 300)

# --- Cross-validation error across replicate runs ---
# seed_summary.tsv columns: K, seed, log-likelihood, CV error ("-nan" -> NA)
rep_runs <- read.table("seed_summary.tsv", sep = "\t", header = FALSE,
                       col.names = c("K", "seed", "loglik", "CV"),
                       stringsAsFactors = FALSE, fill = TRUE, na.strings = c("-nan", "NA", ""))
rep_runs$CV <- suppressWarnings(as.numeric(rep_runs$CV))
rep_runs <- rep_runs[rep_runs$K %in% 2:5, ]

best_runs <- rep_runs %>% group_by(K) %>% slice_max(loglik, n = 1) %>% ungroup()
cv_valid  <- rep_runs[is.finite(rep_runs$CV), ]

cv_plot <- ggplot(cv_valid, aes(x = K, y = CV)) +
  geom_jitter(width = 0.08, height = 0, size = 2, alpha = 0.5, color = "grey40") +
  stat_summary(fun = median, geom = "line", linewidth = 0.6) +
  geom_point(data = best_runs[is.finite(best_runs$CV), ], size = 4, color = "firebrick") +
  theme_minimal(base_size = 14) +
  labs(x = "K", y = "Cross-Validation Error (2-fold)") +
  scale_x_continuous(breaks = 2:5)

cv_plot
ggsave("admixture_cv_error.pdf", cv_plot, width = 7, height = 5)

cat("Best-likelihood run per K:\n"); print(best_runs)
cat("Replicates with finite CV per K:\n"); print(table(cv_valid$K))
cat("Replicates with CV = nan per K:\n"); print(table(rep_runs$K[!is.finite(rep_runs$CV)]))