# PSMC demographic history plot: main estimate + 100 bootstrap replicates
# per individual, scaled to real time (years) and effective population size.
# Scaling follows PSMC/beta-PSMC manual Appendix II:
#   N0 = theta0 / (4*mu) / s
#   T_k = 2 * N0 * t_k * g   (years)
#   N_k = N0 * lambda_k
#
# mu = 7.97e-9 (Bergeron et al. 2023), g = 1 year (Zorte'a 2003), s = 100
# (default fq2psmcfa bin size; confirm no -s flag was used upstream).

setwd("/run/media/rocamontes/EcosurMutica/PSMC/graphic")


library(ggplot2)
library(scales)


MU <- 7.97e-9
G <- 1
S <- 100
# Number of earliest and latest time bins to drop per run (standard PSMC
# practice -- these bins are the least constrained by data and dominate
# bootstrap noise at the extremes of the trajectory).
TRIM_N <- 2

# PSMC has well-documented poor resolution in the most recent time bins,
# where too few recombination events have accumulated to constrain the
# coalescent model. However, a blanket exclusion of this zone would also
# remove genuine signal for samples where it IS informative (e.g. IslNa's
# reduced recent Ne, independently corroborated by its low genome-wide
# theta and elevated ROH). Rather than deleting this data, it is shaded
# on the plot as a visual low-confidence zone; all data is retained in
# main_df/boot_df and in the summary statistics below. Treat any value
# from within this window with appropriate caution, cross-checked against
# independent evidence (theta, ROH, etc.) before interpretation.
LOW_CONFIDENCE_MAX_YEARS <- 10000

parse_psmc_file <- function(filepath, trim_n = TRIM_N) {
  lines <- readLines(filepath)
  run_starts <- grep("^MM\\tVersion:", lines)
  run_ends <- c(run_starts[-1] - 1, length(lines))

  results <- list()
  for (i in seq_along(run_starts)) {
    run_lines <- lines[run_starts[i]:run_ends[i]]
    rd_starts <- grep("^RD\\t", run_lines)
    if (length(rd_starts) == 0) next

    last_rd_start <- rd_starts[length(rd_starts)]
    final_block <- run_lines[last_rd_start:length(run_lines)]

    tr_line <- grep("^TR\\t", final_block, value = TRUE)
    theta0 <- as.numeric(strsplit(tr_line[1], "\t")[[1]][2])

    rs_lines <- grep("^RS\\t", final_block, value = TRUE)
    rs_split <- strsplit(rs_lines, "\t")
    t_k <- sapply(rs_split, function(x) as.numeric(x[3]))
    lambda_k <- sapply(rs_split, function(x) as.numeric(x[4]))

    # Trim the first and last trim_n bins
    n <- length(t_k)
    if (n > 2 * trim_n) {
      keep <- (trim_n + 1):(n - trim_n)
      t_k <- t_k[keep]
      lambda_k <- lambda_k[keep]
    }

    results[[i]] <- list(theta0 = theta0, t_k = t_k, lambda_k = lambda_k)
  }
  return(results)
}

scale_run <- function(run, mu = MU, g = G, s = S) {
  N0 <- run$theta0 / (4 * mu) / s
  data.frame(
    years = 2 * N0 * run$t_k * g,
    Ne = N0 * run$lambda_k
  )
}

# Sample files -> publication codes
sample_files <- c(
  "concha_bootstrap.psmc"    = "ConPa",
  "hobonil_bootstrap.psmc"   = "HobYu",
  "ixtapa_bootstrap.psmc"    = "IxtPa",
  "nayarit_bootstrap.psmc"   = "IslNa",
  "nizanda_bootstrap.psmc"   = "NizOa",
  "poana_bootstrap.psmc"     = "PoaGo",
  "quilamula_bootstrap.psmc" = "QuiGo"
)

main_df <- data.frame()
boot_df <- data.frame()

for (file in names(sample_files)) {
  code <- sample_files[file]
  runs <- parse_psmc_file(file)

  # First run in the concatenated file = main estimate; remaining = bootstraps
  main_scaled <- scale_run(runs[[1]])
  main_scaled$Sample <- code
  main_df <- rbind(main_df, main_scaled)

  for (b in 2:length(runs)) {
    boot_scaled <- scale_run(runs[[b]])
    boot_scaled$Sample <- code
    boot_scaled$Rep <- b - 1
    boot_df <- rbind(boot_df, boot_scaled)
  }
}

# 7 distinct, colorblind-safe colors (Okabe-Ito, extended)
sample_colors <- c(
  "ConPa" = "#0072B2",
  "IxtPa" = "#56B4E9",
  "PoaGo" = "#D55E00",
  "QuiGo" = "#E69F00",
  "HobYu" = "#009E73",
  "NizOa" = "#F0E442",
  "IslNa" = "#CC79A7"
)

p <- ggplot() +
  annotate("rect",
           xmin = min(c(main_df$years, boot_df$years), na.rm = TRUE) * 0.9,
           xmax = LOW_CONFIDENCE_MAX_YEARS, ymin = 1e2, ymax = 1e9,
           fill = "grey50", alpha = 0.15) +
  geom_step(data = boot_df,
            aes(x = years, y = Ne, group = interaction(Sample, Rep), color = Sample),
            alpha = 0.05, linewidth = 0.3) +
  geom_step(data = main_df,
            aes(x = years, y = Ne, color = Sample),
            linewidth = 1) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  coord_cartesian(ylim = c(1e2, 1e9)) +
  scale_color_manual(values = sample_colors) +
  labs(x = "Years ago", y = expression("Effective population size (" * N[e] * ")"),
       color = "Sample") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")

p

ggsave("psmc_demographic_history.pdf", p, width = 10, height = 7)

# --- Summary statistics per sample, for reproducible reporting in Results ---
# Reports min/max Ne (and the corresponding time) from the MAIN estimate only
# (not bootstrap replicates), across the full trajectory and within the
# 10,000-50,000 years ago window. All data is included; each extreme is
# flagged if it falls within the low-confidence zone (< LOW_CONFIDENCE_MAX_YEARS)
# so that values needing a reliability caveat are identified explicitly rather
# than silently included or excluded.

summarize_sample <- function(df, sample_code) {
  d <- df[df$Sample == sample_code, ]
  d <- d[order(d$years), ]

  overall_min <- d[which.min(d$Ne), ]
  overall_max <- d[which.max(d$Ne), ]

  window <- d[d$years >= 10000 & d$years <= 50000, ]
  if (nrow(window) > 0) {
    window_min <- window[which.min(window$Ne), ]
    window_max <- window[which.max(window$Ne), ]
  } else {
    window_min <- data.frame(years = NA, Ne = NA)
    window_max <- data.frame(years = NA, Ne = NA)
  }

  data.frame(
    Sample = sample_code,
    Ne_min_overall = overall_min$Ne,
    years_at_Ne_min = overall_min$years,
    Ne_min_low_confidence = overall_min$years < LOW_CONFIDENCE_MAX_YEARS,
    Ne_max_overall = overall_max$Ne,
    years_at_Ne_max = overall_max$years,
    Ne_max_low_confidence = overall_max$years < LOW_CONFIDENCE_MAX_YEARS,
    Ne_min_10k_50k = window_min$Ne,
    years_at_min_10k_50k = window_min$years,
    Ne_max_10k_50k = window_max$Ne,
    years_at_max_10k_50k = window_max$years
  )
}

summary_table <- do.call(rbind, lapply(names(sample_files), function(f) {
  summarize_sample(main_df, sample_files[f])
}))

print(summary_table)
write.csv(summary_table, "psmc_summary_stats.csv", row.names = FALSE)
