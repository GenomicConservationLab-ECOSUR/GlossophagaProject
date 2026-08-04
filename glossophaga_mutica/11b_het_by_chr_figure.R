library(data.table)
library(dplyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)

# Config ------------------------------------------------------------------
# NOTE: set base_dir below to your own local ROHan output directory before
# running. This should be the parent folder containing your Results8e-4/
# subfolder (produced by run_rohan_all_samples.sh) and your autosomes list.
base_dir <- "/path/to/your/rohan/output/directory"
results_dir <- file.path(base_dir, "Results8e-4")
autosomes_file <- file.path(base_dir, "autosomes.txt")

# --- Task 1: Genome-Wide Heterozygosity Table ---
summary_files <- list.files(path = results_dir,
                            pattern = "\\.summary\\.txt$",
                            recursive = TRUE,
                            full.names = TRUE)

# FIX: sample name parsed from the FILENAME prefix, not the parent
# directory -- our files are flat (Results8e-4/<Sample>_8e-4_tstv.summary.txt),
# not nested in per-sample subfolders. This pattern must match the -o
# naming convention used in run_rohan_all_samples.sh.
names(summary_files) <- sub("_8e-4_tstv\\.summary\\.txt$", "", basename(summary_files))

parse_rohan_summary <- function(file_path) {
  lines <- readLines(file_path)
  het_line <- grep("Genome-wide theta outside ROH:", lines, value = TRUE)
  if (length(het_line) == 0) {
    warning(paste("Could not find 'Genome-wide theta outside ROH' in:", file_path))
    return(data.table(HetEst = NA_real_))
  }
  het_value_str <- sub(".*:\\s*(\\d+\\.\\d+).*", "\\1", het_line[1])
  return(data.table(HetEst = as.numeric(het_value_str)))
}

genome_wide_het <- map_dfr(summary_files, parse_rohan_summary, .id = "Sample")

output_csv <- file.path(base_dir, "all_samples_genome_wide_het.csv")
fwrite(genome_wide_het, output_csv)

# --- Task 2/3: Per-chromosome heterozygosity ---
hest_files <- list.files(path = results_dir,
                         pattern = "\\.hEst\\.gz$",
                         recursive = TRUE,
                         full.names = TRUE)
# Same filename-prefix fix as above
names(hest_files) <- sub("_8e-4_tstv\\.hEst\\.gz$", "", basename(hest_files))

all_het_data <- map_dfr(hest_files, ~ fread(.), .id = "Sample")

all_het_data <- all_het_data %>%
  rename(Chr = `#CHROM`, Pos = BEGIN, End = END, Het = h)

autosomes <- fread(autosomes_file, header = FALSE)$V1

plot_data <- all_het_data %>%
  filter(Chr %in% autosomes)

chr_layout <- plot_data %>%
  group_by(Chr) %>%
  summarize(max_pos = max(End)) %>%
  ungroup() %>%
  mutate(Chr = factor(Chr, levels = autosomes)) %>%
  arrange(Chr) %>%
  mutate(cum_end = cumsum(as.numeric(max_pos)),
         cum_start = cum_end - max_pos)

axis_labels <- chr_layout %>%
  mutate(mid_point = cum_start + (max_pos / 2),
         chr_simple = paste0("Chr", 1:n()))

chr_label_map <- axis_labels %>% select(Chr, chr_simple)

# --- Bar plot: average heterozygosity per chromosome, with genome-wide
# theta values baked directly into the legend labels ---

avg_het_data <- all_het_data %>%
  filter(Chr %in% autosomes) %>%
  mutate(Sample = ifelse(Sample == "Gleachii", "G. leachii", Sample)) %>%
  group_by(Sample, Chr) %>%
  summarize(mean_het = mean(Het, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(chr_label_map, by = "Chr")

sample_het_order <- avg_het_data %>%
  group_by(Sample) %>%
  summarize(overall_mean_het = mean(mean_het, na.rm = TRUE))

ingroup_order <- sample_het_order %>%
  filter(Sample != "G. leachii") %>%
  arrange(overall_mean_het) %>%
  pull(Sample)

new_sample_levels <- c(ingroup_order, "G. leachii")
avg_het_data$Sample <- factor(avg_het_data$Sample, levels = new_sample_levels)
avg_het_data$chr_simple <- factor(avg_het_data$chr_simple, levels = chr_label_map$chr_simple)

num_ingroup <- length(ingroup_order)
ingroup_colors <- brewer.pal(num_ingroup, "Set2")
new_color_values <- c(ingroup_colors, "#808080")
names(new_color_values) <- new_sample_levels

# Build legend labels with genome-wide theta appended, matching the
# genome_wide_het table computed above (Task 1). All labels are built via
# plotmath expressions (bquote) for consistent theta-symbol rendering --
# mixing a literal Unicode theta character with plotmath expressions
# causes the Unicode character to render incorrectly (font-dependent) in
# some graphics devices, even though the plotmath version renders fine.
het_lookup <- genome_wide_het %>%
  mutate(Sample = ifelse(Sample == "Gleachii", "G. leachii", Sample))

label_expr_for <- function(s, italic_name = FALSE) {
  val <- het_lookup$HetEst[het_lookup$Sample == s]
  val_str <- sprintf("%.4f", val)
  if (italic_name) {
    bquote(italic(.(s)) ~ "(" * theta * "=" * .(val_str) * ")")
  } else {
    bquote(.(s) ~ "(" * theta * "=" * .(val_str) * ")")
  }
}

new_legend_labels <- do.call(expression, c(
  lapply(ingroup_order, label_expr_for, italic_name = FALSE),
  list(label_expr_for("G. leachii", italic_name = TRUE))
))

p_bars <- ggplot(avg_het_data, aes(x = chr_simple, y = mean_het, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = new_color_values, labels = new_legend_labels) +
  labs(
    title = NULL,
    x = NULL,
    y = expression("Genetic variation (" * theta * ")"),
    fill = "Sample"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, angle = 0, vjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 10)),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(15, 15, 15, 15)
  )

output_plot_bars <- file.path(base_dir, "avg_heterozygosity_by_chr_bars.pdf")
ggsave(output_plot_bars, p_bars, width = 16, height = 7, limitsize = FALSE)
