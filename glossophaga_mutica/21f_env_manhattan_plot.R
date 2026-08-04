library(tidyverse)
library(ggrepel)

# Manhattan plots for environmentally-defined mainland clusters with
# gene-overlapping candidate windows (Concha, Hobo_Ixta, Nizanda). Poana
# and Quilamula are skipped automatically -- confirmed genuine gene
# deserts, no candidate windows overlap an annotated gene model.
#
# Mirrors the island Manhattan script (only genuine, significant,
# non-TE, protein-coding candidates are plotted at all), with one addition:
# gene-family clusters (e.g. KRTAP, olfactory receptors, zinc-fingers) are
# collapsed into a single combined label, since labeling each paralog
# separately reads as clutter despite being a real signal. Every individual
# point remains plotted; only text labels are consolidated.
#
# IMPORTANT: collapsing is based on shared gene-name root AND multiple
# DISTINCT gene symbols -- not just multiple rows. A gene spanning two
# adjacent, overlapping candidate windows (e.g. TAF15 appearing at two
# nearby positions) is the same gene, not a paralog family, and must not
# be mislabeled as a "cluster".

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_environmental")

CLUSTERS <- c("Concha", "Hobo_Ixta", "Nizanda")

chr_map <- c(
  "CM077303.1" = "1",  "CM077304.1" = "2",  "CM077305.1" = "3",
  "CM077306.1" = "4",  "CM077307.1" = "5",  "CM077308.1" = "6",
  "CM077309.1" = "7",  "CM077310.1" = "8",  "CM077311.1" = "9",
  "CM077312.1" = "10", "CM077313.1" = "11", "CM077314.1" = "12",
  "CM077315.1" = "13", "CM077316.1" = "14", "CM077317.1" = "15"
)

for (cluster in CLUSTERS) {

  dxy_file <- paste0("Cluster_", cluster, "_vs_others_dxy.txt")
  anno_file <- paste0("Validated_", cluster, "_Annotations_Corrected.csv")

  if (!file.exists(anno_file)) {
    cat(sprintf("Skipping %s -- no corrected annotation file.\n", cluster))
    next
  }

  # 1. Load background data
  dxy_data <- read_tsv(dxy_file, show_col_types = FALSE) %>%
    filter(!is.na(avg_dxy)) %>%
    mutate(CHR_ID = as.numeric(recode(chromosome, !!!chr_map))) %>%
    filter(!is.na(CHR_ID)) %>%
    arrange(CHR_ID, window_pos_1)

  # 2. Cumulative genome-wide coordinates
  data_cum <- dxy_data %>%
    group_by(CHR_ID) %>%
    summarise(max_bp = max(window_pos_1)) %>%
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>%
    select(CHR_ID, bp_add)

  plot_data <- dxy_data %>%
    inner_join(data_cum, by = "CHR_ID") %>%
    mutate(bp_cum = window_pos_1 + bp_add)

  axis_set <- plot_data %>%
    group_by(CHR_ID) %>%
    summarize(center = mean(bp_cum))

  # 3. Load corrected annotations, genuine protein-coding candidates only
  highlight_genes <- read_csv(anno_file, show_col_types = FALSE) %>%
    filter(Is_Significant, !Is_TE_Hit, !No_Blast_Hit) %>%
    mutate(
      CHR_ID = as.numeric(recode(chromosome, !!!chr_map)),
      Plot_Label = if_else(Is_Named_Gene_Mismatch, UniProt_GN, Gene_Symbol_Final)
    ) %>%
    inner_join(data_cum, by = "CHR_ID") %>%
    mutate(bp_cum = window_pos_1 + bp_add)

  if (nrow(highlight_genes) == 0) {
    cat(sprintf("Skipping %s -- no genuine candidates after filtering.\n", cluster))
    next
  }

  # 3b. Generalized gene-family cluster collapsing.
  # First, suppress labels for locus tags that passed significance
  # filtering but never resolved to a real gene name (no GN= field in
  # their matching UniProt entry) -- the point stays plotted, only the
  # uninformative raw tag text is dropped (matches the island analysis).
  highlight_genes <- highlight_genes %>%
    mutate(Plot_Label = if_else(str_detect(Plot_Label, "^AAES06"), NA_character_, Plot_Label))

  # family_root: uppercase gene symbol with trailing digits/hyphens
  # stripped (e.g. "OR2A14" -> "OR2A", "KRTAP10-7" -> "KRTAP",
  # "ZNF432" -> "ZNF"). Only collapse when a root has multiple DISTINCT
  # gene symbols (true paralogs) -- not merely multiple rows, which can
  # happen when one gene spans two adjacent, overlapping windows.
  label_data <- highlight_genes %>%
    filter(!is.na(Plot_Label)) %>%
    mutate(family_root = str_remove(str_to_upper(Plot_Label), "[0-9-]+$"))

  family_sizes <- label_data %>%
    distinct(family_root, Plot_Label) %>%
    count(family_root) %>%
    filter(n > 1)

  if (nrow(family_sizes) > 0) {
    collapsed_families <- label_data %>%
      filter(family_root %in% family_sizes$family_root) %>%
      group_by(family_root) %>%
      summarise(
        bp_cum = median(bp_cum),
        Max_Dxy = max(avg_dxy),
        Plot_Label = sprintf("%s cluster (n=%d)", first(family_root), n_distinct(Plot_Label)),
        .groups = "drop"
      )
    individual_labels <- label_data %>%
      filter(!family_root %in% family_sizes$family_root) %>%
      mutate(Max_Dxy = avg_dxy) %>%
      as_tibble()
    label_data <- bind_rows(
      individual_labels[, c("bp_cum", "Max_Dxy", "Plot_Label")],
      collapsed_families[, c("bp_cum", "Max_Dxy", "Plot_Label")]
    )
  } else {
    label_data <- label_data %>%
      mutate(Max_Dxy = avg_dxy) %>%
      as_tibble()
    label_data <- label_data[, c("bp_cum", "Max_Dxy", "Plot_Label")]
  }

  # 4. Generate the plot
  p <- ggplot(plot_data, aes(x = bp_cum, y = avg_dxy)) +
    geom_point(aes(color = as.factor(CHR_ID %% 2)), alpha = 0.25, size = 0.8) +
    scale_color_manual(values = c("grey75", "grey45")) +

    geom_point(data = highlight_genes, aes(x = bp_cum, y = avg_dxy),
               color = "firebrick", size = 2.5, shape = 18) +

    geom_text_repel(data = label_data, aes(x = bp_cum, y = Max_Dxy, label = Plot_Label),
                    color = "black", fontface = "bold.italic", size = 2.8,
                    box.padding = 0.8,
                    point.padding = 0.4,
                    max.overlaps = 50,
                    min.segment.length = 0,
                    segment.color = "grey30",
                    segment.size = 0.3,
                    force = 4,
                    nudge_y = 0.05,
                    direction = "both") +

    scale_x_continuous(label = axis_set$CHR_ID, breaks = axis_set$center) +
    scale_y_continuous(limits = c(0, max(plot_data$avg_dxy, na.rm = TRUE) + 0.05), expand = c(0, 0)) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 14, face = "bold")
    ) +
    labs(x = "Chromosome ID", y = expression(bold(D[XY])))

  print(p)
  ggsave(paste0("Manhattan_", cluster, ".pdf"), p, width = 12, height = 6)
  cat(sprintf("Saved: Manhattan_%s.pdf\n", cluster))
}

cat("\nAll clusters processed.\n")
