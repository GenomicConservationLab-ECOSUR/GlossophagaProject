library(tidyverse)
library(ggrepel)

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_island")

chr_map <- c(
  "CM077303.1" = "1",  "CM077304.1" = "2",  "CM077305.1" = "3",
  "CM077306.1" = "4",  "CM077307.1" = "5",  "CM077308.1" = "6",
  "CM077309.1" = "7",  "CM077310.1" = "8",  "CM077311.1" = "9",
  "CM077312.1" = "10", "CM077313.1" = "11", "CM077314.1" = "12",
  "CM077315.1" = "13", "CM077316.1" = "14", "CM077317.1" = "15"
)

# 1. Load background data
dxy_data <- read_tsv("Island_vs_Mainland_dxy.txt", show_col_types = FALSE) %>%
  filter((pop1 == "Island" & pop2 == "Mainland") | (pop1 == "Mainland" & pop2 == "Island")) %>%
  filter(!is.na(avg_dxy)) %>%
  mutate(CHR_ID = as.numeric(recode(chromosome, !!!chr_map))) %>%
  filter(!is.na(CHR_ID)) %>%
  arrange(CHR_ID, window_pos_1)

# 2. Calculate cumulative genome-wide coordinates
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

# 3. Load corrected annotations, keeping ONLY genuine protein-coding
# candidates -- i.e. exactly the same filter used to build the publication
# table. TE hits, below-significance hits, and non-coding (lncRNA) windows
# are excluded entirely from the plot, not just from labeling, so the
# figure and the table always show the same 1:1 set of candidates.
highlight_genes <- read_csv("Validated_Island_Annotations_Corrected.csv", show_col_types = FALSE) %>%
  filter(Is_Significant, !Is_TE_Hit, !No_Blast_Hit) %>%
  mutate(
    CHR_ID = as.numeric(str_remove(Chromosome, "Chr")),
    window_pos_1 = Pos_Mb * 1000000,
    Plot_Label = if_else(Is_Named_Gene_Mismatch, UniProt_GN, Gene_Symbol_Final)
  ) %>%
  inner_join(data_cum, by = "CHR_ID") %>%
  mutate(bp_cum = window_pos_1 + bp_add)

# 3a. Suppress labels for locus tags that passed significance filtering
# but never resolved to a real gene name (no GN= field in their matching
# UniProt entry) -- the point stays plotted, only the uninformative raw
# tag text is dropped (matches the environmental clusters' analysis).
highlight_genes <- highlight_genes %>%
  mutate(Plot_Label = if_else(str_detect(Plot_Label, "^AAES06"), NA_character_, Plot_Label))

# 3b. Collapse the KRTAP paralog cluster into a single combined label, since
# labeling each paralog separately reads as clutter despite being a real,
# biologically meaningful signal. Every individual KRTAP point remains
# plotted; only the repeated text labels are consolidated.
label_data <- highlight_genes %>% filter(!is.na(Plot_Label))
krtap_cluster <- label_data %>% filter(str_detect(Plot_Label, "^KRTAP"))

if (nrow(krtap_cluster) > 1) {
  krtap_summary <- krtap_cluster %>%
    summarise(
      bp_cum = median(bp_cum),
      Max_Dxy = max(Max_Dxy),
      Plot_Label = sprintf("KRTAP10 cluster (n=%d)", n())
    )
  label_data <- label_data %>%
    filter(!str_detect(Plot_Label, "^KRTAP")) %>%
    bind_rows(krtap_summary)
  message(sprintf("KRTAP collapse applied: %d individual labels -> 1 combined label", nrow(krtap_cluster)))
} else {
  message(sprintf("KRTAP collapse NOT applied (%d matching row(s) found -- need >1 to collapse)", nrow(krtap_cluster)))
}

# 4. Generate the plot
p <- ggplot(plot_data, aes(x = bp_cum, y = avg_dxy)) +
  geom_point(aes(color = as.factor(CHR_ID %% 2)), alpha = 0.25, size = 0.8) +
  scale_color_manual(values = c("grey75", "grey45")) +

  geom_point(data = highlight_genes, aes(x = bp_cum, y = Max_Dxy),
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
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  labs(
    x = "Chromosome ID",
    y = expression(bold(D[XY]))
  )

print(p)
ggsave("Manhattan_island.pdf", p, width = 12, height = 6)
