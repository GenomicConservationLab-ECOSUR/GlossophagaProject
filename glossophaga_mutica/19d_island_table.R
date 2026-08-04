library(tidyverse)

# Build the base candidate-window annotation table (position, gene symbol,
# Dxy, pi, window count) from pixy output and the intersected GFF windows.
# Anonymous NCBI automated locus tags (AAES06_xxxxxx) are retained here
# rather than filtered out, since 19e_island_annotate_uniprot_function.R
# may "rescue" them with a real gene identity via BLAST homology.
# This is the input to 19e, which merges in BLAST/UniProt identity and
# function -- do not merge BLAST results here (see 19e for the corrected
# merge logic: proper accession parsing and TE/named-gene-mismatch
# handling).

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_island")

chr_map <- c(
  "CM077303.1" = "1",  "CM077304.1" = "2",  "CM077305.1" = "3",
  "CM077306.1" = "4",  "CM077307.1" = "5",  "CM077308.1" = "6",
  "CM077309.1" = "7",  "CM077310.1" = "8",  "CM077311.1" = "9",
  "CM077312.1" = "10", "CM077313.1" = "11", "CM077314.1" = "12",
  "CM077315.1" = "13", "CM077316.1" = "14", "CM077317.1" = "15"
)

dxy_stats <- read_tsv("Island_vs_Mainland_dxy.txt", show_col_types = FALSE) %>%
  filter((pop1 == "Island" & pop2 == "Mainland") | (pop1 == "Mainland" & pop2 == "Island")) %>%
  filter(!is.na(avg_dxy)) %>%
  select(chromosome, window_pos_1, avg_dxy)

pi_stats <- read_tsv("Island_vs_Mainland_pi.txt", show_col_types = FALSE) %>%
  filter(pop == "Island") %>%
  filter(!is.na(avg_pi)) %>%
  select(chromosome, window_pos_1, avg_pi)

raw_bed <- read_tsv("annotated_island_windows.txt", col_names = FALSE, show_col_types = FALSE)

tidy_annotations <- raw_bed %>%
  select(chromosome = X1, window_pos_1 = X2, window_pos_2 = X3, metadata = X13) %>%
  mutate(
    Gene_Symbol = case_when(
      str_detect(metadata, "Name=") ~ str_match(metadata, "Name=([^;]+)")[, 2],
      str_detect(metadata, "locus_tag=") ~ str_match(metadata, "locus_tag=([^;]+)")[, 2],
      TRUE ~ "Intergenic"
    ),
    Description = case_when(
      str_detect(metadata, "description=") ~ str_match(metadata, "description=([^;]+)")[, 2],
      str_detect(metadata, "gene_biotype=lncRNA") ~ "Long non-coding RNA element",
      Gene_Symbol == "Intergenic" ~ "Cis-regulatory domain or uncharacterized intergenic space",
      TRUE ~ "Protein-coding gene model"
    )
  ) %>%
  distinct(chromosome, window_pos_1, Gene_Symbol, Description)

table1_island <- tidy_annotations %>%
  inner_join(dxy_stats, by = c("chromosome", "window_pos_1")) %>%
  inner_join(pi_stats, by = c("chromosome", "window_pos_1")) %>%
  mutate(
    CHR_NUM = as.numeric(recode(chromosome, !!!chr_map)),
    Chromosome = paste0("Chr", CHR_NUM),
    Position_Mb = round(window_pos_1 / 1000000, 3)
  ) %>%
  filter(!is.na(CHR_NUM)) %>%
  group_by(CHR_NUM, Chromosome, Gene_Symbol, Description) %>%
  summarise(
    Pos_Mb = min(Position_Mb),
    Max_Dxy = round(max(avg_dxy), 4),
    Avg_Pi = round(mean(avg_pi), 4),
    Block_Width_Windows = n(),
    .groups = "drop"
  ) %>%
  select(Chromosome, Pos_Mb, Gene_Symbol, Description, Max_Dxy, Avg_Pi, Block_Width_Windows, CHR_NUM) %>%
  arrange(CHR_NUM, Pos_Mb) %>%
  select(-CHR_NUM)

write_csv(table1_island, "Table_Island_vs_Mainland_Annotations.csv")
print(as.data.frame(table1_island), row.names = FALSE)
