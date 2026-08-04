library(tidyverse)
library(httr)
library(jsonlite)

# Annotate candidate windows for all five environmental clusters with real
# UniProt Function text. Applies the same corrections validated on the
# Island vs. Mainland analysis:
# 1. Proper UniProt accession parsing (sp|ACCESSION|NAME split on "|",
#    not a fixed 6-character regex, which truncates newer 10-character
#    accessions).
# 2. TE hits flagged by scanning UniProt_Description for "LINE-1" /
#    "retrotransposable" text -- not by matching GN=Pol/L1RE1, which fails
#    whenever the matching UniProt entry lacks a GN= field at all (in which
#    case the hit was previously silently mislabeled with the nearby gene's
#    symbol instead of being flagged as a TE).
# 3. Locus-tag rescue (anonymous AAES06_xxxxxx tags renamed via confident
#    BLAST hits) vs. named-gene mismatch (a reference-annotated gene whose
#    best hit is a DIFFERENT named gene -- both identities retained).
# 4. Significance filtering: a hit must pass BOTH an E-value and a minimum
#    alignment length threshold to be trusted (E-value alone is not
#    sufficient -- see Discussion of the Island analysis for the specific
#    cases this caught: a short immunoglobulin-domain fragment, a
#    bacteriophage protein, and an Arabidopsis protein, all with weak or
#    spurious support despite superficially plausible top hits).
# 5. Duplicate windows: candidate windows overlapping multiple GFF gene
#    models that resolve to the same final identity (matched
#    case-insensitively) are collapsed, keeping the strongest hit.

setwd("/run/media/rocamontes/EcosurMutica/Pixy/Pixy_2026/output_environmental")

CLUSTERS <- c("Concha", "Hobo_Ixta", "Nizanda", "Poana", "Quilamula")

# Shared cache across all clusters (and reusable across the Island analysis
# too, if pointed at the same file) -- avoids re-querying UniProt for
# accessions that recur across multiple clusters (e.g. common LINE-1 hits).
CACHE_FILE <- "uniprot_function_cache.csv"
REQUEST_DELAY_SEC <- 0.5

EVALUE_THRESHOLD <- 1e-5
MIN_ALIGN_LENGTH <- 35

fetch_uniprot_function <- function(accession) {
  url <- sprintf("https://rest.uniprot.org/uniprotkb/%s.json", accession)
  resp <- tryCatch(GET(url), error = function(e) NULL)
  if (is.null(resp) || status_code(resp) != 200) return(NA_character_)
  parsed <- tryCatch(fromJSON(content(resp, "text", encoding = "UTF-8")), error = function(e) NULL)
  if (is.null(parsed) || is.null(parsed$comments)) return(NA_character_)
  func_comments <- parsed$comments[parsed$comments$commentType == "FUNCTION", ]
  if (nrow(func_comments) == 0) return(NA_character_)
  paste(func_comments$texts[[1]]$value, collapse = " ")
}

if (file.exists(CACHE_FILE)) {
  function_cache <- read_csv(CACHE_FILE, show_col_types = FALSE)
} else {
  function_cache <- tibble(UniProt_Accession = character(), UniProt_Function = character())
}

for (cluster in CLUSTERS) {

  table_file  <- paste0("Cluster_", cluster, "_Robust_Candidates_Top200.csv")
  blast_file  <- paste0("specific_genes_blast_", cluster, ".txt")
  anno_file   <- paste0("annotated_", cluster, "_windows.txt")
  output_file <- paste0("Validated_", cluster, "_Annotations_Corrected.csv")

  if (!(file.exists(table_file) && file.exists(blast_file) && file.exists(anno_file))) {
    cat("\nSkipping", cluster, "- Missing CSV, BLAST, or Annotation file.\n")
    next
  }

  cat("\n=== Processing:", cluster, "===\n")

  cluster_table <- read_csv(table_file, show_col_types = FALSE) %>%
    mutate(window_pos_1 = as.numeric(window_pos_1))

  raw_anno <- read_tsv(anno_file, col_names = FALSE, show_col_types = FALSE)
  if (nrow(raw_anno) == 0) {
    cat("  - No gene annotations found - candidate windows fall in intergenic regions. Skipping.\n")
    next
  }
  colnames(raw_anno)[c(1, 2, 3, 13)] <- c("chromosome", "window_pos_1", "window_pos_2", "attributes")

  gene_map <- raw_anno %>%
    select(chromosome, window_pos_1, window_pos_2, attributes) %>%
    mutate(
      Gene_Symbol = case_when(
        str_detect(attributes, "Name=") ~ str_match(attributes, "Name=([^;]+)")[, 2],
        str_detect(attributes, "locus_tag=") ~ str_match(attributes, "locus_tag=([^;]+)")[, 2],
        TRUE ~ "Uncharacterized"
      ),
      Description = case_when(
        str_detect(attributes, "description=") ~ str_match(attributes, "description=([^;]+)")[, 2],
        str_detect(attributes, "gene_biotype=lncRNA") ~ "Long non-coding RNA element",
        TRUE ~ "Protein-coding gene model"
      )
    ) %>%
    filter(Gene_Symbol != "Uncharacterized") %>%
    distinct(chromosome, window_pos_1, Gene_Symbol, Description) %>%
    mutate(window_pos_1 = as.numeric(window_pos_1))

  annotated_table <- cluster_table %>%
    left_join(gene_map, by = c("chromosome", "window_pos_1")) %>%
    filter(!is.na(Gene_Symbol))

  if (nrow(annotated_table) == 0) {
    cat("  - Candidate windows exist but no gene overlap found. Skipping.\n")
    next
  }

  blast_raw <- read_tsv(blast_file,
                        col_names = c("qseqid", "sseqid", "pident", "length",
                                      "evalue", "bitscore", "stitle"),
                        show_col_types = FALSE)

  if (nrow(blast_raw) == 0) {
    cat("  - No BLAST hits found. Retaining original NCBI labels with no annotation.\n")
    final_table <- annotated_table %>%
      mutate(Gene_Symbol_Final = Gene_Symbol,
             Function_Note = "No BLASTX hit obtained for this locus.",
             Is_TE_Hit = FALSE, Is_Significant = FALSE, No_Blast_Hit = TRUE,
             Is_Named_Gene_Mismatch = FALSE)
    write_csv(final_table, output_file)
    next
  }

  blast_parsed <- blast_raw %>%
    mutate(UniProt_Accession = map_chr(str_split(sseqid, "\\|"), ~ .x[2])) %>%
    select(qseqid, UniProt_Accession, UniProt_Description = stitle,
           E_value = evalue, Align_Length = length) %>%
    group_by(qseqid) %>%
    slice_min(E_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      Is_Significant = (E_value < EVALUE_THRESHOLD) & (Align_Length >= MIN_ALIGN_LENGTH),
      # Broadened beyond "LINE-1|retrotransposable" to catch other
      # repeat-/retroelement-derived hits (e.g. retroviral Gag-Pro-Pol
      # polyproteins), which are not host-gene divergence signals either.
      Is_TE_Hit = Is_Significant & str_detect(UniProt_Description,
        regex("LINE-1|retrotransposable|Gag-Pro-Pol|polyprotein|retrovirus|reverse transcriptase|transposon",
              ignore_case = TRUE)),
      UniProt_GN = str_match(UniProt_Description, "GN=(\\S+)")[, 2],
      Is_Locus_Tag = str_detect(qseqid, "^AAES06"),
      Is_Named_Gene_Mismatch = Is_Significant & !Is_TE_Hit & !Is_Locus_Tag &
        !is.na(UniProt_GN) & (str_to_upper(UniProt_GN) != str_to_upper(qseqid)),
      Gene_Symbol_Final = case_when(
        !Is_Significant ~ qseqid,
        Is_TE_Hit ~ qseqid,
        Is_Named_Gene_Mismatch ~ sprintf("%s (annotated as %s in reference genome)", UniProt_GN, qseqid),
        Is_Locus_Tag & !is.na(UniProt_GN) ~ UniProt_GN,
        TRUE ~ qseqid
      )
    )

  n_below <- sum(!blast_parsed$Is_Significant)
  cat(sprintf("  - %d best-hit records: %d below threshold, %d TE hits, %d named-gene mismatches\n",
              nrow(blast_parsed), n_below, sum(blast_parsed$Is_TE_Hit), sum(blast_parsed$Is_Named_Gene_Mismatch)))

  # Fetch UniProt Function for any new, significant, non-TE accessions
  accessions_to_fetch <- blast_parsed %>%
    filter(Is_Significant & !Is_TE_Hit) %>%
    distinct(UniProt_Accession) %>%
    filter(!UniProt_Accession %in% function_cache$UniProt_Accession) %>%
    pull(UniProt_Accession)

  if (length(accessions_to_fetch) > 0) {
    cat(sprintf("  - Fetching UniProt Function for %d new accessions...\n", length(accessions_to_fetch)))
    for (acc in accessions_to_fetch) {
      func_text <- fetch_uniprot_function(acc)
      function_cache <- bind_rows(function_cache, tibble(UniProt_Accession = acc, UniProt_Function = func_text))
      Sys.sleep(REQUEST_DELAY_SEC)
    }
    write_csv(function_cache, CACHE_FILE)
  }

  blast_annotated <- blast_parsed %>%
    left_join(function_cache, by = "UniProt_Accession") %>%
    mutate(
      Function_Note = case_when(
        !Is_Significant ~ sprintf(
          "No significant BLAST hit for %s (best available hit did not meet significance thresholds: E<%.0e, length>=%d aa).",
          qseqid, EVALUE_THRESHOLD, MIN_ALIGN_LENGTH),
        Is_TE_Hit ~ sprintf(
          "Divergent region overlaps %s; best BLAST hit is a LINE-1 retrotransposon protein, consistent with repetitive-element content rather than host-gene divergence.",
          qseqid),
        Is_Named_Gene_Mismatch ~ sprintf(
          "Reference genome annotates this locus as %s; best BLAST hit is %s. %s",
          qseqid, UniProt_GN, coalesce(UniProt_Function, "No curated Function annotation available for this UniProt entry.")),
        !is.na(UniProt_Function) ~ UniProt_Function,
        TRUE ~ "No curated Function annotation available for this UniProt entry."
      )
    )

  final_table <- annotated_table %>%
    left_join(blast_annotated, by = c("Gene_Symbol" = "qseqid")) %>%
    mutate(
      No_Blast_Hit = is.na(Is_TE_Hit),
      Gene_Symbol_Final = if_else(No_Blast_Hit, Gene_Symbol, Gene_Symbol_Final),
      Function_Note = case_when(
        No_Blast_Hit & str_detect(Description, regex("non-coding|lncRNA", ignore_case = TRUE)) ~
          "Long non-coding RNA element; no protein-coding sequence, so no BLASTX hit is expected.",
        No_Blast_Hit ~ "No BLASTX hit obtained for this locus.",
        TRUE ~ Function_Note
      )
    )

  # Collapse duplicate windows sharing the same resolved identity
  # (case-insensitive), keeping the strongest (lowest E-value) hit.
  n_before <- nrow(final_table)
  excluded_rows <- final_table %>% filter(No_Blast_Hit | !Is_Significant)
  hit_rows <- final_table %>%
    filter(!No_Blast_Hit & Is_Significant) %>%
    mutate(.dedup_key = str_to_upper(Gene_Symbol_Final)) %>%
    group_by(chromosome, window_pos_1, .dedup_key) %>%
    slice_min(E_value, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-.dedup_key)
  final_table <- bind_rows(hit_rows, excluded_rows)
  n_after <- nrow(final_table)
  if (n_before != n_after) {
    cat(sprintf("  - Collapsed %d duplicate window(s) (%d -> %d rows)\n", n_before - n_after, n_before, n_after))
  }

  write_csv(final_table, output_file)
  cat("  - Saved:", output_file, "\n")
}

cat("\nAll environmental clusters processed.\n")
