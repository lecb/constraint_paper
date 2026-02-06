#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ============================================================
# gene_level_analysis.R (fixed)
#
# Your run log shows this script already:
# - reads ./gnomAD_genes/*.csv
# - writes Table1_*.csv, Supplementary_Table1_full.csv
# - writes compare_discordant_vs_rest_constraint_gene_level.csv
# - writes Figure2_discordant_vs_rest_constraint.png
# Then it crashes because gene_level_compare is not defined,
# but is referenced ONLY for end-of-script QC prints.
#
# This fixed version ensures gene_level_compare exists by
# loading it from the CSV you already wrote (if needed),
# so the script completes cleanly.
# ============================================================

# ----------------------------
# User-configurable paths
# ----------------------------
GENE_DIR <- "gnomAD_genes"
COMPARE_GENE_LEVEL_CSV <- "compare_discordant_vs_rest_constraint_gene_level.csv"
FIG2_PNG <- "Figure2_discordant_vs_rest_constraint.png"

# ----------------------------
# Basic preflight
# ----------------------------
message("Working directory: ", getwd())

if (!dir.exists(GENE_DIR)) {
  stop("Missing folder: ", GENE_DIR)
}

gene_files <- list.files(GENE_DIR, pattern = "\\.csv$", full.names = TRUE)
message("Found ", length(gene_files), " gene CSV files in ./'", GENE_DIR, "'")
if (length(gene_files) == 0) stop("No gene CSV files found in ", GENE_DIR)

message("Example files: ", paste(head(basename(gene_files), 5), collapse = ", "))

# ============================================================
# IMPORTANT:
# This script assumes you already have your existing analysis
# code that produces all the outputs listed in your log.
#
# Since I don't have your full original script content here,
# I will NOT try to re-implement your whole logic from scratch.
#
# Instead, you should KEEP your existing analysis code
# exactly as-is, and add the PATCH BLOCK below right after
# you write COMPARE_GENE_LEVEL_CSV (and before QC prints).
#
# To make this file "working" as a standalone drop-in,
# I include a small guard that requires the compare CSV
# and Figure2 to already exist (i.e., produced earlier in
# your script). In your real file, those are produced above.
# ============================================================

# ------------------------------------------------------------------
# >>>>>>>>>>>>>>>>>>>>>>>>>>> PATCH BLOCK (START) <<<<<<<<<<<<<<<<<<<<
# Put this block immediately AFTER you write:
#   compare_discordant_vs_rest_constraint_gene_level.csv
# and BEFORE any line that references gene_level_compare.
# ------------------------------------------------------------------

if (!exists("gene_level_compare")) {
  if (file.exists(COMPARE_GENE_LEVEL_CSV)) {
    message("[FIX] gene_level_compare not found; loading from ", COMPARE_GENE_LEVEL_CSV)
    gene_level_compare <- read_csv(COMPARE_GENE_LEVEL_CSV, show_col_types = FALSE)
  } else {
    message("[FIX] gene_level_compare not found AND compare CSV not found.")
    message("[FIX] This means the compare table was not created earlier as expected.")
    stop("Cannot continue without gene_level_compare or ", COMPARE_GENE_LEVEL_CSV)
  }
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>> PATCH BLOCK (END) <<<<<<<<<<<<<<<<<<<<

# ============================================================
# End-of-script QC prints (these were crashing you)
# ============================================================

message("\n[QC] gene_level_compare summary:")
if ("group" %in% names(gene_level_compare)) {
  print(table(gene_level_compare$group, useNA = "ifany"))
} else {
  message("[QC] Column 'group' not present in gene_level_compare (available columns: ",
          paste(names(gene_level_compare), collapse = ", "), ")")
}

if ("trunc_n_variants" %in% names(gene_level_compare)) {
  print(summary(gene_level_compare$trunc_n_variants))
  message("[QC] NAs in trunc_n_variants: ", sum(is.na(gene_level_compare$trunc_n_variants)))
} else {
  message("[QC] Column 'trunc_n_variants' not present in gene_level_compare (available columns: ",
          paste(names(gene_level_compare), collapse = ", "), ")")
}

message("\n[DONE] gene_level_analysis.R completed cleanly.")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# Inputs
UNIVERSE_CONSTRAINT <- "gnomad_lof_discordance_out/loeuf_lt_0.2_precision_filtered.csv"
DISEASE_EVIDENCE    <- "tables/Supp_annotated_universe.csv"

# Output
OUT_COMBINED <- "tables/Supp_combined_universe_flags_disease.csv"

# ---- load ----
u <- read_csv(UNIVERSE_CONSTRAINT, show_col_types = FALSE)
e <- read_csv(DISEASE_EVIDENCE,    show_col_types = FALSE)

# ---- standardise gene key ----
# constraint universe typically uses 'gene'
stopifnot("gene" %in% names(u))

# disease evidence uses 'Gene'
stopifnot("Gene" %in% names(e))

u2 <- u %>%
  mutate(Gene = as.character(gene),
         Gene_key = str_to_upper(str_trim(Gene)))

e2 <- e %>%
  mutate(Gene = as.character(Gene),
         Gene_key = str_to_upper(str_trim(Gene)))

# ---- join ----
# keep all LOEUF<0.2 universe genes (left join)
# ---- join ----
combined <- u2 %>%
  left_join(
    e2 %>% select(Gene_key, discordant, PanelApp_green, ClinGen_HI, HI_score, Any_disease_evidence),
    by = "Gene_key"
  ) %>%
  mutate(
    discordant = ifelse(is.na(discordant), FALSE, as.logical(discordant)),
    constraint_flags = coalesce(as.character(constraint_flags), ""),
    constraint_outlier_flag_present = ifelse(
      str_trim(constraint_flags) == "" | constraint_flags %in% c("[]", "[ ]"),
      FALSE, TRUE
    )
  ) %>%
  select(
    Gene,  # <- this exists (from u2)
    discordant,
    constraint_flags,
    constraint_outlier_flag_present,
    PanelApp_green,
    ClinGen_HI,
    HI_score,
    Any_disease_evidence,
    everything()
  ) %>%
  select(-Gene_key)


# ---- sanity checks ----
# Discordant genes should all be present if the two files correspond to same universe
n_discordant <- sum(combined$discordant, na.rm = TRUE)
message("[QC] Combined universe rows: ", nrow(combined))
message("[QC] Discordant count: ", n_discordant)

# Any NAs in evidence columns? that's okay if evidence file had fewer rows than universe,
# but usually should match. Still, keep explicit NAs rather than silently changing.
write_csv(combined, OUT_COMBINED, na = "")

message("[DONE] Wrote: ", OUT_COMBINED)

