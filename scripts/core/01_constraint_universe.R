#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(stringr)
  library(officer)
  library(flextable)
})

# ============================================================
# gnomAD constraint: LOEUF<0.2 + high LoF o/e point estimate
# LOEUF defined as lof.oe_ci.upper
# Produces:
#   - tables/Table1_word_ready_top15.csv
#   - tables/Table1_word_ready_top15.docx
#   - tables/Supplementary_Table1_discordant_genes.csv
#   - gnomad_lof_discordance_out/cache/gene_features.rds   (NEW)
# KEEP FOR PAPER
# ============================================================

# ---- helpers (define EARLY!) ----
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# ---- gene-level row selector (constraint TSV may be transcript-level) ----
pick_one_tx_per_gene <- function(df_tx) {
  # Prefer MANE Select, then canonical, otherwise highest oe_lof (stable, deterministic)
  df_tx %>%
    arrange(
      desc(!is.na(mane_select) & mane_select %in% c(TRUE, "TRUE", "true", 1, "1")),
      desc(!is.na(canonical)   & canonical   %in% c(TRUE, "TRUE", "true", 1, "1")),
      desc(oe_lof)
    ) %>%
    slice(1)
}

# ----------------------------
# USER INPUTS
# ----------------------------
CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"
VARFILE <- file.path(Sys.getenv("PIPELINE_ROOT","."), "data", "gnomad", "gnomad_variants_all.csv")

OUTDIR  <- "gnomad_lof_discordance_out"
TABDIR  <- "tables"

LOEUF_CUTOFF <- 0.2
DROP_WORST_CI_WIDTH_FRAC <- 0.25
TOP_OE_FRAC <- 0.05
TOP_N <- 15L

WRITE_INTERMEDIATES <- FALSE  # set TRUE if you want the intermediate CSVs + plots in OUTDIR

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TABDIR, showWarnings = FALSE, recursive = TRUE)

CACHEDIR <- file.path(OUTDIR, "cache")
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Load + standardise constraint table
#    + include genome-wide observed counts from TSV
# ----------------------------
df <- read_tsv(CONSTRAINT_TSV, show_col_types = FALSE, progress = FALSE)

# Guard: ensure columns exist
need_constraint_cols <- c("gene","gene_id","transcript","canonical","mane_select","cds_length",
                          "constraint_flags","lof.oe","lof.oe_ci.lower","lof.oe_ci.upper",
                          "syn.obs","lof.obs")
missing_constraint_cols <- setdiff(need_constraint_cols, names(df))
if (length(missing_constraint_cols) > 0) {
  stop("Constraint TSV missing required columns: ",
       paste(missing_constraint_cols, collapse = ", "))
}

d <- df %>%
  transmute(
    gene = gene,
    gene_id = gene_id,
    transcript = transcript,
    canonical = canonical,
    mane_select = mane_select,
    cds_length = cds_length,
    constraint_flags = constraint_flags,
    
    # genome-wide observed counts (from TSV)
    syn_obs = as.numeric(`syn.obs`),
    lof_obs = as.numeric(`lof.obs`),
    
    # constraint metrics
    oe_lof = as.numeric(`lof.oe`),
    oe_lof_lower = as.numeric(`lof.oe_ci.lower`),
    oe_lof_upper = as.numeric(`lof.oe_ci.upper`),
    
    # LOEUF = LoF o/e upper CI (your schema)
    loeuf = as.numeric(`lof.oe_ci.upper`)
  ) %>%
  filter(!is.na(gene), !is.na(oe_lof), !is.na(oe_lof_lower), !is.na(oe_lof_upper), !is.na(loeuf))

stopifnot(nrow(d) > 0)

# ----------------------------
# 1b) Collapse to ONE row per gene (constraint TSV may be transcript-level)
# ----------------------------
constraint_genelevel <- d %>%
  group_by(gene) %>%
  group_modify(~ pick_one_tx_per_gene(.x)) %>%
  ungroup() %>%
  mutate(ci_width = oe_lof_upper - oe_lof_lower) %>%
  filter(!is.na(ci_width), ci_width >= 0)

stopifnot(nrow(constraint_genelevel) > 0)

# ----------------------------
# 1c) One-time build: gene_features.rds (genome-wide; cheap to reuse)
#     Note: trunc_obs here is gnomAD lof.obs (LOFTEE-based LoF), not strictly stop+frameshift.
# ----------------------------
gene_features <- constraint_genelevel %>%
  transmute(
    gene,
    gene_id,
    transcript,
    canonical,
    mane_select,
    cds_length,
    constraint_flags,
    
    loeuf,
    oe_lof,
    oe_lof_lower,
    oe_lof_upper,
    ci_width,
    
    # burdens (genome-wide from TSV)
    syn_obs = syn_obs,
    trunc_obs = lof_obs
  )

saveRDS(gene_features, file.path(CACHEDIR, "gene_features.rds"))

if (isTRUE(WRITE_INTERMEDIATES)) {
  write_csv(gene_features, file.path(OUTDIR, "gene_features.csv"))
}

message("[CACHE] Wrote gene_features.rds to: ", file.path(CACHEDIR, "gene_features.rds"))

# ----------------------------
# 2) Filter LOEUF<cutoff + precision metrics (uses transcript-level d, as before)
# ----------------------------
d0 <- d %>%
  filter(loeuf < LOEUF_CUTOFF) %>%
  mutate(ci_width = oe_lof_upper - oe_lof_lower) %>%
  filter(ci_width >= 0)

stopifnot(nrow(d0) > 0)

ci_cut <- as.numeric(quantile(d0$ci_width, probs = 1 - DROP_WORST_CI_WIDTH_FRAC, na.rm = TRUE))
d1 <- d0 %>% filter(ci_width <= ci_cut)
stopifnot(nrow(d1) > 0)

oe_cut <- as.numeric(quantile(d1$oe_lof, probs = 1 - TOP_OE_FRAC, na.rm = TRUE))

discordant <- d1 %>%
  filter(oe_lof >= oe_cut) %>%
  arrange(desc(oe_lof))

# collapse to one row per gene (constraint file may be transcript-level)
discordant_genelevel <- discordant %>%
  group_by(gene) %>%
  slice_max(order_by = oe_lof, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(oe_lof))

# ----------------------------
# Optional: write intermediates
# ----------------------------
if (isTRUE(WRITE_INTERMEDIATES)) {
  write_csv(d0, file.path(TABDIR, "loeuf_lt_0.2_all_genes_with_ci_metrics.csv"))
  write_csv(d1, file.path(OUTDIR, "loeuf_lt_0.2_precision_filtered.csv"))
  write_csv(discordant, file.path(TABDIR, "discordant_genes_top_oe_lof.csv"))
  
  qc <- tibble(
    loeuf_cutoff = LOEUF_CUTOFF,
    drop_worst_ci_width_frac = DROP_WORST_CI_WIDTH_FRAC,
    top_oe_frac = TOP_OE_FRAC,
    n_loeuf = nrow(d0),
    ci_width_cut = ci_cut,
    n_precision = nrow(d1),
    oe_cut = oe_cut,
    n_discordant = nrow(discordant_genelevel)
  )
  write_csv(qc, file.path(TABDIR, "qc_summary.csv"))
  
  p1 <- ggplot(d0, aes(x = loeuf, y = oe_lof)) +
    geom_point(alpha = 0.35, size = 1) +
    geom_point(data = discordant_genelevel, size = 1.8) +
    labs(
      title = "LOEUF<0.2 genes: LoF o/e point estimate vs LOEUF (upper CI)",
      subtitle = paste0(
        "Precision: CI width <= ", round(ci_cut, 3),
        " (drop worst ", DROP_WORST_CI_WIDTH_FRAC * 100, "%); ",
        "Discordant = top ", TOP_OE_FRAC * 100, "% by lof.oe (cut=", round(oe_cut, 3), ")"
      ),
      x = "LOEUF (lof.oe_ci.upper)",
      y = "LoF observed/expected (lof.oe point estimate)"
    ) +
    theme_minimal()
  ggsave(file.path(OUTDIR, "plot_oe_lof_vs_loeuf.png"), p1, width = 7, height = 5, dpi = 300)
  
  p2 <- ggplot(d0, aes(x = ci_width)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = ci_cut, linetype = "dashed") +
    labs(
      title = "CI width distribution (LOEUF<0.2 genes)",
      subtitle = paste0("Dashed cutoff = ", round(ci_cut, 3),
                        " (drop worst ", DROP_WORST_CI_WIDTH_FRAC * 100, "%)"),
      x = "CI width (lof.oe_ci.upper - lof.oe_ci.lower)",
      y = "Number of genes"
    ) +
    theme_minimal()
  ggsave(file.path(OUTDIR, "plot_ci_width_hist.png"), p2, width = 7, height = 5, dpi = 300)
}

# ----------------------------
# 3) Build Supplementary Table 1 (all discordant genes)
# ----------------------------
supp_table1 <- discordant_genelevel %>%
  transmute(
    Gene = gene,
    gene_id,
    transcript,
    canonical,
    mane_select,
    cds_length,
    constraint_flags,
    LOEUF = round(loeuf, 3),
    `LoF o/e` = round(oe_lof, 3),
    CI_lower = round(oe_lof_lower, 3),
    CI_upper = round(oe_lof_upper, 3),
    ci_width = round(ci_width, 3)
  )

write_csv(supp_table1, file.path(TABDIR, "Supplementary_Table1_discordant_genes.csv"), na = "")
if (isTRUE(WRITE_INTERMEDIATES)) {
  write_csv(supp_table1, file.path(OUTDIR, "Supplementary_Table1_discordant_genes.csv"), na = "")
}

# ----------------------------
# 4) Auto-fill QC columns for Table 1 from per-variant file (VARFILE)
#    NOTE: This is only used for Table 1 QC columns, NOT for gene_features.rds
# ----------------------------
stopifnot(file.exists(VARFILE))
v0 <- read_csv(VARFILE, show_col_types = FALSE)

has_flags      <- "flags" %in% names(v0)
has_variant_id <- "variant_id" %in% names(v0)

need_cols <- c("gene_symbol","consequence","chrom","pos","ref","alt")
missing_need <- setdiff(need_cols, names(v0))
if (length(missing_need) > 0) {
  stop("gnomad_variants_all.csv is missing required columns: ",
       paste(missing_need, collapse = ", "))
}

v_all <- v0 %>%
  mutate(
    gene_u = toupper(trimws(gene_symbol)),
    consequence_l = tolower(consequence %||% ""),
    
    # pLoF definition for QC purposes
    is_stop = str_detect(consequence_l, "stop_gained"),
    is_fs   = str_detect(consequence_l, "frameshift_variant"),
    is_plof = is_stop | is_fs,
    
    # segdup/LCR from flags (if present)
    flags_l = if (has_flags) tolower(flags %||% "") else "",
    in_segdup = if (has_flags) str_detect(flags_l, "segdup") else FALSE,
    in_lcr    = if (has_flags) str_detect(flags_l, "\\blcr\\b") else FALSE,
    in_rep    = in_segdup | in_lcr,
    
    # stable variant id
    var_id = if (has_variant_id) {
      as.character(variant_id)
    } else {
      paste0(gsub("^chr","", as.character(chrom)), ":", as.character(pos), ":", ref, ">", alt)
    }
  ) %>%
  filter(!is.na(gene_u), gene_u != "", !is.na(var_id), var_id != "")

v_plof <- v_all %>% filter(is_plof)

gene_qc <- v_plof %>%
  group_by(gene_u) %>%
  summarise(
    Observed_pLoF_classes = paste(
      c(if (any(is_stop, na.rm = TRUE)) "stop" else NULL,
        if (any(is_fs,   na.rm = TRUE)) "frameshift" else NULL),
      collapse = ", "
    ),
    n_distinct_pLoF = dplyr::n_distinct(var_id),
    `Multiple independent pLoF variants` = ifelse(n_distinct_pLoF >= 2, "Yes", "No"),
    frac_in_rep = mean(in_rep, na.rm = TRUE),
    `segdup/LCR flags` = case_when(
      is.nan(frac_in_rep) ~ "",
      frac_in_rep == 0    ~ "None",
      frac_in_rep < 0.5   ~ "Some",
      TRUE                ~ "Predominant"
    ),
    .groups = "drop"
  )

# ----------------------------
# 5) Table 1 (Word-ready top N) + write CSV + DOCX
# ----------------------------
table1_word_ready <- discordant_genelevel %>%
  slice_head(n = TOP_N) %>%
  mutate(gene_u = toupper(trimws(gene))) %>%
  left_join(gene_qc, by = "gene_u") %>%
  transmute(
    Gene = gene,
    LOEUF = round(loeuf, 3),
    `LoF o/e` = round(oe_lof, 3),
    `90% CI (LoF o/e)` = paste0(
      "(", sprintf("%.3f", oe_lof_lower), "–", sprintf("%.3f", oe_lof_upper), ")"
    ),
    `Observed pLoF classes` = coalesce(Observed_pLoF_classes, ""),
    `Multiple independent pLoF variants` = coalesce(`Multiple independent pLoF variants`, ""),
    `segdup/LCR flags` = coalesce(`segdup/LCR flags`, ""),
    `gnomAD constraint flags` = case_when(
      is.na(constraint_flags) ~ "None",
      str_trim(constraint_flags) == "" ~ "None",
      constraint_flags %in% c("[]", "[ ]") ~ "None",
      TRUE ~ "Present"
    )
  )

# helpful warnings if the QC join didn't populate anything
if (all(table1_word_ready$`Observed pLoF classes` == "")) {
  warning("Observed pLoF classes are empty for all rows. Check that consequence contains stop_gained/frameshift_variant in gnomad_variants_all.csv.")
}
if (all(table1_word_ready$`Multiple independent pLoF variants` == "")) {
  warning("Multiple independent pLoF variants empty for all rows. Check var_id construction / presence of multiple distinct variants per gene.")
}

OUT_TABLE1_CSV  <- file.path(TABDIR, "Table1_word_ready_top15.csv")
OUT_TABLE1_DOCX <- file.path(TABDIR, "Table1_word_ready_top15.docx")

write_csv(table1_word_ready, OUT_TABLE1_CSV, na = "")

ft <- flextable(table1_word_ready) %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Table 1. Genes with discordant LoF constraint (LOEUF < 0.2)", style = "heading 2") %>%
  body_add_flextable(ft)

print(doc, target = OUT_TABLE1_DOCX)

message("[DONE] Wrote:\n  - ", OUT_TABLE1_CSV,
        "\n  - ", OUT_TABLE1_DOCX,
        "\n  - ", file.path(TABDIR, "Supplementary_Table1_discordant_genes.csv"),
        "\n  - ", file.path(CACHEDIR, "gene_features.rds"))
