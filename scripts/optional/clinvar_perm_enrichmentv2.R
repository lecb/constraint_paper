# ============================================================
# FULL PIPELINE (amended + cache-friendly):
# ClinVar VCF → (discordant genes) → truncating variants → P/LP
# + ClinVar gold-star mapping + ≥2-star subset
# + exclude conflicting assertions
# + split Pathogenic vs Likely Pathogenic
# + per-gene enrichment vs expectation (binomial null weighted by gene-specific truncation opportunity)
#
# Outputs (TSV):
#  - discordant_genes_from_filenames.txt
#
#  - clinvar_trunc_discordant_overall.tsv
#  - clinvar_trunc_discordant_by_gene.tsv
#  - clinvar_trunc_discordant_variant_level.tsv
#
#  - clinvar_plp_trunc_discordant_variant_level.tsv
#  - clinvar_plp_trunc_discordant_overall.tsv
#  - clinvar_plp_trunc_discordant_by_gene.tsv
#  - clinvar_plp_trunc_discordant_by_reviewstatus.tsv
#  - clinvar_plp_trunc_discordant_by_stars.tsv
#  - clinvar_plp_trunc_discordant_by_significance.tsv
#
#  - clinvar_plp_trunc_discordant_ge2stars_variant_level.tsv
#  - clinvar_plp_trunc_discordant_ge2stars_overall.tsv
#  - clinvar_plp_trunc_discordant_ge2stars_by_gene.tsv
#
#  - clinvar_plp_trunc_discordant_ge2stars_enrichment_by_gene.tsv
#
#  - Supplementary_Table_ClinVar_PLP_truncating_discordant_genes.tsv
#  - Supplementary_Table_ClinVar_PLP_truncating_discordant_genes_GE2stars.tsv
#
# Notes:
#  - "Gold stars" are inferred from CLNREVSTAT using ClinVar conventions.
#  - Conflicting assertions are excluded by CLNSIG content (and optionally CLNREVSTAT).
#  - P vs LP is based on CLNSIG terms.
#  - Enrichment expectation: total qualifying variants distributed across genes
#    proportional to each gene’s *overall truncating ClinVar variant count*
#    (MC-defined), acting as an “opportunity” weight.
# KEEP FOR PAPER
# ============================================================

suppressPackageStartupMessages({
  library(VariantAnnotation)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(readr)
})

# ----------------------------
# CONFIG
# ----------------------------
DISCORDANT_DIR <- "/Users/eseaby/Documents/PhD/Project/constraint/gnomAD_genes"
VCF_PATH <- file.path("clinvar","clinvar.vcf.gz")
OUTDIR <- "clinvar_bridge"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
TABDIR <- "tables"
dir.create(TABDIR, showWarnings = FALSE, recursive = TRUE)

CACHE_VCF_RDS              <- file.path(OUTDIR, "cache_clinvar_vcf.rds")
CACHE_TRUNC_ALL_RDS        <- file.path(OUTDIR, "cache_vcf_trunc_all.rds")
CACHE_TRUNC_DISC_RDS       <- file.path(OUTDIR, "cache_vcf_trunc_discordant.rds")
CACHE_TRUNC_DISC_TBL_RDS   <- file.path(OUTDIR, "cache_trunc_discordant_variant_tbl.rds")

# Set TRUE to force rebuilding caches
FORCE <- FALSE

DISCORDANT35_PATH <- file.path("tables","Supplementary_Table1_discordant_genes.csv")
TOP_N_DISCORDANT <- 35L
# ----------------------------
# 1) Discordant gene list from filenames
# ----------------------------
# ----------------------------
# 1) Discordant gene list (TOP 35) from Supplementary Table 1
# ----------------------------
stopifnot(file.exists(DISCORDANT35_PATH))

disc_tbl <- suppressMessages(readr::read_csv(DISCORDANT35_PATH, show_col_types = FALSE))

# Find likely gene column
gene_col <- intersect(names(disc_tbl), c("gene", "Gene", "SYMBOL", "symbol"))
stopifnot(length(gene_col) >= 1)

discordant_genes <- disc_tbl[[gene_col[1]]] |>
  as.character() |>
  unique() |>
  sort()

# Enforce TOP_N_DISCORDANT
if (length(discordant_genes) < TOP_N_DISCORDANT) {
  stop("Discordant gene list has < ", TOP_N_DISCORDANT, " genes. Check file contents.")
}
discordant_genes <- discordant_genes[1:TOP_N_DISCORDANT]

message("[OK] discordant genes (top", TOP_N_DISCORDANT, "): ", length(discordant_genes))
writeLines(discordant_genes, file.path(OUTDIR, paste0("discordant_genes_top", TOP_N_DISCORDANT, ".txt")))

# ----------------------------
# 2) Read ClinVar VCF (cached)
# ----------------------------
stopifnot(file.exists(VCF_PATH))

if (!FORCE && file.exists(CACHE_VCF_RDS)) {
  message("[CACHE HIT] ClinVar VCF: ", CACHE_VCF_RDS)
  vcf <- readRDS(CACHE_VCF_RDS)
} else {
  message("[CACHE BUILD] Reading ClinVar VCF (heavy step) ...")
  vcf <- readVcf(VCF_PATH, genome = "hg38")
  saveRDS(vcf, CACHE_VCF_RDS)
  message("[OK] VCF cached: ", CACHE_VCF_RDS)
}
message("[OK] VCF loaded: ", nrow(vcf), " variants")

# ----------------------------
# 3) Helpers
# ----------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

extract_symbols <- function(geneinfo) {
  if (is.na(geneinfo) || geneinfo == "") return(character(0))
  parts <- unlist(strsplit(as.character(geneinfo), "\\|"))
  syms <- sub(":.*$", "", parts)
  syms <- syms[syms != ""]
  as.character(unique(syms))
}

collapse_list_field <- function(x) {
  # x is often a CharacterList (list-like) INFO field
  vapply(x, function(z) {
    if (length(z) == 0) return(NA_character_)
    paste(as.character(z), collapse = "|")
  }, character(1))
}

split_terms <- function(x) {
  # split on pipes/commas/semicolons; normalise to lower snake-ish
  if (is.na(x) || x == "") return(character(0))
  y <- tolower(as.character(x))
  y <- gsub("[\\s]+", "_", y)
  y <- gsub("[\\(\\)\\[\\]]", "", y)
  parts <- unlist(strsplit(y, "[\\|,;]"))
  parts <- parts[parts != ""]
  unique(parts)
}

# ---- ClinVar CLNSIG parsing (P vs LP vs conflict) ----
classify_clnsig <- function(clnsig_str) {
  terms <- split_terms(clnsig_str)
  
  # Common conflict indicator
  has_conflict_term <- any(grepl("conflict", terms)) ||
    any(grepl("conflicting_interpretations", terms)) ||
    any(terms %in% c("conflicting_interpretations_of_pathogenicity"))
  
  has_path <- any(terms %in% c("pathogenic"))
  has_lp   <- any(terms %in% c("likely_pathogenic", "likely-pathogenic", "likelypathogenic"))
  has_ben  <- any(terms %in% c("benign", "likely_benign", "likely-benign", "likelybenign"))
  
  # If both pathogenic-ish and benign-ish appear, treat as conflict even if "conflict" not explicit
  has_pathlike <- has_path || has_lp
  has_benignlike <- has_ben || any(terms %in% c("uncertain_significance", "vus"))
  
  conflict <- has_conflict_term || (has_pathlike && has_ben)
  
  sig_class <- dplyr::case_when(
    conflict ~ "Conflicting",
    has_path && has_lp ~ "Pathogenic_or_Likely_pathogenic",
    has_path ~ "Pathogenic",
    has_lp ~ "Likely_pathogenic",
    TRUE ~ "Other_or_non_PLP"
  )
  
  list(
    sig_class = sig_class,
    is_plp = sig_class %in% c("Pathogenic", "Likely_pathogenic", "Pathogenic_or_Likely_pathogenic"),
    is_conflicting = conflict,
    is_pathogenic = sig_class %in% c("Pathogenic", "Pathogenic_or_Likely_pathogenic"),
    is_likely_pathogenic = sig_class %in% c("Likely_pathogenic", "Pathogenic_or_Likely_pathogenic")
  )
}

# ---- Truncating subset using ClinVar MC (VEP SO terms embedded) ----
is_trunc_mc_field <- function(mc_field) {
  if (length(mc_field) == 0) return(FALSE)
  entries <- as.character(mc_field)
  terms <- sub("^.*\\|", "", entries)
  terms <- tolower(terms)
  
  any(terms %in% c(
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "exon_loss_variant"
    # intentionally excluding: "start_lost", "stop_lost"
  ))
}


# ---- Gold-star mapping from CLNREVSTAT ----
# ClinVar convention (approximate but standard):
# 4 stars: practice guideline
# 3 stars: reviewed by expert panel
# 2 stars: criteria provided, multiple submitters, no conflicts
# 1 star : criteria provided, single submitter  (or multiple submitters without the "no conflicts" qualifier)
# 0 stars: no assertion criteria / no interpretation / other low evidence states

revstat_to_stars <- function(clnrevstat_str) {
  if (is.na(clnrevstat_str) || clnrevstat_str == "") return(0L)
  
  t <- tolower(as.character(clnrevstat_str))
  
  # Your VCF has "|_single_submitter" etc.
  # Convert "|_" to "|" so patterns match.
  t <- gsub("\\|_", "|", t)
  
  # Normalise whitespace/punctuation
  t <- gsub("[[:space:]]+", "_", t)
  t <- gsub("[^a-z0-9\\|_]", "", t)
  
  # 4 stars
  if (grepl("practice_guideline", t)) return(4L)
  
  # 3 stars
  if (grepl("reviewed_by_expert_panel", t)) return(3L)
  
  # 2 stars (must have no_conflicts)
  if (grepl("criteria_provided\\|multiple_submitters\\|no_conflicts", t)) return(2L)
  
  # 1 star
  if (grepl("criteria_provided\\|single_submitter", t)) return(1L)
  if (grepl("criteria_provided\\|multiple_submitters", t)) return(1L)
  
  # 0 stars
  0L
}


# ----------------------------
# 4) Build global "truncating" VCF (cached)
# ----------------------------
if (!FORCE && file.exists(CACHE_TRUNC_ALL_RDS)) {
  message("[CACHE HIT] vcf_trunc_all: ", CACHE_TRUNC_ALL_RDS)
  vcf_trunc_all <- readRDS(CACHE_TRUNC_ALL_RDS)
} else {
  message("[CACHE BUILD] Truncating subset over ALL ClinVar using MC ...")
  mc_all <- info(vcf)$MC
  is_trunc_all <- vapply(mc_all, is_trunc_mc_field, logical(1))
  vcf_trunc_all <- vcf[is_trunc_all]
  message("[OK] truncating (MC-defined) across all genes: ", length(vcf_trunc_all))
  saveRDS(vcf_trunc_all, CACHE_TRUNC_ALL_RDS)
  message("[OK] cached: ", CACHE_TRUNC_ALL_RDS)
}

# ----------------------------
# 5) Truncating subset restricted to discordant genes (cached)
# ----------------------------
if (!FORCE && file.exists(CACHE_TRUNC_DISC_RDS)) {
  message("[CACHE HIT] vcf_trunc_disc: ", CACHE_TRUNC_DISC_RDS)
  vcf_trunc_disc <- readRDS(CACHE_TRUNC_DISC_RDS)
} else {
  message("[CACHE BUILD] Restrict truncating to discordant genes using GENEINFO ...")
  geneinfo_all <- info(vcf_trunc_all)$GENEINFO
  in_discordant <- vapply(geneinfo_all, function(g) {
    any(extract_symbols(g) %in% discordant_genes)
  }, logical(1))
  vcf_trunc_disc <- vcf_trunc_all[in_discordant]
  message("[OK] truncating variants with any discordant gene in GENEINFO: ", length(vcf_trunc_disc))
  saveRDS(vcf_trunc_disc, CACHE_TRUNC_DISC_RDS)
  message("[OK] cached: ", CACHE_TRUNC_DISC_RDS)
}

# ----------------------------
# 6) Build a variant-level table for truncating discordant set (cached)
# ----------------------------
if (!FORCE && file.exists(CACHE_TRUNC_DISC_TBL_RDS)) {
  message("[CACHE HIT] trunc_discordant_variant_tbl: ", CACHE_TRUNC_DISC_TBL_RDS)
  trunc_variant_tbl <- readRDS(CACHE_TRUNC_DISC_TBL_RDS)
} else {
  message("[CACHE BUILD] Building trunc_discordant variant-level table ...")
  
  trunc_variant_tbl <- tibble(
    chrom = as.character(seqnames(rowRanges(vcf_trunc_disc))),
    pos   = start(rowRanges(vcf_trunc_disc)),
    ref   = as.character(rowRanges(vcf_trunc_disc)$REF),
    alt   = vapply(rowRanges(vcf_trunc_disc)$ALT,
                   function(a) paste(as.character(a), collapse=","), character(1)),
    CLNSIG     = collapse_list_field(info(vcf_trunc_disc)$CLNSIG),
    CLNREVSTAT = collapse_list_field(info(vcf_trunc_disc)$CLNREVSTAT),
    GENEINFO   = as.character(info(vcf_trunc_disc)$GENEINFO),
    MC         = vapply(info(vcf_trunc_disc)$MC,
                        function(z) paste(as.character(z), collapse=","), character(1)),
    CLNVC      = as.character(info(vcf_trunc_disc)$CLNVC)
  ) %>%
    mutate(
      stars = vapply(CLNREVSTAT, revstat_to_stars, integer(1))
    )
  
  saveRDS(trunc_variant_tbl, CACHE_TRUNC_DISC_TBL_RDS)
  message("[OK] cached: ", CACHE_TRUNC_DISC_TBL_RDS)
}

print(trunc_variant_tbl %>% count(CLNREVSTAT, stars, sort = TRUE))

# ----------------------------
# 7) Expand to gene-level rows (discordant genes only)
# ----------------------------
trunc_gene_df <- trunc_variant_tbl %>%
  mutate(gene = lapply(GENEINFO, extract_symbols)) %>%
  unnest(gene) %>%
  mutate(gene = as.character(gene)) %>%
  filter(gene %in% discordant_genes)

# ----------------------------
# 8) TRUNCATING summaries (all truncating, regardless of ClinVar significance)
# ----------------------------
trunc_overall <- tibble(
  trunc_variants = nrow(trunc_variant_tbl),
  trunc_genes = n_distinct(trunc_gene_df$gene)
)

trunc_by_gene <- trunc_gene_df %>%
  count(gene, sort = TRUE, name = "trunc_variants")

# ----------------------------
# 9) P/LP calls + conflict exclusion + Pathogenic vs Likely Pathogenic
# ----------------------------
sig_parsed <- lapply(trunc_variant_tbl$CLNSIG, classify_clnsig)
sig_df <- tibble(
  sig_class = vapply(sig_parsed, `[[`, character(1), "sig_class"),
  is_plp = vapply(sig_parsed, `[[`, logical(1), "is_plp"),
  is_conflicting = vapply(sig_parsed, `[[`, logical(1), "is_conflicting"),
  is_pathogenic = vapply(sig_parsed, `[[`, logical(1), "is_pathogenic"),
  is_likely_pathogenic = vapply(sig_parsed, `[[`, logical(1), "is_likely_pathogenic")
)

plp_tbl_all <- bind_cols(trunc_variant_tbl, sig_df) %>%
  # Exclude conflicts explicitly
  filter(is_plp, !is_conflicting)

# Split “Pathogenic vs Likely Pathogenic” into a single label for counting
plp_tbl_all <- plp_tbl_all %>%
  mutate(
    plp_label = case_when(
      is_pathogenic & is_likely_pathogenic ~ "Pathogenic_or_Likely_pathogenic",
      is_pathogenic ~ "Pathogenic",
      is_likely_pathogenic ~ "Likely_pathogenic",
      TRUE ~ "Other_or_non_PLP"
    )
  )

# Gene-expanded P/LP table (still discordant only)
plp_gene_df <- plp_tbl_all %>%
  mutate(gene = lapply(GENEINFO, extract_symbols)) %>%
  unnest(gene) %>%
  mutate(gene = as.character(gene)) %>%
  filter(gene %in% discordant_genes)

plp_overall <- tibble(
  plp_trunc_variants = nrow(plp_tbl_all),
  plp_trunc_genes = n_distinct(plp_gene_df$gene)
)

plp_by_gene <- plp_gene_df %>% count(gene, sort = TRUE, name = "plp_trunc_variants")

plp_by_review <- plp_tbl_all %>% count(CLNREVSTAT, sort = TRUE, name = "n")
plp_by_stars  <- plp_tbl_all %>% count(stars, sort = TRUE, name = "n")
plp_by_sig    <- plp_tbl_all %>% count(plp_label, sort = TRUE, name = "n")

# ----------------------------
# 10) ≥2-star subset (gold stars inferred from CLNREVSTAT)
# ----------------------------
plp_tbl_ge2 <- plp_tbl_all %>% filter(stars >= 2L)

plp_gene_df_ge2 <- plp_tbl_ge2 %>%
  mutate(gene = lapply(GENEINFO, extract_symbols)) %>%
  unnest(gene) %>%
  mutate(gene = as.character(gene)) %>%
  filter(gene %in% discordant_genes)

plp_overall_ge2 <- tibble(
  plp_trunc_variants_ge2stars = nrow(plp_tbl_ge2),
  plp_trunc_genes_ge2stars = n_distinct(plp_gene_df_ge2$gene)
)

plp_by_gene_ge2 <- plp_gene_df_ge2 %>%
  count(gene, sort = TRUE, name = "plp_trunc_variants_ge2stars")

# ----------------------------
# 11) Per-gene enrichment vs expectation (≥2-star P/LP truncating)
#     Expectation uses gene “opportunity” proportional to ALL ClinVar truncating variants (MC-defined)
#     in that gene (within discordant set).
# ----------------------------
# Opportunity: truncating counts per gene (discordant genes) from *all truncating* (not filtered by CLNSIG)
# ----------------------------
# 11) Per-gene enrichment vs expectation (≥2-star P/LP truncating)
#     Expectation uses gene “opportunity” proportional to ALL ClinVar truncating variants (MC-defined)
# ----------------------------

opportunity <- trunc_by_gene %>%
  rename(trunc_opportunity = trunc_variants)

total_qual <- nrow(plp_tbl_ge2)

observed <- plp_by_gene_ge2 %>%
  rename(observed_qual = plp_trunc_variants_ge2stars)

enrichment <- opportunity %>%
  left_join(observed, by = "gene") %>%
  mutate(observed_qual = coalesce(observed_qual, 0L)) %>%
  mutate(
    opp_total = sum(trunc_opportunity),
    p_null = ifelse(opp_total > 0, trunc_opportunity / opp_total, NA_real_),
    expected_qual = total_qual * p_null,
    enrichment_obs_over_exp = ifelse(expected_qual > 0, observed_qual / expected_qual, NA_real_),
    p_value = ifelse(
      is.na(p_null) | total_qual == 0,
      NA_real_,
      stats::pbinom(observed_qual - 1L, size = total_qual, prob = p_null, lower.tail = FALSE)
    )
  ) %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "BH")
  ) %>%
  select(
    gene,
    trunc_opportunity,
    observed_qual,
    expected_qual,
    enrichment_obs_over_exp,
    p_value,
    p_adj_fdr
  ) %>%
  arrange(p_adj_fdr, desc(enrichment_obs_over_exp))


# ----------------------------
# 12) Write outputs
# ----------------------------

# Trunc-only outputs (discordant, MC-defined)
fwrite(trunc_overall, file.path(OUTDIR, "clinvar_trunc_discordant_overall.tsv"), sep = "\t")
fwrite(trunc_by_gene, file.path(OUTDIR, "clinvar_trunc_discordant_by_gene.tsv"), sep = "\t")
fwrite(trunc_variant_tbl, file.path(OUTDIR, "clinvar_trunc_discordant_variant_level.tsv"), sep = "\t")

# P/LP trunc outputs (discordant, conflict-excluded)
fwrite(plp_tbl_all, file.path(OUTDIR, "clinvar_plp_trunc_discordant_variant_level.tsv"), sep = "\t")
fwrite(plp_overall, file.path(OUTDIR, "clinvar_plp_trunc_discordant_overall.tsv"), sep = "\t")
fwrite(plp_by_gene, file.path(OUTDIR, "clinvar_plp_trunc_discordant_by_gene.tsv"), sep = "\t")
fwrite(plp_by_review, file.path(OUTDIR, "clinvar_plp_trunc_discordant_by_reviewstatus.tsv"), sep = "\t")
fwrite(plp_by_stars, file.path(OUTDIR, "clinvar_plp_trunc_discordant_by_stars.tsv"), sep = "\t")
fwrite(plp_by_sig, file.path(OUTDIR, "clinvar_plp_trunc_discordant_by_significance.tsv"), sep = "\t")

# ≥2-star subset outputs
fwrite(plp_tbl_ge2, file.path(OUTDIR, "clinvar_plp_trunc_discordant_ge2stars_variant_level.tsv"), sep = "\t")
fwrite(plp_overall_ge2, file.path(OUTDIR, "clinvar_plp_trunc_discordant_ge2stars_overall.tsv"), sep = "\t")
fwrite(plp_by_gene_ge2, file.path(OUTDIR, "clinvar_plp_trunc_discordant_ge2stars_by_gene.tsv"), sep = "\t")

# Enrichment table (≥2-star)
fwrite(enrichment, file.path(OUTDIR, "clinvar_plp_trunc_discordant_ge2stars_enrichment_by_gene.tsv"), sep = "\t")

# Supplementary tables (paper-facing)
readr::write_tsv(
  plp_tbl_all,
  file.path(TABDIR, "Supplementary_Table_ClinVar_PLP_truncating_discordant_genes.tsv")
)
readr::write_tsv(
  plp_tbl_ge2,
  file.path(TABDIR, "Supplementary_Table_ClinVar_PLP_truncating_discordant_genes_GE2stars.tsv")
)

message("[WROTE] ", OUTDIR)

message("\n=== TRUNC (all) discordant summary ===")
print(trunc_overall)
print(head(trunc_by_gene, 15))

message("\n=== P/LP TRUNC (conflict-excluded) summary ===")
print(plp_overall)
print(head(plp_by_gene, 15))
print(head(plp_by_review, 10))
print(plp_by_stars)
print(plp_by_sig)

message("\n=== P/LP TRUNC ≥2-stars summary ===")
print(plp_overall_ge2)
print(head(plp_by_gene_ge2, 15))

message("\n=== Enrichment (≥2-stars) top hits (lowest FDR) ===")
print(head(enrichment, 15))
