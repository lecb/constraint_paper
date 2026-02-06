#!/usr/bin/env Rscript

# ============================================================
# Domain-aware truncation contrast:
#   Discordant genes vs matched non-discordant genes
#
# Uses gnomad_variants_all.csv + InterPro domains.
# IMPORTANT: protein-domain analysis requires AA positions,
# so we restrict to truncating consequences that reliably
# have HGVS protein strings (stop_gained, frameshift_variant)
# and require non-missing hgvsp.
#
# Outputs:
#   matched_nondiscordant_genes.csv
#   domain_contrast_variant_level_QC.csv
#   SuppTable_domain_contrast_discordant_vs_matched_nondiscordant.csv
#   domain_contrast_tests.txt
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

# ---------------------------
# User-set inputs
# ---------------------------
DISCORDANT_RDS <- "discordant_genes.rds"
VARIANTS_ALL   <- "gnomad_variants_all.csv"

# Domain file (edit if your filename differs)
DOMAIN_PATH <- {
  candidates <- c(
    "domains_interpro_biomart.csv",
    "domains_interpro.csv",
    "interpro_domains.csv",
    "biomart_interpro_domains.csv"
  )
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) stop("No InterPro domain CSV found. Set DOMAIN_PATH manually.")
  candidates[[1]]
}

# Optional constraint metrics file for LOEUF matching
CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"

# Optional protein length cache (if available; otherwise matching falls back to LOEUF or random)
LEN_CACHE_RDS <- "ensembl_translation_length_cache.rds"

# Matching: number of non-discordant genes per discordant gene
K_MATCH <- 3

# Minimum number of truncations per gene to include in gene-level tests
MIN_TRUNC_PER_GENE <- 5

# Truncating classes to include for domain-coordinate analysis
TRUNC_CSQ_KEEP <- c("stop_gained", "frameshift_variant")

# Outputs
OUT_MATCHED    <- "matched_nondiscordant_genes.csv"
OUT_VAR_QC     <- "domain_contrast_variant_level_QC.csv"
OUT_GENE_TABLE <- "SuppTable_domain_contrast_discordant_vs_matched_nondiscordant.csv"
OUT_TESTS      <- "domain_contrast_tests.txt"

# ---------------------------
# Helpers
# ---------------------------
norm_gene <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- str_replace_all(x, "[^A-Z0-9]", "")
  x
}

parse_aa_pos_from_hgvsp <- function(hgvsp) {
  # Parse the first AA position in HGVS protein strings.
  # Examples:
  #   p.Glu24Ter         -> 24
  #   p.Ala11GlyfsTer61  -> 11 (frameshift start position)
  suppressWarnings(as.integer(str_match(hgvsp, "^p\\.[A-Za-z\\*]+(\\d+)")[,2]))
}

# ---------------------------
# 1) Load discordant genes (HGNC symbols in tibble$Gene)
# ---------------------------
dg <- readRDS(DISCORDANT_RDS)
stopifnot(is.data.frame(dg), "Gene" %in% names(dg))

discordant <- unique(norm_gene(dg$Gene))
message("[INFO] discordant genes loaded: ", length(discordant))
message("[INFO] first few discordant genes: ", paste(head(discordant, 10), collapse = ", "))

# ---------------------------
# 2) Load variants + set explicit column names (your schema)
# ---------------------------
va <- readr::read_csv(VARIANTS_ALL, show_col_types = FALSE)

# From your header:
gene_col  <- "gene_symbol"
hgvsp_col <- "hgvsp"
cons_col  <- "consequence"

stopifnot(all(c(gene_col, hgvsp_col, cons_col) %in% names(va)))

# ---------------------------
# 3) Filter to truncations with protein coordinates (stop/frameshift + hgvsp)
# ---------------------------
va2 <- va %>%
  transmute(
    gene_raw = as.character(.data[[gene_col]]),
    gene     = norm_gene(.data[[gene_col]]),
    hgvsp    = as.character(.data[[hgvsp_col]]),
    consequence = tolower(as.character(.data[[cons_col]]))
  ) %>%
  filter(
    !is.na(gene), gene != "",
    !is.na(hgvsp), hgvsp != "",
    consequence %in% TRUNC_CSQ_KEEP
  ) %>%
  mutate(
    aa_pos = parse_aa_pos_from_hgvsp(hgvsp)
  ) %>%
  filter(!is.na(aa_pos), aa_pos > 0)

message("[INFO] stop/frameshift with hgvsp+aa_pos retained: ", nrow(va2),
        " | genes: ", n_distinct(va2$gene))

# QC: how many stop/frameshift rows have hgvsp?
qc_summary <- va %>%
  transmute(
    consequence = tolower(as.character(.data[[cons_col]])),
    hgvsp = as.character(.data[[hgvsp_col]])
  ) %>%
  filter(consequence %in% TRUNC_CSQ_KEEP) %>%
  summarise(
    n_total_trunc = n(),
    n_hgvsp_present = sum(!is.na(hgvsp) & hgvsp != ""),
    frac_hgvsp_present = n_hgvsp_present / n_total_trunc
  )
message("[QC] stop/frameshift with hgvsp present fraction:")
print(qc_summary)

# ---- Prefer core positional dataset if available (has prot_len_aa + group) ----
PREFERRED_POS <- "trunc_position_domain_nmd_variants.csv"
if (file.exists(PREFERRED_POS)) {
  message("[INFO] Overriding va2 using preferred positional input: ", PREFERRED_POS)
  va2 <- readr::read_csv(PREFERRED_POS, show_col_types = FALSE) %>%
    dplyr::transmute(
      gene_raw = as.character(.data[["Gene"]]),
      gene     = norm_gene(.data[["Gene_norm"]]),
      aa_pos   = as.integer(.data[["aa_pos"]]),
      prot_len_aa = as.integer(.data[["prot_len_aa"]]),
      rel_pos  = as.numeric(.data[["rel_pos"]]),
      group    = as.character(.data[["group"]]),
      consequence = tolower(as.character(.data[["consequence"]]))
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(aa_pos), aa_pos > 0)

  message("[INFO] Preferred va2 rows: ", nrow(va2), " | genes: ", dplyr::n_distinct(va2$gene))

# Normalise group labels to exactly two groups: Discordant vs Background
if ("group" %in% names(va2)) {
  va2$group <- ifelse(va2$group %in% c("Background","background","Non-discordant","Non-discordant candidate","Non-discordant candidates"),
                     "Background", va2$group)
  va2$group <- ifelse(va2$group %in% c("Discordant","discordant"), "Discordant", va2$group)
}

  message("[INFO] prot_len_aa non-missing: ", sum(!is.na(va2$prot_len_aa)))
}

# ---------------------------
# 4) Protein length per gene (optional; used for matching)
# ---------------------------
gene_len <- tibble(gene = unique(va2$gene), prot_len_aa = NA_integer_)

# If prot_len_aa is already present (e.g. from trunc_position_domain_nmd_variants.csv), use it directly.
if ("prot_len_aa" %in% names(va2)) {
  gene_len <- va2 %>%
    transmute(gene = gene, prot_len_aa = as.integer(prot_len_aa)) %>%
    group_by(gene) %>%
    summarise(prot_len_aa = max(prot_len_aa, na.rm = TRUE), .groups = "drop")
}

if (file.exists(LEN_CACHE_RDS) && all(is.na(gene_len$prot_len_aa))) {
  lc <- readRDS(LEN_CACHE_RDS)
  
  if (!is.null(names(lc))) {
    # cache may be:
    #  - named numeric vector keyed by gene
    #  - named list keyed by protein_id (ENSP...) with numeric values
    if (is.list(lc) && !is.data.frame(lc)) {
      lc <- unlist(lc, use.names = TRUE)
    }

    if (is.numeric(lc)) {
      nms <- names(lc)

      # If names look like ENSP protein IDs, map ENSP -> gene using VEP map cache
      if (all(grepl("^ENSP", nms))) {
        MAP_RDS <- file.path("cache","vep_input_to_proteinpos_map.rds")
        if (!file.exists(MAP_RDS)) stop("Missing protein_id->gene map: ", MAP_RDS)

        pm <- readRDS(MAP_RDS)
        cn <- tolower(names(pm))
        pid_col <- names(pm)[match(TRUE, cn %in% c("protein_id","ensp","protein"))]
        gcol    <- names(pm)[match(TRUE, cn %in% c("gene","symbol","hgnc_symbol","gene_norm"))]
        if (is.na(pid_col) || is.na(gcol)) {
          stop("Could not find protein_id/gene columns in ", MAP_RDS, " | cols: ", paste(names(pm), collapse=", "))
        }

        pm2 <- pm %>%
          transmute(
            protein_id = as.character(.data[[pid_col]]),
            gene = norm_gene(.data[[gcol]])
          ) %>%
          filter(!is.na(protein_id), protein_id != "", !is.na(gene), gene != "") %>%
          distinct()

        len_pid <- tibble(protein_id = as.character(names(lc)),
                          prot_len_aa = as.integer(lc)) %>% distinct()

        tmp <- pm2 %>%
          left_join(len_pid, by = "protein_id") %>%
          filter(!is.na(prot_len_aa)) %>%
          group_by(gene) %>%
          summarise(prot_len_aa = max(prot_len_aa, na.rm = TRUE), .groups = "drop")

      } else {
        # Otherwise assume names are gene symbols
        tmp <- tibble(gene = norm_gene(names(lc)),
                      prot_len_aa = as.integer(lc)) %>% distinct()
      }

      gene_len <- gene_len %>%
        left_join(tmp, by = "gene") %>%
        transmute(gene, prot_len_aa = coalesce(prot_len_aa.y, prot_len_aa.x))
    }
  } else if (is.data.frame(lc)) {
    cn <- tolower(names(lc))
    gcol <- names(lc)[match(TRUE, cn %in% c("gene","hgnc_symbol","symbol","gene_norm"))]
    lcol <- names(lc)[match(TRUE, cn %in% c("prot_len_aa","protein_length","len_aa","translation_length"))]
    if (!is.na(gcol) && !is.na(lcol)) {
      tmp <- tibble(gene = norm_gene(lc[[gcol]]), prot_len_aa = as.integer(lc[[lcol]])) %>% distinct()
      gene_len <- gene_len %>%
        left_join(tmp, by = "gene") %>%
        transmute(gene, prot_len_aa = coalesce(prot_len_aa.y, prot_len_aa.x))
    }
  }
}

message("[INFO] genes with non-missing prot_len_aa: ", sum(!is.na(gene_len$prot_len_aa)))

# ---------------------------
# 5) Optional LOEUF for matching
# ---------------------------
loeuf_tbl <- NULL
if (file.exists(CONSTRAINT_TSV)) {
  ct <- readr::read_tsv(CONSTRAINT_TSV, show_col_types = FALSE)
  gcol <- names(ct)[match(TRUE, tolower(names(ct)) %in% c("gene","hgnc_symbol","symbol"))]
  
  # Try a few plausible loeuf columns; adjust if yours differs.
  lcol <- names(ct)[match(TRUE, tolower(names(ct)) %in% c("loeuf","lof_oe_ub","lof_oe","lof_oe_ci_upper"))]
  
  if (!is.na(gcol) && !is.na(lcol)) {
    loeuf_tbl <- ct %>%
      transmute(
        gene = norm_gene(.data[[gcol]]),
        loeuf = suppressWarnings(as.numeric(.data[[lcol]]))
      ) %>%
      filter(!is.na(gene), gene != "") %>%
      distinct()
    message("[INFO] constraint table loaded; loeuf rows: ", nrow(loeuf_tbl),
            " | using column: ", lcol)
  } else {
    message("[WARN] constraint table present but could not find gene/loeuf columns.")
  }
} else {
  message("[INFO] constraint table not found; matching will use protein length only (if available).")
}

# ---------------------------
# 6) Define analyzable universe and match non-discordant genes
# ---------------------------
genes_all <- unique(va2$gene)

disc_genes <- intersect(discordant, genes_all)
if (length(disc_genes) == 0) stop("No discordant genes found among analyzable genes in gnomad_variants_all.csv")

non_disc_candidates <- setdiff(genes_all, discordant)

gene_meta <- tibble(gene = genes_all) %>%
  left_join(gene_len, by = "gene")

if (!is.null(loeuf_tbl)) gene_meta <- gene_meta %>% left_join(loeuf_tbl, by = "gene")

disc_meta <- gene_meta %>% filter(gene %in% disc_genes)
cand_meta <- gene_meta %>% filter(gene %in% non_disc_candidates)

message("[INFO] analyzable discordant genes: ", nrow(disc_meta))
message("[INFO] analyzable non-discordant candidate genes: ", nrow(cand_meta))

match_one_gene <- function(target_gene, k = 3) {
  trow <- disc_meta %>% filter(gene == target_gene)
  if (nrow(trow) != 1) return(character(0))
  
  cand <- cand_meta
  
  # scoring:
  # - primary: abs protein length diff (requires prot_len_aa)
  # - secondary: loeuf diff if available (lower weight by scaling)
  cand <- cand %>%
    mutate(
      d_len = if (!is.na(trow) & !is.na(prot_len_aa))
        abs(prot_len_aa - trow$prot_len_aa) else NA_real_,
      d_loeuf = if (!is.null(loeuf_tbl) && !is.na(trow$loeuf) && !is.na(loeuf))
        abs(loeuf - trow$loeuf) else 0,
      score = coalesce(d_len, 1e12) + 1e4 * d_loeuf
    ) %>%
    arrange(score)
  
  head(cand$gene, k)
}

matched_non_disc <- base::setdiff(unique(va2), disc_genes)  # pipeline: use all background genes (no matching)

# ---- Pipeline mode: use ALL background genes (not matched) for S9 ----
FORCE_ALL_BACKGROUND <- TRUE
if (isTRUE(FORCE_ALL_BACKGROUND)) {
  matched_non_disc <- base::setdiff(unique(va2$gene), disc_genes)
  message("[INFO] FORCE_ALL_BACKGROUND=TRUE | using background genes: ", length(matched_non_disc))
}

write_csv(tibble(gene = matched_non_disc), OUT_MATCHED)
message("[WROTE] matched non-discordant gene list: ", OUT_MATCHED, " (n=", length(matched_non_disc), ")")

genes_use <- unique(c(disc_genes, matched_non_disc))
message("[INFO] total genes to analyse (discordant + matched): ", length(genes_use))

# ---------------------------
# 7) Load InterPro domains + compute domain metrics
# ---------------------------
domains_raw <- readr::read_csv(DOMAIN_PATH, show_col_types = FALSE)
nm <- names(domains_raw); nm_l <- tolower(nm)

col_gene  <- nm[match(TRUE, nm_l %in% c("gene","hgnc_symbol","symbol","hgnc"))]
col_start <- nm[match(TRUE, nm_l %in% c("interpro_start","start","dom_start"))]
col_end   <- nm[match(TRUE, nm_l %in% c("interpro_end","end","dom_end"))]
col_id    <- nm[match(TRUE, nm_l %in% c("interpro","domain_id","dom_id"))]

if (any(is.na(c(col_gene, col_start, col_end)))) {
  stop("Domain file missing gene/start/end columns. Found: ", paste(nm, collapse=", "))
}

domains <- domains_raw %>%
  transmute(
    gene = norm_gene(.data[[col_gene]]),
    dom_id = if (!is.na(col_id)) as.character(.data[[col_id]]) else NA_character_,
    dom_start = suppressWarnings(as.integer(.data[[col_start]])),
    dom_end   = suppressWarnings(as.integer(.data[[col_end]]))
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(dom_start), !is.na(dom_end)) %>%
  mutate(
    dom_lo = pmin(dom_start, dom_end),
    dom_hi = pmax(dom_start, dom_end)
  ) %>%
  distinct()

dom_gene <- domains %>%
  group_by(gene) %>%
  summarise(
    n_domains = n_distinct(dom_id, na.rm = TRUE),
    last_domain_end_aa = max(dom_hi, na.rm = TRUE),
    .groups = "drop"
  )

dom_by_gene <- split(domains, domains$gene)

calc_metrics <- function(g, pos) {
  d <- dom_by_gene[[g]]
  if (is.null(d) || nrow(d) == 0) {
    return(list(dist=NA_integer_, within=NA, after_last=NA))
  }
  dist <- min(abs(pos - d$dom_lo), abs(pos - d$dom_hi), na.rm = TRUE)
  within <- any(pos >= d$dom_lo & pos <= d$dom_hi, na.rm = TRUE)
  after_last <- pos > max(d$dom_hi, na.rm = TRUE)
  list(dist = as.integer(dist), within = within, after_last = after_last)
}

# Variant-level subset to analyse
va_use <- va2 %>%
  filter(gene %in% genes_use) %>%
  mutate(group = if_else(gene %in% disc_genes, "Discordant", "Non-discordant")) %>%
  left_join(dom_gene, by = "gene") %>%
  rowwise() %>%
  mutate(
    tmp = list(calc_metrics(gene, aa_pos)),
    dist_to_boundary_aa = tmp$dist,
    within_any_domain = tmp$within,
    after_last_domain = tmp$after_last
  ) %>%
  ungroup() %>%
  select(-tmp)

write_csv(va_use, OUT_VAR_QC)
message("[WROTE] variant-level QC: ", OUT_VAR_QC)

# ---------------------------
# 8) Gene-level summaries and tests
# ---------------------------
gene_sum <- va_use %>%
  group_by(gene, group) %>%
  summarise(
    n_trunc = n(),
    n_domains = suppressWarnings(max(n_domains, na.rm = TRUE)),
    prop_within_domain = mean(within_any_domain, na.rm = TRUE),
    median_dist_to_boundary_aa = median(dist_to_boundary_aa, na.rm = TRUE),
    prop_after_last_domain = mean(after_last_domain, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(group == "Discordant"), gene)

write_csv(gene_sum, OUT_GENE_TABLE)
write_csv(gene_sum, file.path("tables","SuppTable_domain_structured_truncation_features.csv"), na = "")
message("[WROTE] gene-level table: ", OUT_GENE_TABLE)

# Filter for tests (avoid tiny-n genes dominating)
gene_sum_filt <- gene_sum %>% filter(n_trunc >= MIN_TRUNC_PER_GENE)

# Wilcoxon tests across genes
w1 <- tryCatch(wilcox.test(prop_within_domain ~ group, data = gene_sum_filt), error=function(e)e)
w2 <- tryCatch(wilcox.test(median_dist_to_boundary_aa ~ group, data = gene_sum_filt), error=function(e)e)

sink(OUT_TESTS)
cat("Domain-aware truncation contrast: discordant vs matched non-discordant\n")
cat("=====================================================================\n\n")
cat("Truncation classes included: ", paste(TRUNC_CSQ_KEEP, collapse = ", "), "\n\n", sep = "")
cat("Minimum truncations per gene for tests: ", MIN_TRUNC_PER_GENE, "\n\n", sep = "")

cat("Genes per group (gene_sum):\n")
print(gene_sum %>% distinct(gene, group) %>% count(group, name = "n_genes"))
cat("\nGenes per group (gene_sum_filt):\n")
print(gene_sum_filt %>% distinct(gene, group) %>% count(group, name = "n_genes"))
cat("\n")

cat("Test 1: prop_within_domain ~ group\n")
if (inherits(w1,"error")) cat("ERROR: ", w1$message, "\n\n") else { print(w1); cat("\n\n") }

cat("Test 2: median_dist_to_boundary_aa ~ group\n")
if (inherits(w2,"error")) cat("ERROR: ", w2$message, "\n\n") else { print(w2); cat("\n\n") }

cat("Summary medians (gene_sum_filt):\n")
print(gene_sum_filt %>% group_by(group) %>% summarise(
  median_prop_within = median(prop_within_domain, na.rm=TRUE),
  median_dist_to_boundary = median(median_dist_to_boundary_aa, na.rm=TRUE),
  median_prop_after_last = median(prop_after_last_domain, na.rm=TRUE),
  .groups = "drop"
))
sink()
message("[WROTE] tests: ", OUT_TESTS)

# Console summaries
message("\n[SUMMARY] medians by group (gene_sum_filt):")
print(gene_sum_filt %>% group_by(group) %>% summarise(
  median_prop_within = median(prop_within_domain, na.rm=TRUE),
  median_dist_to_boundary = median(median_dist_to_boundary_aa, na.rm=TRUE),
  median_prop_after_last = median(prop_after_last_domain, na.rm=TRUE),
  .groups = "drop"
))



# ---- Write canonical Supplementary Table S9 (pipeline output) ----
if (exists("gene_sum")) {
  dir.create("tables", showWarnings = FALSE, recursive = TRUE)
  out_s9 <- file.path("tables", "SuppTable_domain_structured_truncation_features.csv")
  
  # Sanitize domainless genes: avoid -Inf from max(..., na.rm=TRUE) on empty sets
  if (all(c("n_domains","median_dist_to_boundary_aa","prop_within_domain","prop_after_last_domain") %in% names(gene_sum))) {
    gene_sum <- gene_sum %>%
      dplyr::mutate(
        n_domains = dplyr::if_else(is.finite(n_domains) & n_domains >= 0, as.integer(n_domains), 0L),
        median_dist_to_boundary_aa = dplyr::if_else(is.finite(median_dist_to_boundary_aa), median_dist_to_boundary_aa, NA_real_),
        prop_within_domain = dplyr::if_else(is.finite(prop_within_domain), prop_within_domain, NA_real_),
        prop_after_last_domain = dplyr::if_else(is.finite(prop_after_last_domain), prop_after_last_domain, NA_real_)
      )
  }
readr::write_csv(gene_sum, out_s9, na = "")
  message("[WROTE] ", out_s9)
}

# ---- More informative summary: only genes with >=1 domain ----
if (exists("gene_sum") && all(c("group","n_domains","prop_within_domain","median_dist_to_boundary_aa","prop_after_last_domain") %in% names(gene_sum))) {
  message("\n[SUMMARY] medians by group among genes with >=1 domain (n_domains > 0):")
  tmp <- gene_sum %>%
    dplyr::filter(n_domains > 0) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      median_prop_within = stats::median(prop_within_domain, na.rm = TRUE),
      median_dist_to_boundary = stats::median(median_dist_to_boundary_aa, na.rm = TRUE),
      median_prop_after_last = stats::median(prop_after_last_domain, na.rm = TRUE),
      .groups = "drop"
    )
  print(tmp)
}

# ---- Robust cleanup: n_domains as safe integer without warnings ----
if (exists("gene_sum") && "n_domains" %in% names(gene_sum)) {
  nd <- suppressWarnings(as.numeric(gene_sum$n_domains))
  nd[!is.finite(nd) | nd < 0] <- 0
  nd <- pmin(nd, .Machine$integer.max)
  gene_sum$n_domains <- as.integer(nd)
}
