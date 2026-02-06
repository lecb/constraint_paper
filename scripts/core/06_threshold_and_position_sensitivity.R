#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(tibble)
})

# ============================================================
# Threshold sensitivity + continuous discordance score
# OPTION A: Residualise trunc burden using ONLY syn_obs
#   (avoids cds_length missingness; keeps nearly all genes)
# Writes:
#   - gnomad_lof_discordance_out/threshold_sensitivity_summary.csv
#   - gnomad_lof_discordance_out/sensitivity_counts.csv
#   - gnomad_lof_discordance_out/discordance_score_summary.csv
#   - gnomad_lof_discordance_out/discordance_score_top50_genes.csv
# ============================================================


# ----------------------------
# CONFIG
# ----------------------------
TABLE_DIR <- "tables"
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)

# outputs go here
OUTDIR <- TABLE_DIR

# inputs come from the main pipeline cache (NOT tables/cache)
CACHEDIR <- file.path("gnomad_lof_discordance_out", "cache")
IN_RDS   <- file.path(CACHEDIR, "gene_features.rds")
stopifnot(file.exists(IN_RDS))


# ----------------------------
# Load gene features (one row per gene)
# ----------------------------
gene_df <- readRDS(IN_RDS)

# Required columns (NOTE: cds_length can be present but is not required for Option A)
need <- c("gene","loeuf","oe_lof","oe_lof_lower","oe_lof_upper","ci_width","syn_obs","trunc_obs")
missing <- setdiff(need, names(gene_df))
if (length(missing) > 0) stop("gene_features.rds missing columns: ", paste(missing, collapse = ", "))

# Ensure numeric where needed
gene_df <- gene_df %>%
  mutate(
    loeuf = as.numeric(loeuf),
    oe_lof = as.numeric(oe_lof),
    oe_lof_lower = as.numeric(oe_lof_lower),
    oe_lof_upper = as.numeric(oe_lof_upper),
    ci_width = as.numeric(ci_width),
    syn_obs = as.numeric(syn_obs),
    trunc_obs = as.numeric(trunc_obs)
  )

# ---------- helper: residualised trunc burden (OPTION A: syn only) ----------
add_trunc_resid <- function(df) {
  # Adjust for mutability using synonymous burden only (keeps almost all genes)
  fit <- lm(log1p(trunc_obs) ~ log1p(syn_obs), data = df)
  df %>% mutate(trunc_resid = resid(fit))
}

# ---------- classify discordant under a parameter set ----------
classify_paramset <- function(gene_df,
                              loeuf_cut = 0.2,
                              drop_ci_frac = 0.25,
                              top_oe_frac = 0.05) {
  
  # Stage 1: define universe using constraint metrics only
  u0 <- gene_df %>%
    filter(!is.na(loeuf), !is.na(oe_lof), !is.na(ci_width)) %>%
    filter(loeuf < loeuf_cut)
  
  if (nrow(u0) < 50) return(NULL)
  
  # Drop widest CI fraction within LOEUF subset
  ci_thr <- as.numeric(quantile(u0$ci_width, probs = 1 - drop_ci_frac, na.rm = TRUE))
  u1 <- u0 %>% filter(ci_width <= ci_thr)
  
  if (nrow(u1) < 30) return(NULL)
  
  # Discordant = top oe_lof fraction within remaining subset
  oe_thr <- as.numeric(quantile(u1$oe_lof, probs = 1 - top_oe_frac, na.rm = TRUE))
  u1 <- u1 %>% mutate(is_disc = oe_lof >= oe_thr)
  
  # Stage 2: burden tests need only syn + trunc
  t1 <- u1 %>% filter(!is.na(trunc_obs), !is.na(syn_obs))
  
  disc_n <- sum(t1$is_disc, na.rm = TRUE)
  rest_n <- sum(!t1$is_disc, na.rm = TRUE)
  if (disc_n < 5 || rest_n < 10) return(NULL)
  
  # Residualised trunc burden within current universe
  t1 <- add_trunc_resid(t1)
  
  disc <- t1 %>% filter(is_disc)
  rest <- t1 %>% filter(!is_disc)
  
  # Wilcoxon p-values
  p_raw   <- suppressWarnings(wilcox.test(disc$trunc_obs, rest$trunc_obs)$p.value)
  p_resid <- suppressWarnings(wilcox.test(disc$trunc_resid, rest$trunc_resid)$p.value)
  
  tibble(
    loeuf_cut = loeuf_cut,
    drop_ci_frac = drop_ci_frac,
    top_oe_frac = top_oe_frac,
    n_loeuf = nrow(u0),
    ci_thr = ci_thr,
    n_universe = nrow(u1),
    oe_thr = oe_thr,
    n_disc = sum(u1$is_disc, na.rm = TRUE),
    
    med_trunc_disc = median(disc$trunc_obs, na.rm = TRUE),
    med_trunc_rest = median(rest$trunc_obs, na.rm = TRUE),
    med_trunc_diff = med_trunc_disc - med_trunc_rest,
    
    med_resid_disc = median(disc$trunc_resid, na.rm = TRUE),
    med_resid_rest = median(rest$trunc_resid, na.rm = TRUE),
    med_resid_diff = med_resid_disc - med_resid_rest,
    
    p_trunc_raw = p_raw,
    p_trunc_resid = p_resid
  )
}

# ---------- run grid ----------
run_threshold_grid <- function(gene_df,
                               loeuf_cuts = c(0.10, 0.15, 0.20, 0.25),
                               drop_ci_fracs = c(0.10, 0.25, 0.40),
                               top_oe_fracs = c(0.025, 0.05, 0.10)) {
  
  grid <- tidyr::crossing(
    loeuf_cut = loeuf_cuts,
    drop_ci_frac = drop_ci_fracs,
    top_oe_frac = top_oe_fracs
  )
  
  # Track sizes/status for transparency
  counts <- purrr::pmap_dfr(grid, \(loeuf_cut, drop_ci_frac, top_oe_frac) {
    
    u0 <- gene_df %>%
      filter(!is.na(loeuf), !is.na(oe_lof), !is.na(ci_width)) %>%
      filter(loeuf < loeuf_cut)
    
    n0 <- nrow(u0)
    if (n0 < 50) {
      return(tibble(loeuf_cut, drop_ci_frac, top_oe_frac, n0, n1 = NA_integer_, n_disc = NA_integer_, status = "too_few_after_loeuf"))
    }
    
    ci_thr <- as.numeric(quantile(u0$ci_width, probs = 1 - drop_ci_frac, na.rm = TRUE))
    u1 <- u0 %>% filter(ci_width <= ci_thr)
    n1 <- nrow(u1)
    
    if (n1 < 30) {
      return(tibble(loeuf_cut, drop_ci_frac, top_oe_frac, n0, n1, n_disc = NA_integer_, status = "too_few_after_ci"))
    }
    
    oe_thr <- as.numeric(quantile(u1$oe_lof, probs = 1 - top_oe_frac, na.rm = TRUE))
    n_disc <- sum(u1$oe_lof >= oe_thr, na.rm = TRUE)
    
    tibble(loeuf_cut, drop_ci_frac, top_oe_frac, n0, n1, n_disc, status = "ok")
  })
  
  # Main results
  res <- purrr::pmap_dfr(grid, \(loeuf_cut, drop_ci_frac, top_oe_frac) {
    classify_paramset(gene_df, loeuf_cut, drop_ci_frac, top_oe_frac)
  })
  
  list(results = res, counts = counts)
}

# ----------------------------
# 2) Threshold sensitivity (categorical)
# ----------------------------
out <- run_threshold_grid(gene_df)

sens   <- out$results
counts <- out$counts

OUT_SENS <- file.path(OUTDIR, "threshold_sensitivity_summary.csv")
OUT_COUNTS <- file.path(OUTDIR, "sensitivity_counts.csv")

readr::write_csv(sens, OUT_SENS, na = "")
readr::write_csv(counts, OUT_COUNTS, na = "")

message("[WROTE] ", normalizePath(OUT_SENS, mustWork = FALSE))
message("[WROTE] ", normalizePath(OUT_COUNTS, mustWork = FALSE))
message("[QC] ok parameter sets: ", sum(counts$status == "ok"), " / ", nrow(counts))
message("[QC] sens rows written: ", nrow(sens))

# ----------------------------
# 3) Continuous discordance score (Option A)
#    No thresholds: discordance = excess truncation burden
#    beyond expectation given synonymous variation (mutability proxy)
# ----------------------------

DISC_LOEUF_MAX <- 0.20
DISC_DROP_CI_FRAC <- 0.25   # optional: match your main stability filter

d_score <- gene_df %>%
  filter(
    !is.na(loeuf), !is.na(oe_lof), !is.na(ci_width),
    !is.na(trunc_obs), !is.na(syn_obs)
  ) %>%
  filter(loeuf < DISC_LOEUF_MAX) %>%
  filter(syn_obs >= 50)   # <<<<<< ADD HERE


message("[DISC SCORE A] Genes with LOEUF < ", DISC_LOEUF_MAX, ": ", nrow(d_score))

# Optional: apply the same CI-width precision filter used in main analysis
ci_thr <- quantile(d_score$ci_width, probs = 1 - DISC_DROP_CI_FRAC, na.rm = TRUE)

d_score <- d_score %>%
  filter(ci_width <= ci_thr)

message("[DISC SCORE A] After dropping worst ", DISC_DROP_CI_FRAC*100,
        "% CI widths: ", nrow(d_score))

# ------------------------------------------------------------
# Continuous discordance score = excess truncation burden
# (residual of trunc_obs ~ syn_obs)
# ------------------------------------------------------------

fit_trunc <- lm(log1p(trunc_obs) ~ log1p(syn_obs), data = d_score)

d_score <- d_score %>%
  mutate(
    disc_score_resid = resid(fit_trunc),
    disc_score_rstud = rstudent(fit_trunc)   # reviewer-friendly: standardised outlier score
  )

# ------------------------------------------------------------
# Association with constraint metrics (sanity check)
# ------------------------------------------------------------

cor_score_oe <- suppressWarnings(
  cor.test(d_score$disc_score_rstud, d_score$oe_lof, method = "spearman")
)

cor_score_loeuf <- suppressWarnings(
  cor.test(d_score$disc_score_rstud, d_score$loeuf, method = "spearman")
)

score_summary <- tibble(
  regime = paste0("LOEUF < ", DISC_LOEUF_MAX,
                  ", drop worst ", DISC_DROP_CI_FRAC*100, "% CI widths"),
  n = nrow(d_score),
  
  spearman_rho_score_vs_oe_lof = unname(cor_score_oe$estimate),
  p_score_vs_oe_lof = cor_score_oe$p.value,
  
  spearman_rho_score_vs_loeuf = unname(cor_score_loeuf$estimate),
  p_score_vs_loeuf = cor_score_loeuf$p.value
)

OUT_SCORE <- file.path(OUTDIR, "discordance_score_summary.csv")
readr::write_csv(score_summary, OUT_SCORE, na = "")
message("[WROTE] ", normalizePath(OUT_SCORE, mustWork = FALSE))

# ------------------------------------------------------------
# Export top genes by continuous discordance score
# ------------------------------------------------------------

top_by_score <- d_score %>%
  arrange(desc(disc_score_rstud)) %>%
  select(
    gene, loeuf, oe_lof, ci_width,
    syn_obs, trunc_obs,
    disc_score_resid, disc_score_rstud
  ) %>%
  slice_head(n = 50)

OUT_TOP <- file.path(OUTDIR, "discordance_score_top50_genes.csv")
readr::write_csv(top_by_score, OUT_TOP, na = "")
message("[WROTE] ", normalizePath(OUT_TOP, mustWork = FALSE))

message("[DONE] Continuous discordance score (Option A: excess truncation) complete.")

######################## FIG

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

OUTDIR <- "tables"
TABDIR <- "tables"


SENS_FILE <- file.path("tables","threshold_sensitivity_summary.csv")

GF_RDS <- file.path("gnomad_lof_discordance_out","cache","gene_features.rds")
DISC_FILE <- file.path(TABDIR, "Supplementary_Table1_discordant_genes.csv")

stopifnot(file.exists(SENS_FILE), file.exists(GF_RDS), file.exists(DISC_FILE))

# ----------------------------
# Load inputs
# ----------------------------
sens <- readr::read_csv(SENS_FILE, show_col_types = FALSE) %>%
  mutate(
    sig = p_trunc_resid < 0.05,
    run_id = row_number(),
    note = ifelse(!sig,
                  paste0("n=", n_disc),
                  NA_character_)
  )

disc <- readr::read_csv(DISC_FILE, show_col_types = FALSE)$Gene
gf   <- readRDS(GF_RDS)

# ----------------------------
# Panel A: threshold robustness (zoomed)
# ----------------------------
sens <- sens %>%
  arrange(med_resid_diff) %>%
  mutate(run_id = row_number())

xmin <- min(sens$med_resid_diff) - 0.02
xmax <- max(sens$med_resid_diff) + 0.02

pA <- ggplot(sens, aes(x = med_resid_diff, y = run_id)) +
  
  # points
  geom_point(size = 2.8) +
  
  # reference at zero (off-scale)
  annotate("segment",
           x = xmin, xend = xmin,
           y = 0, yend = max(sens$run_id),
           linetype = "dashed") +
  
  annotate("text",
           x = xmin,
           y = max(sens$run_id) + 1,
           label = "0 (null)",
           size = 3,
           hjust = -0.1) +
  
  scale_y_continuous(
    breaks = seq(1, max(sens$run_id), by = 5)
  ) +
  
  coord_cartesian(xlim = c(xmin, xmax)) +
  
  labs(
    x = "Median residual truncation difference (discordant − rest)",
    y = "Threshold parameter set"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

# ----------------------------
# Panel B: continuous discordance score
# ----------------------------
DISC_LOEUF_MAX    <- 0.20
DISC_DROP_CI_FRAC <- 0.25
MIN_SYN           <- 50

d_score <- gf %>%
  filter(!is.na(gene),
         !is.na(loeuf),
         !is.na(ci_width),
         !is.na(trunc_obs),
         !is.na(syn_obs)) %>%
  filter(loeuf < DISC_LOEUF_MAX,
         syn_obs >= MIN_SYN)

ci_thr <- quantile(d_score$ci_width,
                   probs = 1 - DISC_DROP_CI_FRAC,
                   na.rm = TRUE)

d_score <- d_score %>% filter(ci_width <= ci_thr)

# Continuous score = studentised residual of trunc_obs ~ syn_obs
fit_trunc <- lm(log1p(trunc_obs) ~ log1p(syn_obs), data = d_score)

d_score <- d_score %>%
  mutate(
    disc_score = rstudent(fit_trunc),
    group = ifelse(toupper(gene) %in% toupper(disc),
                   "Categorical discordant",
                   "Other LOEUF<0.2")
  )

pB <- ggplot(d_score, aes(x = group, y = disc_score)) +
  geom_violin(trim = TRUE, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.size = 0.7) +
  labs(
    x = NULL,
    y = "Studentised excess truncation score"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

# ----------------------------
# Save combined figure (A/B only)
# ----------------------------
OUT_PDF <- file.path(OUTDIR, "SuppFig_robustness_clean.pdf")
OUT_PNG <- file.path(OUTDIR, "SuppFig_robustness_clean.png")

if (requireNamespace("patchwork", quietly = TRUE)) {
  suppressPackageStartupMessages(library(patchwork))
  
  combo <- pA / pB +
    plot_annotation(tag_levels = "A")
  
  ggsave(OUT_PDF, combo, width = 9, height = 7.5, dpi = 300)
  ggsave(OUT_PNG, combo, width = 9, height = 7.5, dpi = 300)
} else {
  ggsave(file.path(OUTDIR, "SuppFigA_robustness.pdf"),
         pA, width = 9, height = 3.5, dpi = 300)
  ggsave(file.path(OUTDIR, "SuppFigB_continuous_score.pdf"),
         pB, width = 9, height = 3.5, dpi = 300)
}

message("[DONE] Wrote cleaned supplementary robustness figure:")
message("  - ", OUT_PDF)
message("  - ", OUT_PNG)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(grid)   # for unit()
})

OUTDIR <- "tables"
TABDIR <- "tables"


GF_RDS <- file.path("gnomad_lof_discordance_out","cache","gene_features.rds")
DISC_FILE <- file.path(TABDIR, "Supplementary_Table1_discordant_genes.csv")
stopifnot(file.exists(GF_RDS), file.exists(DISC_FILE))

gf <- readRDS(GF_RDS)

disc_genes <- readr::read_csv(DISC_FILE, show_col_types = FALSE) %>%
  pull(Gene) %>%
  toupper()

# ----------------------------
# 1) Compute genome-wide truncation-excess score
#    (studentised residual of log1p(trunc_obs) ~ log1p(syn_obs))
# ----------------------------
d <- gf %>%
  filter(!is.na(gene), !is.na(loeuf), !is.na(trunc_obs), !is.na(syn_obs)) %>%
  mutate(gene_u = toupper(gene))

fit <- lm(log1p(trunc_obs) ~ log1p(syn_obs), data = d)

d <- d %>%
  mutate(
    trunc_excess_rstud = rstudent(fit),
    
    # Flip so higher = more excess truncation (intuitive axis)
    disc_score = -trunc_excess_rstud,
    
    is_discordant = gene_u %in% disc_genes,
    group = case_when(
      is_discordant ~ "Discordant (subset of LOEUF<0.2)",
      loeuf < 0.2   ~ "Other LOEUF<0.2",
      TRUE          ~ "Genome-wide (all genes)"
    ),
    group = factor(group, levels = c(
      "Genome-wide (all genes)",
      "Other LOEUF<0.2",
      "Discordant (subset of LOEUF<0.2)"
    ))
  )

# ----------------------------
# 2) Median line (ONLY discordant)
# ----------------------------
disc_median <- d %>%
  filter(group == "Discordant (subset of LOEUF<0.2)") %>%
  summarise(med = median(disc_score, na.rm = TRUE))

# ----------------------------
# 3) Colours (all reds)
# ----------------------------
red_pal <- c(
  "Genome-wide (all genes)"          = "#F4A6A6",  # light
  "Other LOEUF<0.2"                  = "#D55E5E",  # medium
  "Discordant (subset of LOEUF<0.2)" = "#7A0000"   # dark
)

# ----------------------------
# 4) ECDF plot (clean, no grid; discordant median only)
# ----------------------------
p <- ggplot(d, aes(x = disc_score, colour = group)) +
  stat_ecdf(linewidth = 1.2) +
  
  # Rug marks for discordant genes (subtle)
  geom_rug(
    data = d %>% filter(is_discordant),
    aes(x = disc_score),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.35,
    length = unit(0.02, "npc")
  ) +
  
  # ONLY the discordant median reference line
  geom_vline(
    data = disc_median,
    aes(xintercept = med),
    colour = red_pal[["Discordant (subset of LOEUF<0.2)"]],
    linewidth = 0.8,
    linetype = "dashed",
    alpha = 0.85
  ) +
  
  scale_colour_manual(values = red_pal) +
  
  labs(
    x = "Studentised excess truncation score",
    y = "Empirical cumulative distribution (ECDF)",
    colour = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.80, 0.25),
    legend.background = element_blank(),
    axis.line = element_line(linewidth = 0.6)
  )

OUT_PDF <- file.path(OUTDIR, "SuppFig_discScore_ECDF_genome_vs_loeuf_vs_discordant.pdf")
OUT_PNG <- file.path(OUTDIR, "SuppFig_discScore_ECDF_genome_vs_loeuf_vs_discordant.png")

ggsave(OUT_PDF, p, width = 10, height = 6, dpi = 300)
ggsave(OUT_PNG, p, width = 10, height = 6, dpi = 300)

message("[WROTE] ", normalizePath(OUT_PDF, mustWork = FALSE))
message("[WROTE] ", normalizePath(OUT_PNG, mustWork = FALSE))

if (FALSE) {
  # Assume you have per-gene tables:
  # gene_metrics: gene, LOEUF, lof_oe, ci_width, syn_obs, trunc_obs_total
  # and if available: trunc_obs_exome, trunc_obs_genome, syn_obs_exome, syn_obs_genome
  
  library(dplyr)
  
  # Residualised truncation score per subset
  score_subset <- function(df, trunc_col, syn_col) {
    m <- lm(log1p(df[[trunc_col]]) ~ log1p(df[[syn_col]]), data = df)
    df %>% mutate(disc_score = rstudent(m))
  }
  
  # Exomes
  ex <- gene_metrics %>%
    filter(!is.na(trunc_obs_exome), !is.na(syn_obs_exome)) %>%
    score_subset("trunc_obs_exome", "syn_obs_exome")
  
  # Genomes
  gn <- gene_metrics %>%
    filter(!is.na(trunc_obs_genome), !is.na(syn_obs_genome)) %>%
    score_subset("trunc_obs_genome", "syn_obs_genome")
  
  # Compare effect direction in each subset
  test_effect <- function(df, disc_genes) {
    df %>% mutate(is_disc = gene %in% disc_genes) %>%
      summarise(
        p = wilcox.test(disc_score ~ is_disc)$p.value,
        med_disc = median(disc_score[is_disc]),
        med_bg   = median(disc_score[!is_disc])
      )
  }
  
  # Spearman correlation across subsets (matched genes)
  joined <- inner_join(
    ex %>% select(gene, disc_score_ex = disc_score),
    gn %>% select(gene, disc_score_gn = disc_score),
    by = "gene"
  )
  cor.test(joined$disc_score_ex, joined$disc_score_gn, method = "spearman")
}

##############
if (FALSE) {
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# ----------------------------
# USER INPUTS
# ----------------------------
CLINGEN_BED <- "ClinGen_haploinsufficiency_gene_GRCh38.bed"
# uc2 must exist and contain columns: gene, is_disc (logical)
# uc2 should represent your LOEUF<0.2 + precision-filtered universe (e.g., n ~ 410)

# ----------------------------
# 0) Sanity checks
# ----------------------------
if (!exists("uc2")) stop("Object 'uc2' not found. Load/create uc2 before running.")
req_cols <- c("gene", "is_disc")
missing <- setdiff(req_cols, names(uc2))
if (length(missing) > 0) stop("uc2 is missing required columns: ", paste(missing, collapse = ", "))

if (!file.exists(CLINGEN_BED)) {
  stop("ClinGen BED not found at: ", CLINGEN_BED, "\nSet CLINGEN_BED to the correct path.")
}

# ----------------------------
# 1) Load ClinGen BED
# Format: chr, start, end, gene_symbol, HI_score
# ----------------------------
clingen_raw <- read_tsv(
  CLINGEN_BED,
  comment = "track",
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(clingen_raw) < 5) stop("Expected >=5 columns in ClinGen BED. Found: ", ncol(clingen_raw))

clingen_raw <- clingen_raw %>%
  transmute(
    chr = X1,
    start = X2,
    end = X3,
    gene = toupper(as.character(X4)),
    HI = as.numeric(X5)
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(HI))

cat("ClinGen rows loaded:", nrow(clingen_raw), "\n")
cat("ClinGen HI score distribution:\n")
print(table(clingen_raw$HI, useNA = "ifany"))

# ----------------------------
# 2) Define ClinGen gene sets
# HI=3 (sufficient evidence), HI>=2 (some or sufficient evidence)
# ----------------------------
hi3_genes <- clingen_raw %>%
  filter(HI == 3) %>%
  pull(gene) %>%
  unique()

hi2plus_genes <- clingen_raw %>%
  filter(HI %in% c(2, 3)) %>%
  pull(gene) %>%
  unique()

cat("\nClinGen genes:\n")
cat("  HI=3 genes:", length(hi3_genes), "\n")
cat("  HI>=2 genes:", length(hi2plus_genes), "\n")

# ----------------------------
# 3) Join ClinGen scores to uc2 and compute enrichment flags
# ----------------------------
uc2_work <- uc2 %>%
  mutate(geneU = toupper(as.character(gene))) %>%
  left_join(
    clingen_raw %>% distinct(gene, HI),
    by = c("geneU" = "gene")
  ) %>%
  mutate(
    in_clingen_hi3 = !is.na(HI) & HI == 3,
    in_clingen_hi2plus = !is.na(HI) & HI %in% c(2, 3)
  )

cat("\nuc2 size:", nrow(uc2_work), "\n")
cat("Discordant genes in uc2:", sum(uc2_work$is_disc, na.rm = TRUE), "\n")
cat("Genes with any ClinGen HI score (non-NA):", sum(!is.na(uc2_work$HI)), "\n")

# ----------------------------
# 4) Fisher tests (HI>=2 and HI=3)
# ----------------------------
tab_hi2plus <- table(uc2_work$is_disc, uc2_work$in_clingen_hi2plus)
cat("\nContingency table: is_disc x ClinGen HI>=2\n")
print(tab_hi2plus)

if (all(dim(tab_hi2plus) == c(2, 2))) {
  cat("\nFisher test (HI>=2):\n")
  print(fisher.test(tab_hi2plus))
} else {
  cat("\nFisher test (HI>=2) not run: table is not 2x2.\n")
}

tab_hi3 <- table(uc2_work$is_disc, uc2_work$in_clingen_hi3)
cat("\nContingency table: is_disc x ClinGen HI=3\n")
print(tab_hi3)

if (all(dim(tab_hi3) == c(2, 2))) {
  cat("\nFisher test (HI=3):\n")
  print(fisher.test(tab_hi3))
} else {
  cat("\nFisher test (HI=3) not run: table is not 2x2.\n")
}

# ----------------------------
# 5) Optional: Wilcoxon test on HI score distribution
# Note: many genes may have NA (not curated). This test uses non-NA only.
# ----------------------------
uc2_non_na <- uc2_work %>% filter(!is.na(HI))

cat("\nNon-NA HI scores in uc2:", nrow(uc2_non_na), "\n")
if (length(unique(uc2_non_na$is_disc)) == 2) {
  cat("\nWilcoxon test of HI score by discordant status (non-NA only):\n")
  print(wilcox.test(HI ~ is_disc, data = uc2_non_na))
} else {
  cat("\nWilcoxon test not run: only one discordant group present among non-NA HI genes.\n")
}

# ----------------------------
# 6) Helpful summaries for writing
# ----------------------------
summ_hi <- uc2_work %>%
  group_by(is_disc) %>%
  summarise(
    n_genes = n(),
    n_hi_scored = sum(!is.na(HI)),
    prop_hi_scored = n_hi_scored / n_genes,
    n_hi3 = sum(in_clingen_hi3, na.rm = TRUE),
    prop_hi3 = n_hi3 / n_genes,
    n_hi2plus = sum(in_clingen_hi2plus, na.rm = TRUE),
    prop_hi2plus = n_hi2plus / n_genes,
    .groups = "drop"
  )

cat("\nSummary by discordant status:\n")
print(summ_hi)
}

#########
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

INFILE <- "trunc_position_domain_nmd_variants.csv"
OUTDIR <- "tables"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

d <- readr::read_csv(INFILE, show_col_types = FALSE)

DISC_FILE <- file.path(TABDIR, "Supplementary_Table1_discordant_genes.csv")
stopifnot(file.exists(DISC_FILE))

disc35 <- readr::read_csv(DISC_FILE, show_col_types = FALSE) %>%
  dplyr::pull(Gene) %>%
  toupper() %>%
  unique()

message("[positional] Loaded discordant genes: ", length(disc35))

# ---- sanity ----
stopifnot(all(c("Gene","rel_pos") %in% names(d)))

# ============================================================
# Build 2-group table
# Prefer disc35 if it exists; otherwise use d$group if it is 2-level
# ============================================================

dd <- d %>%
  dplyr::transmute(
    geneU = toupper(as.character(Gene)),
    pos   = suppressWarnings(as.numeric(rel_pos))
  ) %>%
  dplyr::filter(is.finite(pos), !is.na(geneU), geneU != "") %>%
  dplyr::mutate(
    group = dplyr::if_else(geneU %in% disc35, "Discordant", "Background"),
    group = factor(group, levels = c("Background", "Discordant"))
  )

# ---- diagnostics ----
cat("\nGroup counts (rows/variants):\n")
print(table(dd$group))

cat("\nUnique genes by group:\n")
print(dd %>% distinct(group, geneU) %>% count(group))

# Ensure both groups present
if (sum(dd$group == "Discordant") < 2 || sum(dd$group == "Background") < 2) {
  message("[positional] Variant rows per group:")
  print(table(dd$group, useNA = "ifany"))
  
  message("[positional] Unique genes per group:")
  print(dd %>% dplyr::distinct(group, geneU) %>% dplyr::count(group))
  
  message("[positional] Example discordant genes not matching file (first 20):")
  print(head(setdiff(disc35, unique(dd$geneU)), 20))
  
  message("[positional] Example file genes (first 20):")
  print(head(unique(dd$geneU), 20))
  stop("One of the groups has <2 data points. Check that disc35 matches Gene symbols, or that 'group' parsing worked.\n",
       "Tip: try using Gene_norm instead of Gene, or print head(disc35) vs head(unique(d$Gene)).")
}

x_bg   <- dd %>% filter(group == "Background") %>% pull(pos)
x_disc <- dd %>% filter(group == "Discordant") %>% pull(pos)

# ============================================================
# (2) KS test
# ============================================================
ks <- stats::ks.test(x_disc, x_bg)

ks_out <- tibble(
  test = "KS",
  n_discordant = length(x_disc),
  n_background = length(x_bg),
  statistic_D = unname(ks$statistic),
  p_value = ks$p.value
)
write_csv(ks_out, file.path(OUTDIR, "KS_test_trunc_position.csv"))
print(ks_out)

# ============================================================
# (3) Logistic regression for pre-specified late third
# ============================================================
dd2 <- dd %>%
  mutate(
    late_third = pos >= (2/3),
    is_disc = group == "Discordant"
  )

fit <- glm(late_third ~ is_disc, data = dd2, family = binomial())

beta <- coef(summary(fit))["is_discTRUE", "Estimate"]
se   <- coef(summary(fit))["is_discTRUE", "Std. Error"]
or   <- exp(beta)
ci   <- exp(beta + c(-1, 1) * 1.96 * se)
pval <- coef(summary(fit))["is_discTRUE", "Pr(>|z|)"]

logit_out <- tibble(
  model = "late_third ~ discordant",
  cutoff = "pos >= 2/3",
  n_discordant = sum(dd2$is_disc),
  n_background = sum(!dd2$is_disc),
  OR = or,
  CI_low = ci[1],
  CI_high = ci[2],
  p_value = pval
)
write_csv(logit_out, file.path(OUTDIR, "logistic_late_third_OR.csv"))
print(logit_out)

# ============================================================
# (4) ECDF + ΔCDF effect-size plots
# ============================================================

p_ecdf <- ggplot(dd, aes(x = pos, linetype = group)) +
  stat_ecdf(linewidth = 0.8) +
  labs(
    x = "Relative truncation position (rel_pos)",
    y = "Empirical CDF",
    linetype = NULL
  ) +
  theme_classic()

ggsave(file.path(OUTDIR, "ECDF_trunc_position.png"), p_ecdf, width = 6.5, height = 4.5, dpi = 300)

grid <- seq(0, 1, by = 0.001)
ecdf_bg   <- stats::ecdf(x_bg)
ecdf_disc <- stats::ecdf(x_disc)

delta <- tibble(
  x = grid,
  cdf_bg = ecdf_bg(grid),
  cdf_disc = ecdf_disc(grid),
  delta_cdf = cdf_disc - cdf_bg
)

max_sep <- delta %>% slice_max(abs(delta_cdf), n = 1)

write_csv(delta,   file.path(OUTDIR, "delta_CDF_table.csv"))
write_csv(max_sep, file.path(OUTDIR, "delta_CDF_max_separation.csv"))

print(max_sep)

p_delta <- ggplot(delta, aes(x = x, y = delta_cdf)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "Relative truncation position (rel_pos)",
    y = "ΔCDF (Discordant − Background)"
  ) +
  theme_classic()

# ============================================================
# (4b) ECDF + ΔCDF as native ggplots + tight patchwork combine
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Make sure dd and delta exist (from your earlier code)
stopifnot(exists("dd"), exists("delta"))
stopifnot(all(c("pos","group") %in% names(dd)))
stopifnot(all(c("x","delta_cdf") %in% names(delta)))



# ----------------------------
# Panel A: ECDF
# ----------------------------
p_ecdf <- ggplot(dd, aes(x = pos, linetype = group)) +
  stat_ecdf(linewidth = 0.8) +
  labs(
    x = "Relative truncation position (rel_pos)",
    y = "Empirical CDF",
    linetype = NULL
  ) +
  theme_classic() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 12),
    legend.position = "right",
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

# ----------------------------
# Panel B: ΔCDF (Discordant - Background)
# ----------------------------
p_delta <- ggplot(delta, aes(x = x, y = delta_cdf)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  labs(
    x = "Relative truncation position (rel_pos)",
    y = expression(Delta*CDF~"(Discordant \u2212 Background)")
  ) +
  theme_classic() +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0, size = 12),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5)
  )

p_ecdf  <- p_ecdf  + labs(title = NULL)
p_delta <- p_delta + labs(title = NULL)

# ----------------------------
# Combine with patchwork
#   - Titles are the filenames
#   - Tag letters A/B are automatic and tight
# ----------------------------
combined <- p_ecdf + p_delta +
  plot_layout(ncol = 2, widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 14, face = "bold"),
      plot.tag.position = c(0.01, 0.99),  # top-left inside each panel
      plot.margin = margin(2, 2, 2, 2)
    )
  )


# Save
out_png <- file.path(OUTDIR, "SuppFig_trunc_position_ECDF_and_DeltaCDF_NATIVE.png")
out_pdf <- file.path(OUTDIR, "SuppFig_trunc_position_ECDF_and_DeltaCDF_NATIVE.pdf")

ggsave(out_png, combined, width = 12, height = 5, dpi = 300)
ggsave(out_pdf, combined, width = 12, height = 5)

message("Native combined figure written to:\n  ", out_png, "\n  ", out_pdf)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

# ============================================================
# Inputs
# ============================================================

OUTDIR <- "tables"
TABLEDIR <- "tables"
dir.create(TABLEDIR, showWarnings = FALSE, recursive = TRUE)

# Files produced by your positional test script
f_ks    <- file.path(OUTDIR, "KS_test_trunc_position.csv")
f_logit <- file.path(OUTDIR, "logistic_late_third_OR.csv")

# ============================================================
# (1) Load test results
# ============================================================

ks <- read_csv(f_ks, show_col_types = FALSE)

logit <- read_csv(f_logit, show_col_types = FALSE)

# ---- Wilcoxon result ----
# You currently have Wilcoxon p=0.084 from earlier.
# Put it explicitly here so it lives in the reproducible table.

wilcox <- tibble(
  test = "Wilcoxon rank-sum",
  comparison = "Median rel_pos (Discordant vs Background)",
  effect = NA_character_,
  statistic = NA_character_,
  p_value = 0.084
)

# ============================================================
# (2) Format KS row
# ============================================================

ks_row <- ks %>%
  transmute(
    test = "Kolmogorov–Smirnov",
    comparison = "Distribution of rel_pos (Discordant vs Background)",
    effect = paste0("D = ", signif(statistic_D, 3)),
    statistic = paste0("n_disc = ", n_discordant,
                       ", n_bg = ", n_background),
    p_value = p_value
  )

# ============================================================
# (3) Format Logistic regression row
# ============================================================

logit_row <- logit %>%
  transmute(
    test = "Logistic regression",
    comparison = "Late truncation (rel_pos ≥ 2/3)",
    effect = paste0(
      "OR = ", signif(OR, 3),
      " (", signif(CI_low, 3), "–", signif(CI_high, 3), ")"
    ),
    statistic = paste0("n_disc = ", n_discordant,
                       ", n_bg = ", n_background),
    p_value = p_value
  )

# ============================================================
# (4) Combine into single supplementary table
# ============================================================

supp_table <- bind_rows(
  wilcox,
  ks_row,
  logit_row
) %>%
  mutate(
    p_value = format.pval(p_value, digits = 3, eps = 1e-10)
  )

# ============================================================
# (5) Write outputs
# ============================================================

outfile_csv <- file.path(TABLEDIR, "SuppTable_truncation_position_tests.csv")
outfile_tsv <- file.path(TABLEDIR, "SuppTable_truncation_position_tests.tsv")

write_csv(supp_table, outfile_csv)
write_tsv(supp_table, outfile_tsv)

message("✅ Supplementary positional test table written to:")
message(" - ", outfile_csv)
message(" - ", outfile_tsv)

# Print preview
print(supp_table)

