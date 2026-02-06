#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)   # <-- ADD THIS
  library(ggplot2)
})

# ============================================================
# 2-panel figure: CDS length vs truncation burden (LOEUF<0.2)
# Panel A: Scatter (CDS kb vs truncating count)
# Panel B: Violin/box (truncating per kb CDS)
#
# Inputs:
# - gnomad.v4.1.constraint_metrics.tsv
# - gnomad_lof_discordance_out/Table1_top_discordant_genes_QC_skeleton.csv
#   (or switch USE_TABLE1_TOP15 <- FALSE to use discordant_genes_top_oe_lof.csv)
#
# Output:
# - SuppFig_length_adjusted_2panel.png / .pdf
# - length_adjusted_stats.csv
# ============================================================

CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"

# ---- Choose discordant source ----
USE_TABLE1_TOP15 <- TRUE
TABLE1_N         <- 15

TABLE1_CSV     <- "gnomad_lof_discordance_out/Table1_top_discordant_genes_QC_skeleton.csv"
DISCORDANT_CSV <- "gnomad_lof_discordance_out/discordant_genes_top_oe_lof.csv"

# ---- Parameters ----
LOEUF_CUTOFF <- 0.20

OUT_STATS <- "length_adjusted_stats.csv"
OUT_PNG   <- "SuppFig_length_adjusted_2panel.png"
OUT_PDF   <- "SuppFig_length_adjusted_2panel.pdf"

# Plot styling
BASE_SIZE <- 14
POINT_GREY_ALPHA <- 0.35
POINT_GREY_SIZE  <- 1.0
POINT_RED_SIZE   <- 2.2

norm_sym <- function(x) toupper(trimws(as.character(x)))

# ---- Require patchwork (cleanest for 2-panel) ----
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
}

stopifnot(file.exists(CONSTRAINT_TSV))
if (USE_TABLE1_TOP15) stopifnot(file.exists(TABLE1_CSV)) else stopifnot(file.exists(DISCORDANT_CSV))

# ----------------------------
# Load discordant genes
# ----------------------------
if (USE_TABLE1_TOP15) {
  t1 <- read_csv(TABLE1_CSV, show_col_types = FALSE)
  stopifnot("Gene" %in% names(t1))
  disc_genes <- t1 %>%
    transmute(Gene = norm_sym(Gene)) %>%
    slice_head(n = TABLE1_N) %>%
    pull(Gene)
  discordant_label <- paste0("Table 1 top ", TABLE1_N)
} else {
  dsrc <- read_csv(DISCORDANT_CSV, show_col_types = FALSE)
  gene_col <- if ("gene" %in% names(dsrc)) "gene" else if ("Gene" %in% names(dsrc)) "Gene" else NULL
  if (is.null(gene_col)) stop("No gene/Gene column found in DISCORDANT_CSV.")
  disc_genes <- dsrc %>%
    transmute(Gene = norm_sym(.data[[gene_col]])) %>%
    distinct() %>%
    pull(Gene)
  discordant_label <- paste0("Discordant set (n=", length(disc_genes), ")")
}

message("[INFO] Discordant source: ", discordant_label)

# ----------------------------
# Load constraint data + pick columns
# ----------------------------
df_raw <- read_tsv(CONSTRAINT_TSV, show_col_types = FALSE)

stopifnot("gene" %in% names(df_raw))
stopifnot("lof.oe_ci.upper" %in% names(df_raw))

lof_obs_col <- if ("lof_obs" %in% names(df_raw)) {
  "lof_obs"
} else if ("lof_hc_lc.obs" %in% names(df_raw)) {
  "lof_hc_lc.obs"
} else if ("lof.obs" %in% names(df_raw)) {
  "lof.obs"
} else {
  stop("No LoF observed-count column found (expected lof_obs or lof_hc_lc.obs or lof.obs).")
}

cds_col <- if ("cds_length" %in% names(df_raw)) "cds_length" else NULL
if (is.null(cds_col)) stop("No cds_length column found in constraint file.")

message("[INFO] using lof_obs column: ", lof_obs_col)
message("[INFO] using cds_length column: ", cds_col)

# ----------------------------
# Collapse to 1 row per gene (min LOEUF per gene)
# ----------------------------
d <- df_raw %>%
  transmute(
    Gene    = norm_sym(gene),
    LOEUF   = suppressWarnings(as.numeric(`lof.oe_ci.upper`)),
    lof_obs = suppressWarnings(as.numeric(.data[[lof_obs_col]])),
    cds_len = suppressWarnings(as.numeric(.data[[cds_col]]))
  ) %>%
  filter(is.finite(LOEUF), is.finite(lof_obs), lof_obs >= 0, is.finite(cds_len), cds_len > 0) %>%
  arrange(Gene, LOEUF) %>%
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

stopifnot(nrow(d) == n_distinct(d$Gene))

# ----------------------------
# Restrict to LOEUF<0.2 universe + compute length-normalised burden
# ----------------------------
u <- d %>%
  filter(LOEUF < LOEUF_CUTOFF) %>%
  mutate(
    discordant = Gene %in% disc_genes,
    cds_kb = cds_len / 1000,
    trunc_per_kb = lof_obs / cds_kb
  )

message("[QC] LOEUF<", LOEUF_CUTOFF, " genes: n=", nrow(u))
message("[QC] discordant within LOEUF<", LOEUF_CUTOFF, ": n=", sum(u$discordant))

# Hard assertion for the Table1-top15 mode (since you define them as LOEUF<0.2)
if (USE_TABLE1_TOP15) {
  stopifnot(all(u$Gene[u$discordant] %in% disc_genes))
  # This additionally ensures the plotted representative row obeys LOEUF<0.2
  stopifnot(all(u$LOEUF[u$discordant] < LOEUF_CUTOFF))
}

# ----------------------------
# Stats: Wilcoxon for trunc_per_kb
# ----------------------------
w <- wilcox.test(trunc_per_kb ~ discordant, data = u, exact = FALSE)

stats <- u %>%
  group_by(discordant) %>%
  summarise(
    n = n(),
    median_trunc_per_kb = median(trunc_per_kb, na.rm = TRUE),
    iqr_trunc_per_kb = IQR(trunc_per_kb, na.rm = TRUE),
    median_cds_kb = median(cds_kb, na.rm = TRUE),
    median_lof_obs = median(lof_obs, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    wilcox_p = w$p.value,
    loeuf_cutoff = LOEUF_CUTOFF,
    discordant_source = ifelse(USE_TABLE1_TOP15, paste0("Table1_top", TABLE1_N), "discordant_genes_top_oe_lof.csv")
  )

write_csv(stats, OUT_STATS)
message("[DONE] wrote: ", OUT_STATS)

# ----------------------------
# Panel A: Scatter (CDS kb vs trunc count)
# ----------------------------
pA <- ggplot(u, aes(x = cds_kb, y = lof_obs)) +
  geom_point(alpha = POINT_GREY_ALPHA, size = POINT_GREY_SIZE, color = "grey70") +
  geom_point(
    data = u %>% filter(discordant),
    size = POINT_RED_SIZE, shape = 21, fill = "red3", color = "black", stroke = 0.35
  ) +
  scale_y_log10() +
  labs(
    title = "Discordant genes are not explained by CDS length alone",
    subtitle = paste0("Universe: LOEUF < ", LOEUF_CUTOFF),
    x = "CDS length (kb)",
    y = "Observed truncating variants per gene (log10)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------
# Panel B: Violin/box (trunc per kb)
# ----------------------------
u2 <- u %>%
  mutate(group = ifelse(discordant, "Discordant", paste0("Other LOEUF<", LOEUF_CUTOFF)))

pB <- ggplot(u2, aes(x = group, y = trunc_per_kb)) +
  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.alpha = 0.25) +
  scale_y_log10() +
  labs(
    title = "Truncating burden remains elevated after CDS-length normalisation",
    subtitle = paste0("Universe: LOEUF < ", LOEUF_CUTOFF, " | Wilcoxon p = ", format(w$p.value, digits = 3)),
    x = NULL,
    y = "Observed truncating variants per kb CDS (log10)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------
# Combine into 2-panel figure
# ----------------------------
combined <- pA + pB +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1)) +
  patchwork::plot_annotation(tag_levels = "A")  # adds A, B tags

ggsave(OUT_PNG, combined, width = 12.8, height = 5.6, units = "in", dpi = 300)
ggsave(OUT_PDF, combined, width = 12.8, height = 5.6, units = "in")

message("[DONE] wrote: ", OUT_PNG)
message("[DONE] wrote: ", OUT_PDF)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Supplementary 2-panel figure:
# (A) Histogram of null medians + observed line
# (B) ECDF of null medians + observed line
#
# Question: Does the discordant gene set have a median truncation burden
# greater than expected under random sampling from LOEUF<0.2 genes?
#
# Outputs:
# - SuppFig_perm_null_2panel.png/.pdf
# - perm_null_trunc_burden_summary.csv
# - perm_null_trunc_burden_null_values.csv
# ============================================================

# ----------------------------
# Inputs
# ----------------------------
CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"

# Discordant set source:
USE_DISCORDANT_RDS <- TRUE
DISCORDANT_RDS <- "discordant_genes.rds"

# Optional alternative:
TABLE1_CSV <- "gnomad_lof_discordance_out/Table1_top_discordant_genes_QC_skeleton.csv"
TABLE1_N   <- 15

# ----------------------------
# Parameters
# ----------------------------
LOEUF_CUTOFF <- 0.20
N_REPS       <- 1000
SEED         <- 1

# Statistic:
#   "raw"    = median(lof_obs) across genes
#   "per_kb" = median(lof_obs / (cds_length/1000)) across genes
STAT_MODE <- "per_kb"  # set to "raw" if you prefer

# Outputs
OUT_PNG <- "SuppFig_perm_null_2panel.png"
OUT_PDF <- "SuppFig_perm_null_2panel.pdf"

OUT_SUMMARY <- "perm_null_trunc_burden_summary.csv"
OUT_NULLCSV <- "perm_null_trunc_burden_null_values.csv"

# Styling
BASE_SIZE <- 14

norm_sym <- function(x) toupper(trimws(as.character(x)))

# ---- Require patchwork (for 2-panel layout) ----
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
}

stopifnot(file.exists(CONSTRAINT_TSV))
if (USE_DISCORDANT_RDS) stopifnot(file.exists(DISCORDANT_RDS)) else stopifnot(file.exists(TABLE1_CSV))

# ----------------------------
# Load discordant genes
# ----------------------------
if (USE_DISCORDANT_RDS) {
  disc <- readRDS(DISCORDANT_RDS)
  stopifnot("Gene" %in% names(disc))
  disc_genes <- disc %>% transmute(Gene = norm_sym(Gene)) %>% distinct() %>% pull(Gene)
  disc_label <- paste0("Discordant set (RDS; n=", length(disc_genes), ")")
} else {
  t1 <- read_csv(TABLE1_CSV, show_col_types = FALSE)
  stopifnot("Gene" %in% names(t1))
  disc_genes <- t1 %>% transmute(Gene = norm_sym(Gene)) %>% slice_head(n = TABLE1_N) %>% pull(Gene)
  disc_label <- paste0("Table1 top ", TABLE1_N)
}
N_SET <- length(disc_genes)

message("[INFO] discordant source: ", disc_label)
message("[INFO] N in observed set: ", N_SET)

# ----------------------------
# Load constraint table
# ----------------------------
df_raw <- read_tsv(CONSTRAINT_TSV, show_col_types = FALSE)

stopifnot("gene" %in% names(df_raw))
stopifnot("lof.oe_ci.upper" %in% names(df_raw))

lof_obs_col <- if ("lof_obs" %in% names(df_raw)) {
  "lof_obs"
} else if ("lof_hc_lc.obs" %in% names(df_raw)) {
  "lof_hc_lc.obs"
} else if ("lof.obs" %in% names(df_raw)) {
  "lof.obs"
} else {
  stop("No LoF observed-count column found (expected lof_obs or lof_hc_lc.obs or lof.obs).")
}

cds_col <- if ("cds_length" %in% names(df_raw)) "cds_length" else NULL
if (STAT_MODE == "per_kb" && is.null(cds_col)) {
  stop("STAT_MODE='per_kb' requires cds_length column, but it was not found.")
}

message("[INFO] using LoF observed column: ", lof_obs_col)
if (!is.null(cds_col)) message("[INFO] using CDS length column: ", cds_col)

# ----------------------------
# Collapse to 1 row per gene (min LOEUF per gene)
# ----------------------------
d <- df_raw %>%
  transmute(
    Gene    = norm_sym(gene),
    LOEUF   = suppressWarnings(as.numeric(`lof.oe_ci.upper`)),
    lof_obs = suppressWarnings(as.numeric(.data[[lof_obs_col]])),
    cds_len = if (!is.null(cds_col)) suppressWarnings(as.numeric(.data[[cds_col]])) else NA_real_
  ) %>%
  filter(is.finite(LOEUF), is.finite(lof_obs), lof_obs >= 0) %>%
  { if (STAT_MODE == "per_kb") filter(., is.finite(cds_len), cds_len > 0) else . } %>%
  arrange(Gene, LOEUF) %>%
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

stopifnot(nrow(d) == n_distinct(d$Gene))

# Universe LOEUF<0.2
u <- d %>% filter(LOEUF < LOEUF_CUTOFF)
message("[QC] universe LOEUF<", LOEUF_CUTOFF, ": n_genes=", nrow(u))

# Ensure discordant genes are in universe; intersect if needed
missing <- setdiff(disc_genes, u$Gene)
if (length(missing) > 0) {
  message("[WARN] discordant genes missing from LOEUF< cutoff universe: n=", length(missing))
  message("[WARN] examples: ", paste(head(missing, 20), collapse = ", "))
  disc_genes <- intersect(disc_genes, u$Gene)
  N_SET <- length(disc_genes)
  message("[WARN] using intersected discordant set size: n=", N_SET)
}
if (N_SET < 5) stop("Too few discordant genes in universe after intersection; check inputs.")

# Burden per gene
u <- u %>%
  mutate(
    burden = if (STAT_MODE == "raw") {
      lof_obs
    } else {
      lof_obs / (cds_len / 1000)
    }
  )

obs_median <- median(u$burden[u$Gene %in% disc_genes], na.rm = TRUE)
message("[QC] observed median burden (", STAT_MODE, "): ", signif(obs_median, 5))

# ----------------------------
# Null distribution via random sampling
# ----------------------------
set.seed(SEED)
gene_pool <- u$Gene

null_medians <- numeric(N_REPS)
for (i in seq_len(N_REPS)) {
  samp <- sample(gene_pool, size = N_SET, replace = FALSE)
  null_medians[i] <- median(u$burden[u$Gene %in% samp], na.rm = TRUE)
}

# Empirical one-sided p: P(null >= observed)
p_emp <- (sum(null_medians >= obs_median) + 1) / (N_REPS + 1)

# Save null values
null_tbl <- tibble(rep = seq_len(N_REPS), null_median = null_medians)
write_csv(null_tbl, OUT_NULLCSV)

# Summary
summary_tbl <- tibble(
  loeuf_cutoff = LOEUF_CUTOFF,
  stat_mode = STAT_MODE,
  discordant_source = disc_label,
  n_discordant_used = N_SET,
  n_universe = nrow(u),
  reps = N_REPS,
  seed = SEED,
  observed_median = obs_median,
  null_median_mean = mean(null_medians),
  null_median_sd = sd(null_medians),
  null_median_q025 = as.numeric(quantile(null_medians, 0.025)),
  null_median_q975 = as.numeric(quantile(null_medians, 0.975)),
  p_emp_one_sided = p_emp
)
write_csv(summary_tbl, OUT_SUMMARY)

message("[DONE] wrote: ", OUT_NULLCSV)
message("[DONE] wrote: ", OUT_SUMMARY)
print(summary_tbl)

# ----------------------------
# Labels
# ----------------------------
xlab <- if (STAT_MODE == "raw") {
  "Median observed truncating variants per gene"
} else {
  "Median truncating variants per kb CDS"
}

main_title <- "Discordant genes exceed null expectation under random sampling"
sub_title  <- paste0(
  "Universe: LOEUF < ", LOEUF_CUTOFF,
  " | N=", N_SET,
  " | reps=", N_REPS,
  " | empirical one-sided p=", format(p_emp, digits = 3)
)

# ----------------------------
# Panel A: Histogram + observed line
# ----------------------------
pA <- ggplot(null_tbl, aes(x = null_median)) +
  geom_histogram(bins = 40, alpha = 0.85) +
  geom_vline(xintercept = obs_median, linewidth = 0.9) +
  labs(
    title = "A. Null distribution of median truncation burden",
    x = xlab,
    y = "Number of resamples"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------
# Panel B: ECDF + observed line
# ----------------------------
pB <- ggplot(null_tbl, aes(x = null_median)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = obs_median, linewidth = 0.9) +
  labs(
    title = "B. Empirical CDF of null medians",
    x = xlab,
    y = "Cumulative probability"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------
# Combine to 2-panel figure
# ----------------------------
combined <- (pA + pB) +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1)) +
  patchwork::plot_annotation(title = main_title, subtitle = sub_title)

ggsave(OUT_PNG, combined, width = 13.2, height = 5.4, units = "in", dpi = 300)
ggsave(OUT_PDF, combined, width = 13.2, height = 5.4, units = "in")

message("[DONE] wrote: ", OUT_PNG)
message("[DONE] wrote: ", OUT_PDF)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Negative control figure: synonymous / missense vs truncation
#
# Panel A: LOEUF vs synonymous observed burden (log scale)
#          Expect discordant genes do NOT stand out here.
#
# Panel B: truncating burden (lof_obs) vs missense constraint (mis.oe or MIS_OE_UPPER)
#          (fallback to missense observed if o/e isn't present)
#
# Inputs:
# - gnomad.v4.1.constraint_metrics.tsv
# - discordant_genes.rds (must contain column "Gene")
#
# Output:
# - SuppFig_negative_control_syn_mis.png / .pdf
# - negctrl_qc_summary.csv
# ============================================================

CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"
DISCORDANT_RDS <- "discordant_genes.rds"

OUT_PNG <- "SuppFig_negative_control_syn_mis.png"
OUT_PDF <- "SuppFig_negative_control_syn_mis.pdf"
OUT_QC  <- "negctrl_qc_summary.csv"

LOEUF_CUTOFF <- 0.20
SEED <- 1

BASE_SIZE <- 14

norm_sym <- function(x) toupper(trimws(as.character(x)))

# patchwork for 2-panel layout
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Package 'patchwork' is required. Install with: install.packages('patchwork')")
}

stopifnot(file.exists(CONSTRAINT_TSV))
stopifnot(file.exists(DISCORDANT_RDS))

# ----------------------------
# Load discordant gene symbols
# ----------------------------
disc <- readRDS(DISCORDANT_RDS)
stopifnot("Gene" %in% names(disc))
disc_genes <- disc %>%
  transmute(Gene = norm_sym(Gene)) %>%
  distinct() %>%
  pull(Gene)

message("[INFO] discordant genes loaded: n=", length(disc_genes))

# ----------------------------
# Load constraint table
# ----------------------------
df_raw <- read_tsv(CONSTRAINT_TSV, show_col_types = FALSE)

stopifnot("gene" %in% names(df_raw))
stopifnot("lof.oe_ci.upper" %in% names(df_raw))

# Helper: pick the first existing column from candidates
pick_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms][1]
  if (is.na(hit) || is.null(hit)) return(NULL)
  hit
}

nms <- names(df_raw)

# LoF observed (truncating burden proxy)
lof_obs_col <- pick_col(nms, c("lof_obs", "lof_hc_lc.obs", "lof.obs"))
if (is.null(lof_obs_col)) stop("No LoF observed-count column found (lof_obs / lof_hc_lc.obs / lof.obs).")

# Synonymous observed
syn_obs_col <- pick_col(nms, c("syn_obs", "syn.obs", "syn_hc_lc.obs", "synonymous.obs", "synonymous_obs"))
if (is.null(syn_obs_col)) {
  stop("No synonymous observed column found. Checked: syn_obs, syn.obs, syn_hc_lc.obs, synonymous.obs, synonymous_obs")
}

# Missense o/e point estimate (preferred)
mis_oe_col <- pick_col(nms, c("mis.oe", "missense.oe", "mis_oe", "missense_oe"))
# Missense o/e upper CI (optional)
mis_oe_upper_col <- pick_col(nms, c("mis.oe_ci.upper", "missense.oe_ci.upper", "mis_oe_ci_upper", "missense_oe_ci_upper"))

# Missense observed (fallback if missense o/e absent)
mis_obs_col <- pick_col(nms, c("mis_obs", "mis.obs", "mis_hc_lc.obs", "missense.obs", "missense_obs"))

message("[INFO] using columns:")
message("  LOEUF: lof.oe_ci.upper")
message("  LoF obs: ", lof_obs_col)
message("  Syn obs: ", syn_obs_col)
message("  Missense oe: ", ifelse(is.null(mis_oe_col), "<none>", mis_oe_col))
message("  Missense oe upper: ", ifelse(is.null(mis_oe_upper_col), "<none>", mis_oe_upper_col))
message("  Missense obs fallback: ", ifelse(is.null(mis_obs_col), "<none>", mis_obs_col))

if (is.null(mis_oe_col) && is.null(mis_obs_col)) {
  stop("Neither missense o/e nor missense observed columns were found; cannot build Panel B.")
}

# ----------------------------
# Collapse to one row per gene (min LOEUF) to avoid transcript multiplicity artefacts
# ----------------------------
d <- df_raw %>%
  transmute(
    Gene  = norm_sym(gene),
    LOEUF = suppressWarnings(as.numeric(`lof.oe_ci.upper`)),
    lof_obs = suppressWarnings(as.numeric(.data[[lof_obs_col]])),
    syn_obs = suppressWarnings(as.numeric(.data[[syn_obs_col]])),
    mis_oe  = if (!is.null(mis_oe_col)) suppressWarnings(as.numeric(.data[[mis_oe_col]])) else NA_real_,
    mis_oe_upper = if (!is.null(mis_oe_upper_col)) suppressWarnings(as.numeric(.data[[mis_oe_upper_col]])) else NA_real_,
    mis_obs = if (!is.null(mis_obs_col)) suppressWarnings(as.numeric(.data[[mis_obs_col]])) else NA_real_
  ) %>%
  filter(is.finite(LOEUF), is.finite(lof_obs), lof_obs >= 0, is.finite(syn_obs), syn_obs >= 0) %>%
  arrange(Gene, LOEUF) %>%
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

stopifnot(nrow(d) == dplyr::n_distinct(d$Gene))

# Restrict to LOEUF<0.2 universe (this is the universe your discordance is defined within)
u <- d %>%
  filter(LOEUF < LOEUF_CUTOFF) %>%
  mutate(discordant = Gene %in% disc_genes)

message("[QC] universe LOEUF<", LOEUF_CUTOFF, ": n=", nrow(u))
message("[QC] discordant matched in universe: n=", sum(u$discordant))

# QC summary table (helpful for rebutting reviewers)
qc <- u %>%
  summarise(
    n_universe = n(),
    n_discordant = sum(discordant),
    syn_median_discordant = median(syn_obs[discordant], na.rm=TRUE),
    syn_median_other      = median(syn_obs[!discordant], na.rm=TRUE),
    lof_median_discordant = median(lof_obs[discordant], na.rm=TRUE),
    lof_median_other      = median(lof_obs[!discordant], na.rm=TRUE),
    # Wilcoxon tests:
    p_syn = suppressWarnings(wilcox.test(syn_obs ~ discordant, exact = FALSE)$p.value),
    p_lof = suppressWarnings(wilcox.test(lof_obs ~ discordant, exact = FALSE)$p.value)
  )

write_csv(qc, OUT_QC)
message("[DONE] wrote: ", OUT_QC)
print(qc)

# ----------------------------
# Panel A: LOEUF vs synonymous observed
# (If discordant genes are just "highly variable genes", they'd pop here too.)
# ----------------------------
pA <- ggplot(u, aes(x = LOEUF, y = syn_obs)) +
  geom_point(alpha = 0.30, size = 0.9, color = "grey75") +
  geom_point(data = u %>% filter(discordant),
             size = 2.2, shape = 21, fill = "red3", color = "black", stroke = 0.35) +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, LOEUF_CUTOFF)) +
  labs(
    title = "A. Synonymous variation is not elevated in discordant genes",
    subtitle = paste0("Universe: LOEUF < ", LOEUF_CUTOFF,
                      " | Wilcoxon p (syn obs) = ", format(qc$p_syn, digits = 3)),
    x = "LOEUF (LoF o/e upper bound; lower = stronger constraint)",
    y = "Observed synonymous variants per gene (log10)"
  ) +
  theme_classic(base_size = BASE_SIZE) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------
# Panel B: Missense constraint vs truncating burden
# Prefer missense o/e (point estimate). If missing, fall back to missense observed.
# ----------------------------
set.seed(SEED)

if (!all(is.na(u$mis_oe))) {
  # Use missense o/e (point estimate); optionally use upper CI if you prefer
  pB <- ggplot(u, aes(x = mis_oe, y = lof_obs)) +
    geom_point(alpha = 0.30, size = 0.9, color = "grey75") +
    geom_point(data = u %>% filter(discordant),
               size = 2.2, shape = 21, fill = "red3", color = "black", stroke = 0.35) +
    scale_y_log10() +
    labs(
      title = "B. Discordance is not a generic 'high variation' signal",
      subtitle = "Missense constraint (x) vs truncating burden (y); LOEUF<0.2 universe",
      x = "Missense observed/expected (mis.oe; higher = more tolerated)",
      y = "Observed truncating variants per gene (log10)"
    ) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
} else {
  # Fallback: missense observed (still a useful negative control)
  u2 <- u %>% filter(is.finite(mis_obs), mis_obs >= 0)
  pB <- ggplot(u2, aes(x = mis_obs, y = lof_obs)) +
    geom_point(alpha = 0.30, size = 0.9, color = "grey75") +
    geom_point(data = u2 %>% filter(discordant),
               size = 2.2, shape = 21, fill = "red3", color = "black", stroke = 0.35) +
    scale_x_log10() +
    scale_y_log10() +
    labs(
      title = "B. Discordance is not a generic 'high variation' signal",
      subtitle = "Missense observed burden (x) vs truncating burden (y); LOEUF<0.2 universe",
      x = "Observed missense variants per gene (log10)",
      y = "Observed truncating variants per gene (log10)"
    ) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

# ----------------------------
# Combine + save
# ----------------------------
main_title <- "Negative control: discordance is specific to truncation biology"
sub_title  <- paste0("Discordant genes (red) are defined within LOEUF<", LOEUF_CUTOFF,
                     "; synonymous burden provides a neutral comparison.")

combined <- (pA + pB) +
  patchwork::plot_layout(ncol = 2, widths = c(1, 1)) +
  patchwork::plot_annotation(title = main_title, subtitle = sub_title)

ggsave(OUT_PNG, combined, width = 13.2, height = 5.4, units = "in", dpi = 300)
ggsave(OUT_PDF, combined, width = 13.2, height = 5.4, units = "in")

message("[DONE] wrote: ", OUT_PNG)
message("[DONE] wrote: ", OUT_PDF)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Canonical vs minimum-LOEUF transcript comparison (Criticism E)
#
# x = canonical transcript LOEUF
# y = minimum LOEUF across all transcripts for that gene
# highlight discordant genes
#
# Robust to canonical flag encoding (logical / numeric / character).
# Falls back to mane_select if canonical yields zero TRUE rows.
#
# Outputs:
# - SuppFig_canonical_vs_min_LOEUF.png/.pdf
# - canonical_vs_min_LOEUF_table.csv
# - canonical_vs_min_LOEUF_qc.csv
# ============================================================

CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"
DISCORDANT_RDS <- "discordant_genes.rds"

OUT_PNG <- "SuppFig_canonical_vs_min_LOEUF.png"
OUT_PDF <- "SuppFig_canonical_vs_min_LOEUF.pdf"
OUT_CSV <- "canonical_vs_min_LOEUF_table.csv"
OUT_QC  <- "canonical_vs_min_LOEUF_qc.csv"

LOEUF_REF_CUTOFF <- 0.20
BASE_SIZE <- 14

norm_sym <- function(x) toupper(trimws(as.character(x)))

stopifnot(file.exists(CONSTRAINT_TSV))
stopifnot(file.exists(DISCORDANT_RDS))

# ----------------------------
# Load discordant genes
# ----------------------------
disc <- readRDS(DISCORDANT_RDS)
stopifnot("Gene" %in% names(disc))
disc_genes <- disc %>%
  transmute(Gene = norm_sym(Gene)) %>%
  distinct() %>%
  pull(Gene)

message("[INFO] discordant genes loaded: n=", length(disc_genes))

# ----------------------------
# Load constraint table
# ----------------------------
df_raw <- read_tsv(CONSTRAINT_TSV, show_col_types = FALSE)

stopifnot("gene" %in% names(df_raw))
stopifnot("lof.oe_ci.upper" %in% names(df_raw))

# ----------------------------
# Robust canonical-flag parser
# ----------------------------
parse_flag <- function(x) {
  # returns logical vector
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  if (is.numeric(x) || is.integer(x)) return(ifelse(is.na(x), FALSE, x != 0))
  if (is.character(x)) {
    y <- tolower(trimws(x))
    return(y %in% c("true","t","1","yes","y"))
  }
  # fallback for other types
  return(rep(FALSE, length(x)))
}

# Try canonical first; if none flagged, try mane_select
flag_cols <- c("canonical", "mane_select")
flag_cols <- flag_cols[flag_cols %in% names(df_raw)]
if (length(flag_cols) == 0) {
  stop("No 'canonical' or 'mane_select' column found in constraint table.")
}

chosen_flag <- NULL
flag_counts <- list()

for (fc in flag_cols) {
  raw_vals <- df_raw[[fc]]
  parsed <- parse_flag(raw_vals)
  n_true <- sum(parsed, na.rm = TRUE)
  flag_counts[[fc]] <- list(
    n_true = n_true,
    uniq_raw = head(sort(unique(as.character(raw_vals))), 30)
  )
  if (is.null(chosen_flag) && n_true > 0) chosen_flag <- fc
}

# Write QC about flag parsing (so you can debug file formats)
qc_flags <- tibble(
  flag_col = names(flag_counts),
  n_true = vapply(flag_counts, function(z) z$n_true, numeric(1)),
  example_raw_values = vapply(flag_counts, function(z) paste(z$uniq_raw, collapse = " | "), character(1))
)
write_csv(qc_flags, OUT_QC)
message("[DONE] wrote: ", OUT_QC)
print(qc_flags)

if (is.null(chosen_flag)) {
  stop(
    "Could not identify any canonical transcripts: parsed 0 TRUE values for canonical/mane_select.\n",
    "Check OUT_QC for the raw values seen in those columns."
  )
}

message("[INFO] using transcript flag column: ", chosen_flag)

# ----------------------------
# Build transcript-level table
# ----------------------------
tx <- df_raw %>%
  transmute(
    Gene = norm_sym(gene),
    LOEUF = suppressWarnings(as.numeric(`lof.oe_ci.upper`)),
    is_canonical = parse_flag(.data[[chosen_flag]])
  ) %>%
  filter(!is.na(Gene), is.finite(LOEUF))

# ----------------------------
# Per-gene summaries
# ----------------------------
min_tbl <- tx %>%
  group_by(Gene) %>%
  summarise(min_LOEUF = min(LOEUF, na.rm = TRUE), .groups = "drop")

canon_tbl <- tx %>%
  filter(is_canonical) %>%
  group_by(Gene) %>%
  # if multiple flagged, keep minimum among flagged
  summarise(canonical_LOEUF = min(LOEUF, na.rm = TRUE), .groups = "drop")

g <- min_tbl %>%
  inner_join(canon_tbl, by = "Gene") %>%
  mutate(discordant = Gene %in% disc_genes)

# Hard check so we never silently draw an empty plot again
if (nrow(g) == 0) {
  stop(
    "After joining canonical + min, there are 0 genes. This should not happen if canonical flags exist.\n",
    "Please inspect OUT_QC and the constraint TSV schema."
  )
}

message("[QC] genes with canonical+min available: n=", nrow(g))
message("[QC] discordant present in plot table: n=", sum(g$discordant))

write_csv(g %>% arrange(desc(discordant), Gene), OUT_CSV)
message("[DONE] wrote: ", OUT_CSV)

# ----------------------------
# Plot
# ----------------------------
max_lim <- max(c(g$canonical_LOEUF, g$min_LOEUF, LOEUF_REF_CUTOFF), na.rm = TRUE)

p <- ggplot(g, aes(x = canonical_LOEUF, y = min_LOEUF)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7) +
  geom_point(alpha = 0.35, size = 1.0, color = "grey70") +
  geom_point(
    data = g %>% filter(discordant),
    size = 2.3, shape = 21, fill = "red3", color = "black", stroke = 0.35
  ) +
  # Reference lines at LOEUF 0.2 (optional but helpful)
  geom_vline(xintercept = LOEUF_REF_CUTOFF, linetype = "dotted", alpha = 0.6) +
  geom_hline(yintercept = LOEUF_REF_CUTOFF, linetype = "dotted", alpha = 0.6) +
  coord_cartesian(xlim = c(0, max_lim), ylim = c(0, max_lim)) +
  labs(
    title = "Canonical vs minimum-transcript LOEUF per gene",
    subtitle = paste0(
      "Flag column: ", chosen_flag,
      " | dashed line = y = x | dotted = LOEUF ", LOEUF_REF_CUTOFF,
      " | discordant genes highlighted (red)"
    ),
    x = "Canonical transcript LOEUF (lof.oe_ci.upper)",
    y = "Minimum LOEUF across transcripts (lof.oe_ci.upper)"
  ) +
  theme_classic(base_size = BASE_SIZE) + 
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(OUT_PNG, p, width = 7.2, height = 6.3, units = "in", dpi = 300)
ggsave(OUT_PDF, p, width = 7.2, height = 6.3, units = "in")

message("[DONE] wrote: ", OUT_PNG)
message("[DONE] wrote: ", OUT_PDF)



#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

# ------------------------------------------------------------
# Figure: LOEUF (x) vs truncation burden (y), highlight discordant genes
# Fixes ARID1B spam by enforcing ONE POINT PER GENE before plotting/labeling
# Inputs:
#   - gnomad.v4.1.constraint_metrics.tsv
#   - gnomad_variants_all.csv
#   - discordant_genes.rds
# Outputs:
#   - Figure_LOEUF_vs_trunc_variants_discordant_highlight.pdf/png
#   - Figure_LOEUF_vs_trunc_AC_discordant_highlight.pdf/png (optional) - potential for figure 1
# ------------------------------------------------------------

CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"
VARIANTS_CSV   <- "gnomad_variants_all.csv"
DISCORD_RDS    <- "discordant_genes.rds"

# Set to TRUE if you want to zoom to LOEUF < 0.2 (your ultra-constrained universe)
ZOOM_ULTRA <- FALSE
ULTRA_CUTOFF <- 0.2

# How many discordant genes to label (top by trunc burden)
N_LABEL <- 10

# -------------------------
# Load constraint metrics
# -------------------------
cm <- readr::read_tsv(CONSTRAINT_TSV, show_col_types = FALSE)

# Auto-detect gene symbol + LOEUF columns
gene_col <- grep("^gene$|gene_name|gene_symbol|symbol", names(cm), ignore.case = TRUE, value = TRUE)[1]
loeuf_col <- grep("loeuf|lof.*oe.*(ci_)?upper|oe.*lof.*upper|lof_oe_ci_upper|oe_lof_upper|lof_oe_upper",
                  names(cm), ignore.case = TRUE, value = TRUE)[1]

if (is.na(gene_col) || !nzchar(gene_col)) {
  stop("Could not find a gene symbol column in constraint metrics. Check names(cm).")
}
if (is.na(loeuf_col) || !nzchar(loeuf_col)) {
  stop("Could not find LOEUF column in constraint metrics. Check names(cm). Look for lof_oe_ci_upper or similar.")
}

message("[USING] gene_col=", gene_col, " | loeuf_col=", loeuf_col)

gene_tbl <- cm %>%
  transmute(
    gene_symbol = as.character(.data[[gene_col]]),
    loeuf = suppressWarnings(as.numeric(.data[[loeuf_col]]))
  ) %>%
  filter(!is.na(gene_symbol), nzchar(gene_symbol), is.finite(loeuf))

# -------------------------
# Load variants + compute trunc burden per gene
# -------------------------
v <- readr::read_csv(VARIANTS_CSV, show_col_types = FALSE)

gene_trunc <- v %>%
  filter(grepl("stop_gained|frameshift_variant", tolower(consequence))) %>%
  group_by(gene_symbol) %>%
  summarise(
    n_trunc  = n_distinct(variant_id),
    trunc_AC = sum(coalesce(exome_ac, 0L) + coalesce(genome_ac, 0L), na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# Load discordant gene list
# -------------------------
discordant <- readRDS(DISCORD_RDS) %>%
  transmute(gene_symbol = as.character(Gene)) %>%
  distinct()

# -------------------------
# Combine gene-level table
# (one row per gene!)
# -------------------------
plot_df <- gene_tbl %>%
  left_join(gene_trunc, by = "gene_symbol") %>%
  mutate(
    n_trunc  = replace_na(n_trunc, 0L),
    trunc_AC = replace_na(trunc_AC, 0),
    is_discordant = gene_symbol %in% discordant$gene_symbol
  ) %>%
  # Ensure uniqueness: one row per gene_symbol
  group_by(gene_symbol) %>%
  summarise(
    loeuf = first(loeuf),
    n_trunc = max(n_trunc),
    trunc_AC = max(trunc_AC),
    is_discordant = any(is_discordant),
    .groups = "drop"
  )

if (ZOOM_ULTRA) {
  plot_df <- plot_df %>% filter(loeuf < ULTRA_CUTOFF)
}

# -------------------------
# Label top discordant genes only
# -------------------------
label_df <- plot_df %>%
  filter(is_discordant) %>%
  arrange(desc(n_trunc)) %>%
  slice_head(n = N_LABEL)

# -------------------------
# Plot 1: LOEUF vs truncating variant count
# -------------------------
p1 <- ggplot(plot_df, aes(x = loeuf, y = n_trunc)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_point(data = subset(plot_df, is_discordant), size = 2.2) +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "LOEUF (LoF o/e upper bound; lower = more constrained)",
    y = "Truncating variants per gene (log10)",
    title = "Constraint and truncation burden are continuous; discordant genes highlighted"
  ) +
  theme_classic(base_size = 14) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25
  )

ggsave("Figure_LOEUF_vs_trunc_variants_discordant_highlight.pdf", p1, width = 7.6, height = 5.2)
ggsave("Figure_LOEUF_vs_trunc_variants_discordant_highlight.png", p1, width = 7.6, height = 5.2, dpi = 300)

# -------------------------
# Plot 2: LOEUF vs truncating allele count (optional companion)
# -------------------------
label_df2 <- plot_df %>%
  filter(is_discordant) %>%
  arrange(desc(trunc_AC)) %>%
  slice_head(n = N_LABEL)

p2 <- ggplot(plot_df, aes(x = loeuf, y = trunc_AC)) +
  geom_point(alpha = 0.25, size = 1) +
  geom_point(data = subset(plot_df, is_discordant), size = 2.2) +
  scale_y_continuous(trans = "log10") +
  labs(
    x = "LOEUF (LoF o/e upper bound; lower = more constrained)",
    y = "Total truncating AC per gene (log10)",
    title = "Constraint and truncating allele burden are continuous; discordant genes highlighted"
  ) +
  theme_classic(base_size = 14) +
  ggrepel::geom_text_repel(
    data = label_df2,
    aes(label = gene_symbol),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25
  )

ggsave("Figure_LOEUF_vs_trunc_AC_discordant_highlight.pdf", p2, width = 7.6, height = 5.2)
ggsave("Figure_LOEUF_vs_trunc_AC_discordant_highlight.png", p2, width = 7.6, height = 5.2, dpi = 300)

message("[DONE] wrote:\n  Figure_LOEUF_vs_trunc_variants_discordant_highlight.(pdf/png)\n  Figure_LOEUF_vs_trunc_AC_discordant_highlight.(pdf/png)")

print(p1)
