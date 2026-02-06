#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# Plot: LOEUF (lof.oe_ci.upper) vs truncating variants per gene
#
# CLEAN VERSION (reviewer-proof):
# - One explicit selection cutoff line at LOEUF = 0.20
# - NO grey shaded LOEUF band (removes threshold confusion)
# - Grey background genes shown across LOEUF 0–0.30 for continuity
# - Red points + labels = exactly Table 1 Top N genes
#   (read from Table1_top_discordant_genes_QC_skeleton.csv)
# - One row per gene for background and red points:
#   pick the minimum LOEUF row per gene to avoid transcript mismatches
# - Smooth expectation fitted on non-Table1 genes only, restricted to LOEUF<=0.5
# - Horizontal dashed line = median truncating burden among non-Table1 genes with LOEUF<=0.20
# - No arrows, no free text on plot
# ============================================================

CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"
TABLE1_CSV     <- "gnomad_lof_discordance_out/Table1_top_discordant_genes_QC_skeleton.csv"

OUT_PDF <- "Figure_LOEUF_vs_truncating_clean_Table1Top15label.pdf"
OUT_PNG <- "Figure_LOEUF_vs_truncating_clean_Table1Top15label.png"

# ---- Parameters ----
TABLE1_N             <- 15
DISCORDANT_LOEUF_MAX <- 0.20
X_LIM                <- c(0, 0.30)   # keep context; set NULL for full range
SMOOTH_MAX_LOEUF     <- 0.50
USE_PSEUDOCOUNT      <- TRUE
LOESS_SPAN           <- 0.75

norm_sym <- function(x) toupper(trimws(as.character(x)))

stopifnot(file.exists(CONSTRAINT_TSV))
stopifnot(file.exists(TABLE1_CSV))

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  stop("Package 'ggrepel' is required for labels. Install with: install.packages('ggrepel')")
}

# ----------------------------
# Load Table 1 genes (truth set)
# ----------------------------
t1 <- read_csv(TABLE1_CSV, show_col_types = FALSE)
stopifnot("Gene" %in% names(t1))

TABLE1_GENES <- t1 %>%
  mutate(Gene = norm_sym(Gene)) %>%
  slice_head(n = TABLE1_N) %>%
  pull(Gene)

message("[INFO] Table1 genes loaded: n=", length(TABLE1_GENES))
message("[INFO] Table1 genes: ", paste(TABLE1_GENES, collapse = ", "))

# ----------------------------
# Load constraint data
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
message("[INFO] using LoF observed column: ", lof_obs_col)

# ----------------------------
# One row per gene, transcript-safe:
# choose the MIN LOEUF row per gene (most constrained).
# This avoids transcript/canonical mismatches and guarantees
# that if a gene had any row with LOEUF<0.2, its gene-level
# representative stays <0.2.
# ----------------------------
df_bg <- df_raw %>%
  transmute(
    Gene    = norm_sym(gene),
    LOEUF   = suppressWarnings(as.numeric(`lof.oe_ci.upper`)),
    lof_obs = suppressWarnings(as.numeric(.data[[lof_obs_col]]))
  ) %>%
  filter(is.finite(LOEUF), is.finite(lof_obs), lof_obs >= 0) %>%
  arrange(Gene, LOEUF) %>%              # min LOEUF first
  group_by(Gene) %>%
  slice(1) %>%
  ungroup()

stopifnot(nrow(df_bg) == dplyr::n_distinct(df_bg$Gene))

# Red/Table1 layer from same representation
df_red <- df_bg %>%
  filter(Gene %in% TABLE1_GENES) %>%
  mutate(discordant = LOEUF < DISCORDANT_LOEUF_MAX)

message("[QC] Table1 LOEUF range (plot rows): ",
        paste(range(df_red$LOEUF, na.rm = TRUE), collapse = " - "))

# Enforce your definition: Table1 genes must all be LOEUF<0.2
stopifnot(all(df_red$LOEUF < DISCORDANT_LOEUF_MAX))

# ----------------------------
# Log-safe y
# ----------------------------
if (USE_PSEUDOCOUNT) {
  df_bg  <- df_bg  %>% mutate(y_plot = lof_obs + 1)
  df_red <- df_red %>% mutate(y_plot = lof_obs + 1)
  y_lab <- "Observed truncating variants per gene (log10; +1 pseudocount)"
} else {
  df_bg  <- df_bg  %>% filter(lof_obs > 0) %>% mutate(y_plot = lof_obs)
  df_red <- df_red %>% filter(lof_obs > 0) %>% mutate(y_plot = lof_obs)
  y_lab <- "Observed truncating variants per gene (log10)"
}

# ----------------------------
# Horizontal reference: median truncating burden among OTHER genes with LOEUF<=0.2
# ----------------------------
ref_df <- df_bg %>% filter(LOEUF <= DISCORDANT_LOEUF_MAX, !(Gene %in% TABLE1_GENES))
y_ref_med <- if (nrow(ref_df) > 0) median(ref_df$y_plot, na.rm = TRUE) else NA_real_

# ----------------------------
# Smooth expectation: non-Table1 genes only; restricted LOEUF to avoid tail weirdness
# ----------------------------
smooth_df <- df_bg %>% filter(!(Gene %in% TABLE1_GENES), LOEUF <= SMOOTH_MAX_LOEUF)

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(df_bg, aes(x = LOEUF, y = y_plot)) +
  
  # Explicit selection cutoff (ONLY cutoff shown)
  geom_vline(
    xintercept = DISCORDANT_LOEUF_MAX,
    linewidth = 0.45,
    linetype = "dashed",
    color = "grey55"
  ) +
  
  # Background points
  geom_point(color = "grey85", alpha = 0.35, size = 0.8) +
  
  # Smooth expectation (non-Table1 genes only)
  geom_smooth(
    data = smooth_df,
    aes(x = LOEUF, y = y_plot),
    method = "loess",
    span = LOESS_SPAN,
    se = FALSE,
    color = "grey60",
    linewidth = 0.6
  ) +
  
  # Red points (Table1 genes)
  geom_point(
    data = df_red,
    shape = 21,
    fill = "red3",
    color = "black",
    stroke = 0.35,
    size = 2.6
  ) +
  
  # Horizontal median reference (within LOEUF<=0.2, non-Table1)
  { if (!is.na(y_ref_med))
    geom_hline(yintercept = y_ref_med, linetype = "dashed", linewidth = 0.45)
  } +
  
  # Labels for red points
  ggrepel::geom_text_repel(
    data = df_red,
    aes(label = Gene),
    size = 3.2,
    min.segment.length = 0,
    segment.alpha = 0.5,
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf
  ) +
  
  scale_y_log10() +
  { if (!is.null(X_LIM)) coord_cartesian(xlim = X_LIM) } +
  
  labs(
    title = "Highly constrained genes can still tolerate structured partial truncation",
    x = "LOEUF (lower = stronger constraint vs true LoF)",
    y = y_lab
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(OUT_PDF, p, width = 9, height = 5.6, units = "in")
ggsave(OUT_PNG, p, width = 9, height = 5.6, units = "in", dpi = 300)

message("[DONE] wrote: ", OUT_PDF)
message("[DONE] wrote: ", OUT_PNG)
