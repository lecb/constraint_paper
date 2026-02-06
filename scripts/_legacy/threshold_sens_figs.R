#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

OUTDIR <- "gnomad_lof_discordance_out"
TABDIR <- "tables"
CACHEDIR <- file.path(OUTDIR, "cache")

GF_RDS    <- file.path(CACHEDIR, "gene_features.rds")
SENS_FILE <- file.path(OUTDIR, "threshold_sensitivity_summary.csv")
DISC_FILE <- file.path(TABDIR, "Supplementary_Table1_discordant_genes.csv")

stopifnot(file.exists(GF_RDS), file.exists(SENS_FILE), file.exists(DISC_FILE))

# Optional: patchwork for multi-panel layout
HAS_PATCHWORK <- requireNamespace("patchwork", quietly = TRUE)
if (HAS_PATCHWORK) suppressPackageStartupMessages(library(patchwork))

# ----------------------------
# Load inputs
# ----------------------------
gf <- readRDS(GF_RDS)

# Expect: gene, loeuf, oe_lof, ci_width, syn_obs, trunc_obs (as in your analysis script)
need <- c("gene","loeuf","oe_lof","ci_width","syn_obs","trunc_obs")
miss <- setdiff(need, names(gf))
if (length(miss) > 0) stop("gene_features.rds missing columns: ", paste(miss, collapse=", "))

disc_genes <- readr::read_csv(DISC_FILE, show_col_types = FALSE) %>%
  pull(Gene) %>%
  toupper() %>%
  unique()

sens <- readr::read_csv(SENS_FILE, show_col_types = FALSE) %>%
  mutate(run_id = row_number())

# ============================================================
# FIGURE 1: SuppFig_robustness_clean (A/B)
# ============================================================

# ----------------------------
# Panel A: threshold robustness
# tight x-range + show x=0 null in left margin
# ----------------------------

sensA <- sens %>%
  arrange(med_resid_diff) %>%
  mutate(run_id = row_number())

# Tight x-limits around observed values (no whitespace)
xmin <- min(sensA$med_resid_diff, na.rm = TRUE) - 0.02
xmax <- max(sensA$med_resid_diff, na.rm = TRUE) + 0.02

pA <- ggplot(sensA, aes(x = med_resid_diff, y = run_id)) +
  geom_point(size = 2.8) +
  
  # True null line at x = 0 (will be outside panel, drawn into margin)
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
  
  # Label for the null line (also outside panel)
  annotate(
    "text",
    x = 0,
    y = max(sensA$run_id) + 1,
    label = "0 (null)",
    size = 3,
    hjust = 0
  ) +
  
  labs(
    x = "Median residual truncation difference (discordant − background)",
    y = "Threshold parameter set"
  ) +
  
  # Keep tight x-range BUT allow drawing outside the panel
  coord_cartesian(xlim = c(xmin, xmax), clip = "off") +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    # Add left/top margin so the off-panel line/label are visible
    plot.margin = margin(t = 8, r = 8, b = 8, l = 28)
  )


# Panel B: continuous score within LOEUF<0.2 + precision filter (match main)
DISC_LOEUF_MAX    <- 0.20
DISC_DROP_CI_FRAC <- 0.25
MIN_SYN           <- 50

d_score <- gf %>%
  filter(!is.na(gene), !is.na(loeuf), !is.na(ci_width), !is.na(trunc_obs), !is.na(syn_obs)) %>%
  filter(loeuf < DISC_LOEUF_MAX, syn_obs >= MIN_SYN)

ci_thr <- quantile(d_score$ci_width, probs = 1 - DISC_DROP_CI_FRAC, na.rm = TRUE)
d_score <- d_score %>% filter(ci_width <= ci_thr)

fit_trunc <- lm(log1p(trunc_obs) ~ log1p(syn_obs), data = d_score)

d_score <- d_score %>%
  mutate(
    disc_score = rstudent(fit_trunc),
    group = ifelse(toupper(gene) %in% disc_genes, "Discordant", "Other LOEUF<0.2")
  )

pB <- ggplot(d_score, aes(x = group, y = disc_score)) +
  geom_violin(trim = TRUE, alpha = 0.25) +
  geom_boxplot(width = 0.2, outlier.size = 0.7) +
  labs(x = NULL, y = "Studentised truncation-excess score") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

out1_pdf <- file.path(OUTDIR, "SuppFig_robustness_clean.pdf")
out1_png <- file.path(OUTDIR, "SuppFig_robustness_clean.png")

if (HAS_PATCHWORK) {
  combo1 <- pA / pB + plot_annotation(tag_levels = "A")
  ggsave(out1_pdf, combo1, width = 9, height = 7.5, dpi = 300)
  ggsave(out1_png, combo1, width = 9, height = 7.5, dpi = 300)
} else {
  ggsave(file.path(OUTDIR, "SuppFig_robustness_clean_A.pdf"), pA, width = 9, height = 3.6, dpi = 300)
  ggsave(file.path(OUTDIR, "SuppFig_robustness_clean_B.pdf"), pB, width = 9, height = 3.6, dpi = 300)
}

message("[WROTE] ", out1_pdf)
message("[WROTE] ", out1_png)

# ============================================================
# FIGURE 2: Missing figure merged into one output
#   SuppFig_continuous_discordance_AC_only_top35 (A/C)
#   (LoF o/e ECDF + ranked tail), based on gene_features.rds
# ============================================================

CONSTRAINT_TSV <- "gnomad.v4.1.constraint_metrics.tsv"
stopifnot(file.exists(CONSTRAINT_TSV))

LOEUF_CUTOFF <- 0.2
DROP_WORST_CI_WIDTH_FRAC <- 0.25
TOP_N <- 35L
WRITE_PDF <- TRUE

raw <- readr::read_tsv(CONSTRAINT_TSV, show_col_types = FALSE, progress = FALSE)

required_cols <- c("gene","lof.oe","lof.oe_ci.lower","lof.oe_ci.upper","syn.obs","cds_length")
missing_cols <- setdiff(required_cols, names(raw))
if (length(missing_cols) > 0) stop("Missing required columns in TSV: ", paste(missing_cols, collapse = ", "))

gene_level <- raw %>%
  group_by(gene) %>%
  summarise(
    LOEUF    = if (all(is.na(`lof.oe_ci.upper`))) NA_real_ else min(`lof.oe_ci.upper`, na.rm = TRUE),
    lof_oe   = if (all(is.na(`lof.oe`)))         NA_real_ else min(`lof.oe`,         na.rm = TRUE),
    ci_lower = if (all(is.na(`lof.oe_ci.lower`))) NA_real_ else min(`lof.oe_ci.lower`, na.rm = TRUE),
    ci_upper = if (all(is.na(`lof.oe_ci.upper`))) NA_real_ else min(`lof.oe_ci.upper`, na.rm = TRUE),
    ci_width = ci_upper - ci_lower,
    syn_obs  = if (all(is.na(`syn.obs`)))        NA_real_ else max(`syn.obs`,        na.rm = TRUE),
    cds_length = if (all(is.na(cds_length)))     NA_real_ else max(cds_length,       na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(LOEUF), is.finite(lof_oe), is.finite(ci_width))

uc <- gene_level %>% filter(LOEUF < LOEUF_CUTOFF)
if (nrow(uc) < TOP_N) stop("LOEUF<", LOEUF_CUTOFF, " yields fewer genes (", nrow(uc), ") than TOP_N=", TOP_N)

ci_cut <- quantile(uc$ci_width, probs = 1 - DROP_WORST_CI_WIDTH_FRAC, na.rm = TRUE)
uc_precise <- uc %>%
  filter(ci_width <= ci_cut) %>%
  arrange(desc(lof_oe))

topN_cont <- uc_precise %>% slice_head(n = TOP_N)
cutoff_oe <- min(topN_cont$lof_oe, na.rm = TRUE)

# --- Panel A: ECDF with cutoff annotation (match old style) ---
pA2 <- ggplot(uc_precise, aes(x = lof_oe)) +
  stat_ecdf() +
  geom_vline(xintercept = cutoff_oe, linetype = "dashed") +
  annotate("text",
           x = cutoff_oe + 0.005, y = 0.88,
           label = paste0("Top ", TOP_N, " cutoff = ", signif(cutoff_oe, 3)),
           hjust = 0, size = 3.2) +
  labs(
    x = "LoF observed/expected (point estimate)",
    y = "ECDF"
  ) +
  theme_minimal(base_size = 12)

# --- Panel C: rank curve with topN highlighted (match old style) ---
uc_rank <- uc_precise %>% mutate(rank = row_number())

pC2 <- ggplot(uc_rank, aes(x = rank, y = lof_oe)) +
  geom_line() +
  geom_point(data = uc_rank %>% slice_head(n = TOP_N),
             colour = "red", size = 2) +
  geom_vline(xintercept = TOP_N, linetype = "dashed") +
  annotate("text",
           x = TOP_N + 5, y = cutoff_oe + 0.005,
           label = paste0("Top ", TOP_N, " cutoff = ", signif(cutoff_oe, 3)),
           hjust = 0, size = 3.2) +
  labs(
    x = paste0("Rank within LOEUF < ", LOEUF_CUTOFF, " genes (1 = highest LoF o/e)"),
    y = "LoF observed/expected"
  ) +
  theme_minimal(base_size = 12)

out2_pdf <- file.path(OUTDIR, "SuppFig_continuous_discordance_AC_only_top35.pdf")
out2_png <- file.path(OUTDIR, "SuppFig_continuous_discordance_AC_only_top35.png")

if (HAS_PATCHWORK) {
  combo2 <- pA2 / pC2 + plot_annotation(tag_levels = "A")
  ggsave(out2_pdf, combo2, width = 8.5, height = 7.5, dpi = 300)
  ggsave(out2_png, combo2, width = 8.5, height = 7.5, dpi = 300)
} else {
  ggsave(file.path(OUTDIR, "SuppFig_continuous_discordance_AC_only_top35_A.pdf"), pA2, width = 8.5, height = 3.6, dpi = 300)
  ggsave(file.path(OUTDIR, "SuppFig_continuous_discordance_AC_only_top35_C.pdf"), pC2, width = 8.5, height = 3.6, dpi = 300)
}

message("[WROTE] ", out2_pdf)
message("[WROTE] ", out2_png)
# ============================================================
# FIGURE 3: SuppFig_discScore_ECDF_genome_vs_loeuf_vs_discordant
# ============================================================

# ----------------------------
# Fig 3: ECDF of truncation-excess score (genome vs LOEUF<0.2 vs discordant)
# ----------------------------

# ----------------------------
# Fig 3: ECDF of truncation-excess score (genome vs LOEUF<0.2 vs discordant)
# ----------------------------

d_all <- gf %>%
  filter(!is.na(gene), !is.na(loeuf), !is.na(trunc_obs), !is.na(syn_obs)) %>%
  mutate(gene_u = toupper(gene))

fit_all <- lm(log1p(trunc_obs) ~ log1p(syn_obs), data = d_all)

# compute discordant n once (for stable label)
disc_u <- toupper(disc_genes)
n_disc <- sum(d_all$gene_u %in% disc_u, na.rm = TRUE)
disc_label <- paste0("Discordant (n=", n_disc, ")")

d_all <- d_all %>%
  mutate(
    trunc_excess_rstud = rstudent(fit_all),
    disc_score = trunc_excess_rstud,  # higher = more excess truncation
    is_discordant = gene_u %in% disc_u,
    group = case_when(
      is_discordant ~ disc_label,
      loeuf < 0.2   ~ "Other LOEUF<0.2",
      TRUE          ~ "Genome-wide"
    ),
    group = factor(group, levels = c("Genome-wide", "Other LOEUF<0.2", disc_label))
  )

# Red palette (light -> dark), built safely
red_pal <- setNames(
  c("#FAD4D4", "#E45757", "#8B0000"),
  c("Genome-wide", "Other LOEUF<0.2", disc_label)
)

# Vertical dashed line at the discordant threshold (min discordant score)
disc_cut <- d_all %>%
  filter(is_discordant) %>%
  summarise(cut = min(disc_score, na.rm = TRUE)) %>%
  pull(cut)

p3 <- ggplot(d_all, aes(x = disc_score, colour = group)) +
  stat_ecdf(linewidth = 1.15) +
  
  # Discordant cutoff (matches discordant colour)
  geom_vline(
    xintercept = disc_cut,
    linetype = "dashed",
    linewidth = 0.8,
    colour = red_pal[[disc_label]],
    alpha = 0.9
  ) +
  
  # Rug marks for discordant genes only (matches discordant colour)
  geom_rug(
    data = d_all %>% filter(is_discordant),
    aes(x = disc_score),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.35,
    linewidth = 0.35,
    colour = red_pal[[disc_label]]
  ) +
  
  scale_colour_manual(values = red_pal, drop = FALSE) +
  
  labs(
    x = "Studentised truncation-excess score",
    y = "ECDF",
    colour = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.78, 0.25),
    legend.background = element_blank(),
    legend.key = element_blank(),
    axis.line = element_line(linewidth = 0.6)
  )

out3_pdf <- file.path(OUTDIR, "SuppFig_discScore_ECDF_genome_vs_loeuf_vs_discordant.pdf")
out3_png <- file.path(OUTDIR, "SuppFig_discScore_ECDF_genome_vs_loeuf_vs_discordant.png")

ggsave(out3_pdf, p3, width = 10, height = 6)
ggsave(out3_png, p3, width = 10, height = 6, dpi = 300)

message("[WROTE] ", out3_pdf)
message("[WROTE] ", out3_png)

# ============================================================
# SUPER FIGURE (A–E): panel tags OUTSIDE plot area; Panel E legend inside E
# Assumes these exist: pA, pB, pA2, pC2, p3
# ============================================================

stopifnot(exists("pA"), exists("pB"), exists("pA2"), exists("pC2"), exists("p3"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# Styling constants
# ----------------------------
base_theme <- theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10)
  )

# Put panel tags OUTSIDE the plotting region (prevents overlap with axes/geoms)
tag_outside_theme <- theme(
  plot.tag = element_text(size = 16, face = "bold"),
  # negative x pushes tag to the left of the panel
  plot.tag.position = c(-0.06, 0.99),
  # give room on the left so the tag doesn't get clipped by device edge
  plot.margin = margin(14, 14, 14, 34)
)

# Helper to apply tag + unclipped drawing region
tag_outside <- function(p, letter) {
  p +
    labs(tag = letter) +
    coord_cartesian(clip = "off") +
    base_theme +
    tag_outside_theme
}

# ----------------------------
# Apply per-panel tags + kill legends in A–D
# ----------------------------
pA_fix <- tag_outside(pA,  "A") + theme(legend.position = "none")
pB_fix <- tag_outside(pB,  "B") + theme(legend.position = "none")
pC_fix <- tag_outside(pA2, "C") + theme(legend.position = "none")
pD_fix <- tag_outside(pC2, "D") + theme(legend.position = "none")

# Panel E: legend INSIDE the plot (so no extra right-hand whitespace)
pE_fix <- tag_outside(p3, "E") +
  theme(
    legend.position = c(0.78, 0.28),        # inside panel E
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 9)
  )

# ----------------------------
# Layout: left column (A over B), right column (C over D over E)
# Add spacing with plot.margin (above) + panel spacing
# ----------------------------
left_col  <- pA_fix / pB_fix + plot_layout(heights = c(1, 1.05))
right_col <- pC_fix / pD_fix / pE_fix + plot_layout(heights = c(0.95, 1.05, 1.15))

super_fig <- left_col | right_col

# IMPORTANT: do NOT collect guides, otherwise legend can drift
super_fig <- super_fig + plot_layout(guides = "keep")

# Slightly increase whitespace between panels globally
super_fig <- super_fig & theme(panel.spacing = unit(10, "pt"))

# ----------------------------
# Save big enough for readability
# ----------------------------
OUTDIR <- ifelse(exists("OUTDIR"), OUTDIR, "gnomad_lof_discordance_out")
OUT_SUPER_PDF <- file.path(OUTDIR, "SuppFig_super_discordance_robustness.pdf")
OUT_SUPER_PNG <- file.path(OUTDIR, "SuppFig_super_discordance_robustness.png")

ggsave(OUT_SUPER_PDF, super_fig, width = 14.5, height = 9.5)
ggsave(OUT_SUPER_PNG, super_fig, width = 14.5, height = 9.5, dpi = 300)

message("[WROTE] ", OUT_SUPER_PDF)
message("[WROTE] ", OUT_SUPER_PNG)
