#!/usr/bin/env Rscript
# ============================================================
# PEXT (GTEx v10) transcription-awareness analysis for discordant genes
# - Streams huge pext TSV.GZ in chunks (no full load)
# - Extracts pLoF-ish variants for:
#     (1) Discordant35 genes
#     (2) Background LOEUF<0.2 genes excluding Discordant35
# - Writes filtered variant-level TSV.GZ files
# - Builds gene-level summaries + statistical tests
# - Produces manuscript-ready figures (including the 2-panel combined)
# - Writes supplementary tables
#
# Outputs (default):
#   pext/pext_variants_discordant35.tsv.gz
#   pext/pext_variants_background_LOEUF02.tsv.gz
#   tables/SuppTable_pext_variant_level_discordant_vs_bg.tsv.gz
#   tables/SuppTable_pext_gene_level_discordant_vs_bg.tsv.gz
#   tables/pext_summary_stats.txt
#   figures/SuppFig_pext_expression_attenuation_discordant35_2panel.png (+ .pdf)
#   figures/SuppFig_pext_ecdf_discordant_vs_bg.png
#   figures/SuppFig_pext_genelevel_lowpext_violin.png
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ----------------------------
# CONFIG
# ----------------------------
INFILE   <- "pext/gnomad.pext.gtex_v10.annotation_level.tsv.gz"
DISC_RDS <- "cache/discordant_35.rds"
CG_RDS   <- "cache/constraint_gene.rds"

OUT_DISC <- "pext/pext_variants_discordant35.tsv.gz"
OUT_BG   <- "pext/pext_variants_background_LOEUF02.tsv.gz"

FIG_DIR  <- "figures"
TAB_DIR  <- "tables"

# pext threshold for "low expression"
PEXT_CUTOFF <- 0.2

# Define consequences to keep (variant-level pext file is already consequence-annotated)
PLOF_CSQS <- c("stop_gained","frameshift_variant","splice_donor_variant","splice_acceptor_variant")

# Chunk streaming settings
CHUNK_SIZE <- 200000L

# Preflight options (fast sanity checks without full streaming)
RUN_PREFLIGHT_ONLY <- FALSE     # TRUE = do not stream full file
PREFLIGHT_CHUNKS   <- 10L       # number of chunks to read in preflight
PREFLIGHT_SEED     <- 1L

# ----------------------------
# Helpers
# ----------------------------
dir.create(dirname(OUT_DISC), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(OUT_BG),   showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,           showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR,           showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(INFILE))
stopifnot(file.exists(DISC_RDS))
stopifnot(file.exists(CG_RDS))

msg <- function(...) message(sprintf(...))

# Read header only
hdr <- names(read_tsv(INFILE, n_max = 0, show_col_types = FALSE))

needed <- c("locus","alleles","gene_symbol","most_severe_consequence","lof","lof_flags","exp_prop_mean")
missing <- setdiff(needed, hdr)
if (length(missing) > 0) {
  stop("Missing expected columns in pext file: ", paste(missing, collapse = ", "))
}
msg("OK: required pext columns present.")

# Column spec (force exp_prop_mean numeric; avoid the <double>/<character> bind error)
col_spec <- cols(.default = col_skip())
col_spec$cols$locus <- col_character()
col_spec$cols$alleles <- col_character()
col_spec$cols$gene_symbol <- col_character()
col_spec$cols$most_severe_consequence <- col_character()
col_spec$cols$lof <- col_character()
col_spec$cols$lof_flags <- col_character()
col_spec$cols$exp_prop_mean <- col_double()

# ----------------------------
# Load gene lists
# ----------------------------
disc_obj  <- readRDS(DISC_RDS)
disc_genes <- unique(as.character(disc_obj$Gene))
stopifnot(length(disc_genes) == 35)
msg("Discordant genes: %d", length(disc_genes))

constraint_gene <- readRDS(CG_RDS)

# Collapse to 1 row per gene (constraint_gene may be transcript-level)
gene_df <- constraint_gene %>%
  group_by(Gene) %>%
  summarise(
    LOEUF = if (all(is.na(LOEUF))) NA_real_ else min(LOEUF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Gene = as.character(Gene))

uc_genes <- gene_df %>%
  filter(is.finite(LOEUF), LOEUF < 0.2) %>%
  pull(Gene) %>%
  unique()

bg_genes <- setdiff(uc_genes, disc_genes)

msg("LOEUF<0.2 genes: %d ; background genes (excluding discordant): %d", length(uc_genes), length(bg_genes))
stopifnot(length(bg_genes) > 0)

target_all <- unique(c(disc_genes, bg_genes))

# ----------------------------
# Chunk-stream extraction (in-memory collection)
# ----------------------------
collect_matches <- function(n_chunks = Inf) {
  # This returns a single tibble of matched rows with a 'group' column.
  # We collect in a list to avoid slow repeated rbind.
  seen_chunks <- 0L
  out_list <- vector("list", 0)
  
  cb <- function(x, pos) {
    seen_chunks <<- seen_chunks + 1L
    # Filter to target genes + pLoF-ish consequences + non-missing pext
    y <- x %>%
      filter(!is.na(gene_symbol)) %>%
      filter(gene_symbol %in% target_all) %>%
      filter(most_severe_consequence %in% PLOF_CSQS) %>%
      filter(!is.na(exp_prop_mean)) %>%
      transmute(
        locus,
        alleles,
        gene_symbol,
        most_severe_consequence,
        lof,
        lof_flags,
        exp_prop_mean,
        group = if_else(gene_symbol %in% disc_genes, "Discordant35", "Background_LOEUF02")
      )
    
    if (nrow(y) > 0) out_list[[length(out_list) + 1L]] <<- y
    
    if (is.finite(n_chunks) && seen_chunks >= n_chunks) {
      # Stop streaming early
      stop("STOP_AFTER_N_CHUNKS")
    }
    invisible()
  }
  
  msg("Streaming pext in chunks (chunk_size=%d)...", CHUNK_SIZE)
  
  # Use chunked reader; abort early in preflight by throwing an error we catch
  res <- try(
    read_tsv_chunked(
      file = INFILE,
      callback = SideEffectChunkCallback$new(cb),
      chunk_size = CHUNK_SIZE,
      col_types = col_spec,
      progress = interactive()
    ),
    silent = TRUE
  )
  
  if (inherits(res, "try-error")) {
    if (!grepl("STOP_AFTER_N_CHUNKS", as.character(res), fixed = TRUE)) {
      stop(res)
    }
  }
  
  v_all <- bind_rows(out_list)
  
  msg("Chunks read: %d", seen_chunks)
  msg("Matched rows: %d (Discordant=%d ; Background=%d)",
      nrow(v_all),
      sum(v_all$group == "Discordant35"),
      sum(v_all$group == "Background_LOEUF02"))
  
  v_all
}

# ----------------------------
# Preflight (fast check)
# ----------------------------
if (RUN_PREFLIGHT_ONLY) {
  set.seed(PREFLIGHT_SEED)
  v_pf <- collect_matches(n_chunks = PREFLIGHT_CHUNKS)
  
  # quick sanity
  if (!all(c("Discordant35","Background_LOEUF02") %in% unique(v_pf$group))) {
    warning("Preflight did not capture both groups. This can happen by chance; increase PREFLIGHT_CHUNKS.")
  }
  msg("Preflight done. Exiting (RUN_PREFLIGHT_ONLY=TRUE).")
  quit(save = "no", status = 0)
}

# ----------------------------
# FULL extraction
# ----------------------------
v_all <- collect_matches(n_chunks = Inf)

# Split and write filtered variant files
disc_v <- v_all %>% filter(group == "Discordant35") %>%
  select(-group)
bg_v   <- v_all %>% filter(group == "Background_LOEUF02") %>%
  select(-group)

if (nrow(disc_v) == 0) stop("No discordant variant rows were extracted — unexpected. Check gene_symbol matching.")
if (nrow(bg_v) == 0)   stop("No background variant rows were extracted — usually a gene_symbol mismatch. Check LOEUF background gene symbols vs pext gene_symbol.")

write_tsv(disc_v, OUT_DISC)
write_tsv(bg_v,   OUT_BG)
msg("Wrote filtered discordant variants to: %s", OUT_DISC)
msg("Wrote filtered background variants to: %s", OUT_BG)

# ----------------------------
# Build tables + stats
# ----------------------------
# Variant-level summary table (both groups, for supplement)
v_all2 <- v_all %>%
  mutate(
    is_low_pext = exp_prop_mean < PEXT_CUTOFF
  )

# Per-group counts and proportions
variant_group_summary <- v_all2 %>%
  group_by(group) %>%
  summarise(
    n_variants = n(),
    frac_low_pext = mean(is_low_pext),
    median_pext = median(exp_prop_mean),
    mean_pext = mean(exp_prop_mean),
    .groups = "drop"
  )

# Gene-level summary (each gene contributes equally)
gene_pext <- v_all2 %>%
  group_by(group, gene_symbol) %>%
  summarise(
    n_variants = n(),
    prop_low_pext = mean(is_low_pext),
    median_pext = median(exp_prop_mean),
    mean_pext = mean(exp_prop_mean),
    .groups = "drop"
  )

# Stats:
# 1) Variant-level shift (Wilcoxon)
w_variant <- wilcox.test(exp_prop_mean ~ group, data = v_all2)
# 2) Gene-level shift in low-pext proportion (Wilcoxon)
w_gene <- wilcox.test(prop_low_pext ~ group, data = gene_pext)

# 3) Variant-level 2x2 enrichment low vs high (Fisher)
tab_low <- table(v_all2$group, v_all2$is_low_pext)
f_low <- fisher.test(tab_low)

# Write tables
SUPP_VAR_TAB <- file.path(TAB_DIR, "SuppTable_pext_variant_level_discordant_vs_bg.tsv.gz")
SUPP_GENE_TAB <- file.path(TAB_DIR, "SuppTable_pext_gene_level_discordant_vs_bg.tsv.gz")
STATS_TXT <- file.path(TAB_DIR, "pext_summary_stats.txt")

write_tsv(v_all2, SUPP_VAR_TAB)
write_tsv(gene_pext, SUPP_GENE_TAB)

# Write stats text
stats_lines <- c(
  "PEXT transcription-awareness summary",
  "===================================",
  sprintf("PEXT cutoff (low expression): %0.3f", PEXT_CUTOFF),
  "",
  "Variant-level summary by group:",
  paste(capture.output(print(variant_group_summary)), collapse = "\n"),
  "",
  "Variant-level Wilcoxon test: exp_prop_mean ~ group",
  sprintf("  W = %s ; p = %.3g", format(w_variant$statistic), w_variant$p.value),
  "",
  "Gene-level Wilcoxon test: prop_low_pext ~ group",
  sprintf("  W = %s ; p = %.3g", format(w_gene$statistic), w_gene$p.value),
  "",
  "Variant-level Fisher test (low vs not low pext):",
  paste(capture.output(print(tab_low)), collapse = "\n"),
  sprintf("  OR = %.3g ; p = %.3g", unname(f_low$estimate), f_low$p.value)
)
writeLines(stats_lines, STATS_TXT)

msg("Wrote tables:")
msg("  %s", SUPP_VAR_TAB)
msg("  %s", SUPP_GENE_TAB)
msg("Wrote stats:")
msg("  %s", STATS_TXT)

# ============================================================
# Figures (FINAL styling tweaks as requested)
# - Legend moved (panel A) so it doesn’t sit in the way
# - Stats text same size for A and B
# - Low-expression threshold label nudged right (no overlap with dashed line)
# - 2-panel combined saved to PNG + PDF
# ============================================================

# ----------------------------
# Figures
# ----------------------------

# Helper: pretty p-values (so you get "< 1e-16" instead of "0")
p_pretty <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-16) return("< 1e-16")
  formatC(p, format = "g", digits = 3)
}

# Panel-level stats strings
p_variant_txt <- paste0("Variant-level Wilcoxon p = ", p_pretty(w_variant$p.value))
p_gene_txt    <- paste0("Gene-level Wilcoxon p = ", p_pretty(w_gene$p.value))

# Low-pext fractions (variant-level) for Panel A subtitle
frac_disc <- mean(v_all2$exp_prop_mean[v_all2$group == "Discordant35"] < PEXT_CUTOFF)
frac_bg   <- mean(v_all2$exp_prop_mean[v_all2$group == "Background_LOEUF02"] < PEXT_CUTOFF)

# Gene-level medians + n genes for Panel B subtitle
med_disc <- median(gene_pext$prop_low_pext[gene_pext$group == "Discordant35"])
med_bg   <- median(gene_pext$prop_low_pext[gene_pext$group == "Background_LOEUF02"])
n_disc_g <- sum(gene_pext$group == "Discordant35")
n_bg_g   <- sum(gene_pext$group == "Background_LOEUF02")

# Consistent subtitle font size across panels
SUBTITLE_SIZE <- 10

# (A) ECDF comparison (Discordant vs Background)
#     - legend moved inside the plot
#     - threshold label nudged right so it doesn’t overlap the dashed line
#     - red shades (two reds)
p_ecdf <- ggplot(v_all2, aes(x = exp_prop_mean, colour = group)) +
  stat_ecdf(geom = "step", linewidth = 0.9) +
  geom_vline(xintercept = PEXT_CUTOFF, linetype = "dashed", linewidth = 0.8) +
  annotate(
    "text",
    x = PEXT_CUTOFF + 0.03,  # << nudged right
    y = 0.12,
    label = paste0("Low-expression threshold (pext < ", PEXT_CUTOFF, ")"),
    hjust = 0,
    size = 3
  ) +
  scale_colour_manual(
    values = c(
      "Background_LOEUF02" = "#8B0000",  # dark red
      "Discordant35"       = "#FF0000"   # bright red
    )
  ) +
  labs(
    title = "A.",
    subtitle = paste0(
      p_variant_txt,
      " | Low-pext (pext < ", PEXT_CUTOFF, "): ",
      sprintf("Discordant %.1f%% vs Background %.1f%%", 100*frac_disc, 100*frac_bg)
    ),
    x = "Mean proportion expression (pext)",
    y = "ECDF",
    colour = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = SUBTITLE_SIZE),
    
    # ✅ Legend bottom-right, out of the ECDF lines
    legend.position = c(0.98, 0.05),
    legend.justification = c(1, 0),
    
    legend.background = element_rect(fill = "white", colour = NA),
    legend.key = element_rect(fill = "white", colour = NA)
  )

# (B) Gene-level distribution (violin + box + jitter)
#     - ensures box/median visible for BOTH groups:
#       draw boxplot AFTER violin with strong outline and opaque fill
p_gene <- ggplot(gene_pext, aes(x = group, y = prop_low_pext)) +
  geom_violin(trim = TRUE, fill = "grey90", colour = "black", linewidth = 0.7) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.35, size = 1.6) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    colour = "black",
    fill = "white",
    alpha = 0.55,
    linewidth = 0.9
  ) +
  labs(
    title = "B.",
    subtitle = paste0(
      p_gene_txt,
      " | Gene medians: ",
      sprintf("Discordant %.2f vs Background %.2f", med_disc, med_bg),
      " | n genes: ", n_disc_g, " vs ", n_bg_g
    ),
    x = NULL,
    y = paste0("Per-gene fraction of variants with pext < ", PEXT_CUTOFF)
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = SUBTITLE_SIZE),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

# Combined 2-panel
p_combined <- p_ecdf / p_gene

OUT2P_PNG <- file.path(FIG_DIR, "SuppFig_pext_expression_attenuation_discordant35_2panel.png")
OUT2P_PDF <- file.path(FIG_DIR, "SuppFig_pext_expression_attenuation_discordant35_2panel.pdf")

ggsave(OUT2P_PNG, p_combined, width = 9.0, height = 6.6, dpi = 300)
ggsave(OUT2P_PDF, p_combined, width = 9.0, height = 6.6)

# Also save the individual panels if you want
ggsave(file.path(FIG_DIR, "SuppFig_pext_ecdf_discordant_vs_bg.png"),
       p_ecdf, width = 9.0, height = 4.3, dpi = 300)

ggsave(file.path(FIG_DIR, "SuppFig_pext_genelevel_lowpext_violin.png"),
       p_gene, width = 9.0, height = 4.3, dpi = 300)

msg("Wrote figures:")
msg("  %s", OUT2P_PNG)
msg("  %s", OUT2P_PDF)
msg("  %s", file.path(FIG_DIR, "SuppFig_pext_ecdf_discordant_vs_bg.png"))
msg("  %s", file.path(FIG_DIR, "SuppFig_pext_genelevel_lowpext_violin.png"))

msg("DONE.")
