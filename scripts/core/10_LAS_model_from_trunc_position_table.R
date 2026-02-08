#!/usr/bin/env Rscript

# ============================================================
# LAS_model_from_trunc_position_table.R (clean)
# - Input: trunc_position_domain_nmd_variants.csv
#   (Gene_norm, variant_id, rel_pos, group)
# - Builds per-gene features:
#     n_trunc_sites (distinct variant sites)
#     median_relpos, IQR_relpos
# - Fits:
#     (1) Sites-only logistic model
#     (2) Position-only logistic model
#     (3) Sites+Position elastic-net (glmnet)
# - Reports stratified CV AUCs + writes outputs
# - Produces ONE 2-panel figure based on OUT-OF-FOLD predictions:
#     SuppFig_LAS_model_2panel.png/.pdf  (A=OOF ROC, B=OOF scores by class)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(rsample)
  library(pROC)
  library(glmnet)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------
# Args
# -----------------------------
option_list <- list(
  make_option(c("--variants"), type="character", default="gnomad_lof_discordance_out/trunc_position_domain_nmd_variants.csv",
              help="Variant-level CSV [default %default]."),
  make_option(c("--outdir"), type="character", default="results/las",
              help="Output directory [default %default]."),
  make_option(c("--seed"), type="integer", default=1, help="Seed [default %default]."),
  make_option(c("--folds"), type="integer", default=5, help="CV folds [default %default]."),
  make_option(c("--min_sites_per_gene"), type="integer", default=1,
              help="Min distinct trunc sites per gene [default %default].")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
set.seed(opt$seed)

# -----------------------------
# Read + validate
# -----------------------------
message("[LAS] Reading: ", opt$variants)

if (!file.exists(opt$variants)) {
  stop("Input variants file not found: ", opt$variants)
}

v <- readr::read_csv(opt$variants, show_col_types = FALSE)

need_cols <- c("Gene_norm","variant_id","rel_pos","group")
miss <- setdiff(need_cols, names(v))
if (length(miss)) stop("Missing required columns in variants file: ", paste(miss, collapse=", "))

is_discordant_from_group <- function(g) {
  gl <- tolower(trimws(as.character(g)))
  ifelse(grepl("^disc", gl), 1L, 0L)
}

dat <- v %>%
  transmute(
    gene = toupper(as.character(Gene_norm)),
    variant_id = as.character(variant_id),
    rel_pos = suppressWarnings(as.numeric(rel_pos)),
    group = as.character(group)
  ) %>%
  mutate(is_discordant = is_discordant_from_group(group)) %>%
  filter(!is.na(gene), gene != "", !is.na(variant_id), variant_id != "")

message("[LAS] Variants: ", nrow(dat))
message("[LAS] Distinct genes: ", n_distinct(dat$gene))
message("[LAS] Discordant genes present: ", n_distinct(dat$gene[dat$is_discordant == 1L]))

# -----------------------------
# Per-gene features
# -----------------------------
gene_features <- dat %>%
  group_by(gene) %>%
  summarise(
    is_discordant = max(is_discordant),
    n_trunc_sites = n_distinct(variant_id),
    median_relpos = median(rel_pos, na.rm = TRUE),
    iqr_relpos = IQR(rel_pos, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_trunc_sites >= opt$min_sites_per_gene) %>%
  mutate(log_trunc_sites = log1p(n_trunc_sites))

# Impute missing relpos summaries if any (rare, but safe)
for (cc in c("median_relpos","iqr_relpos")) {
  m <- median(gene_features[[cc]], na.rm = TRUE)
  if (!is.finite(m)) m <- 0
  gene_features[[cc]] <- ifelse(is.na(gene_features[[cc]]), m, gene_features[[cc]])
}

message("[LAS] Genes for modelling: ", nrow(gene_features),
        " (discordant=", sum(gene_features$is_discordant == 1L), ")")

if (sum(gene_features$is_discordant == 1L) < 5) {
  warning("Very few discordant genes after filtering; CV AUC will be unstable.")
}

# -----------------------------
# Stratified CV AUC helper (GLM)
# -----------------------------
cv_auc_glm <- function(df, formula, v=5, seed=1) {
  set.seed(seed)
  folds <- vfold_cv(df, v = v, strata = is_discordant)
  
  aucs <- map_dbl(folds$splits, function(sp) {
    tr <- analysis(sp); te <- assessment(sp)
    m <- glm(formula, data = tr, family = binomial())
    p <- as.numeric(predict(m, newdata = te, type = "response"))
    roc_obj <- pROC::roc(te$is_discordant, p, quiet = TRUE)
    as.numeric(pROC::auc(roc_obj))
  })
  
  list(mean = mean(aucs), sd = sd(aucs), per_fold = aucs)
}

# -----------------------------
# 1) Sites-only model
# -----------------------------
res_sites <- cv_auc_glm(
  gene_features, is_discordant ~ log_trunc_sites,
  v = opt$folds, seed = opt$seed
)

# -----------------------------
# 2) Position-only model
# -----------------------------
res_pos <- cv_auc_glm(
  gene_features, is_discordant ~ median_relpos + iqr_relpos,
  v = opt$folds, seed = opt$seed
)

# -----------------------------
# 3) Elastic-net (sites + position)
# Outer CV for AUC + out-of-fold predictions for honest ROC/score plots
# -----------------------------
set.seed(opt$seed)
folds <- vfold_cv(gene_features, v = opt$folds, strata = is_discordant)

oof_list <- vector("list", length(folds$splits))
aucs_enet <- numeric(length(folds$splits))

for (i in seq_along(folds$splits)) {
  sp <- folds$splits[[i]]
  tr <- analysis(sp); te <- assessment(sp)
  
  x_tr <- as.matrix(tr %>% select(log_trunc_sites, median_relpos, iqr_relpos))
  y_tr <- tr$is_discordant
  
  x_te <- as.matrix(te %>% select(log_trunc_sites, median_relpos, iqr_relpos))
  y_te <- te$is_discordant
  
  fit_cv <- cv.glmnet(
    x = x_tr, y = y_tr,
    family = "binomial",
    alpha = 0.5,
    nfolds = 5,
    type.measure = "auc"
  )
  
  p <- as.numeric(predict(fit_cv, newx = x_te, s = "lambda.min", type = "response"))
  roc_obj <- pROC::roc(y_te, p, quiet = TRUE)
  aucs_enet[i] <- as.numeric(pROC::auc(roc_obj))
  
  oof_list[[i]] <- tibble(
    gene = te$gene,
    y = y_te,
    p = p
  )
}

oof_df <- bind_rows(oof_list)
res_enet <- list(mean = mean(aucs_enet), sd = sd(aucs_enet), per_fold = aucs_enet)

# -----------------------------
# Final elastic-net fit on full data (DESCRIPTIVE ONLY)
# Used to provide coefficients + a per-gene score column in output tables.
# -----------------------------
x_all <- as.matrix(gene_features %>% select(log_trunc_sites, median_relpos, iqr_relpos))
y_all <- gene_features$is_discordant

final_cv <- cv.glmnet(
  x_all, y_all,
  family = "binomial",
  alpha = 0.5,
  nfolds = 5,
  type.measure = "auc"
)

gene_features <- gene_features %>%
  mutate(LAS_fullfit = as.numeric(predict(final_cv, newx = x_all, s = "lambda.min", type = "response")))

coef_mat <- as.matrix(coef(final_cv, s = "lambda.min"))
coef_df <- tibble(
  term = rownames(coef_mat),
  coef = as.numeric(coef_mat[, 1]),
  odds_ratio = exp(as.numeric(coef_mat[, 1]))
)

# -----------------------------
# Honest 2-panel figure (OUT-OF-FOLD)
# A: OOF ROC, B: OOF predicted probabilities by class
# -----------------------------
roc_oof <- pROC::roc(oof_df$y, oof_df$p, quiet = TRUE)
auc_oof <- as.numeric(pROC::auc(roc_oof))

roc_oof_df <- tibble(
  fpr = 1 - roc_oof$specificities,
  tpr = roc_oof$sensitivities
)

pA <- ggplot(roc_oof_df, aes(x = fpr, y = tpr)) +
  geom_line(linewidth = 0.9) +
  geom_abline(linetype = 2) +
  coord_equal() +
  theme_classic(base_size = 12) +
  labs(
    title = "ROC curve (out-of-fold)",
    subtitle = sprintf("%d-fold CV AUC = %.3f ± %.3f", opt$folds, res_enet$mean, res_enet$sd),
    x = "False positive rate",
    y = "True positive rate"
  )

oof_plot_df <- oof_df %>%
  mutate(class = factor(y, levels = c(0, 1), labels = c("Non-discordant", "Discordant")))

pB <- ggplot(oof_plot_df, aes(x = class, y = p)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA) +
  theme_classic(base_size = 12) +
  labs(
    title = "Predicted probability by class",
    x = NULL,
    y = "Predicted probability"
  )

fig <- (pA | pB) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 16, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

FIG_DIR <- "figures"
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)
out_png <- file.path(FIG_DIR, "SuppFig_LAS_model_2panel.png")
out_pdf <- file.path(FIG_DIR, "SuppFig_LAS_model_2panel.pdf")
ggsave(out_png, fig, width = 10.5, height = 4.8, units = "in", dpi = 600)
ggsave(out_pdf, fig, width = 10.5, height = 4.8, units = "in")
message("[LAS] Wrote: ", out_png)
message("[LAS] Wrote: ", out_pdf)

# -----------------------------
# Write outputs
# -----------------------------
readr::write_tsv(gene_features, file.path(opt$outdir, "LAS_gene_features.tsv"))
readr::write_tsv(coef_df,       file.path(opt$outdir, "LAS_elasticnet_coefficients.tsv"))
readr::write_tsv(oof_df,        file.path(opt$outdir, "LAS_enet_oof_predictions.tsv"))

summary_path <- file.path(opt$outdir, "LAS_model_summary.txt")
writeLines(c(
  "LAS modelling summary",
  "====================",
  sprintf("Input variants: %s", opt$variants),
  sprintf("Genes modelled: %d (discordant=%d)", nrow(gene_features), sum(gene_features$is_discordant == 1L)),
  "",
  sprintf("Sites-only GLM:      %d-fold CV AUC = %.3f ± %.3f", opt$folds, res_sites$mean, res_sites$sd),
  sprintf("Position-only GLM:   %d-fold CV AUC = %.3f ± %.3f", opt$folds, res_pos$mean, res_pos$sd),
  sprintf("Elastic-net (all):   %d-fold CV AUC = %.3f ± %.3f", opt$folds, res_enet$mean, res_enet$sd),
  "",
  "Per-fold AUCs:",
  paste0("  Sites-only:    ", paste(sprintf("%.3f", res_sites$per_fold), collapse = ", ")),
  paste0("  Position-only: ", paste(sprintf("%.3f", res_pos$per_fold), collapse = ", ")),
  paste0("  Elastic-net:   ", paste(sprintf("%.3f", res_enet$per_fold), collapse = ", ")),
  "",
  "Final elastic-net coefficients at lambda.min (full-data fit; descriptive only):",
  paste0("  ", coef_df$term, "\t", sprintf("%.6f", coef_df$coef), collapse = "\n")
), con = summary_path)

message(sprintf("[LAS] Sites-only CV AUC = %.3f ± %.3f", res_sites$mean, res_sites$sd))
message(sprintf("[LAS] Position-only CV AUC = %.3f ± %.3f", res_pos$mean, res_pos$sd))
message(sprintf("[LAS] Elastic-net CV AUC = %.3f ± %.3f", res_enet$mean, res_enet$sd))
message("[LAS] Wrote outputs to: ", opt$outdir)
message("[LAS] Done.")
