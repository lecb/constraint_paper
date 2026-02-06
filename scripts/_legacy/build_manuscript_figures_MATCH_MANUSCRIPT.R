#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)   # <-- ADD THIS
  library(ggbeeswarm)
})

# ============================================================
# Cache utilities
# ============================================================

CACHE_DIR <- "cache"

if (!dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR, recursive = TRUE)
}

cache_path <- function(name) {
  file.path(CACHE_DIR, paste0(name, ".rds"))
}

cache_exists <- function(name) {
  file.exists(cache_path(name))
}

read_cache <- function(name) {
  readRDS(cache_path(name))
}

write_cache <- function(name, obj) {
  saveRDS(obj, cache_path(name))
}


# ============================================================
# build_manuscript_figures_MATCH_MANUSCRIPT.R
# 2026-01-20
#
# Goal: regenerate figures to MATCH manuscriptv2.docx (embedded images).
#
# Outputs:
#   manuscript_figures/main/
#   manuscript_figures/supp/
#
# Run:
#   source("build_manuscript_figures_MATCH_MANUSCRIPT.R", local=TRUE)
#   run_fig_build(mode="all", force=TRUE)
#   run_fig_build(mode="one", fig="Figure_LOEUF_vs_trunc_AC_discordant_highlight", force=TRUE)
# ============================================================

# ----------------------------
# Config
# ----------------------------
CFG <- list(
  CONSTRAINT_TSV   = "gnomad.v4.1.constraint_metrics.tsv",
  DISCORDANT_RDS   = "discordant_genes.rds",          # expects a data frame with a gene symbol column
  TABLE1_CSV       = "Table1_final_main.csv",         # expects column 'Gene'
  VARIANTS_CSV     = "gnomad_variants_all.csv",       # optional; used for trunc_AC aggregation
  TRANSCRIPT_TSV   = "constraint_transcripts.tsv",    # optional; for canonical vs min LOEUF
  OUT_MAIN_DIR     = file.path("manuscript_figures","main"),
  OUT_SUPP_DIR     = file.path("manuscript_figures","supp"),
  CACHE_DIR        = "cache",
  variants_csv = "gnomad_variants_all.csv",
  LOEUF_CUTOFF     = 0.2,
  SEED            = 1,
  base_size = 12
)

# ---- HARD GUARD: ensure output paths exist and are valid strings ----
if (is.null(CFG$out_root) || !nzchar(CFG$out_root)) {
  CFG$out_root <- "manuscript_figures"
}
if (is.null(CFG$out_main) || !nzchar(CFG$out_main)) {
  CFG$out_main <- file.path(CFG$out_root, "main")
}
if (is.null(CFG$out_supp) || !nzchar(CFG$out_supp)) {
  CFG$out_supp <- file.path(CFG$out_root, "supp")
}

stopifnot(is.character(CFG$out_root), length(CFG$out_root) == 1, nzchar(CFG$out_root))
stopifnot(is.character(CFG$out_main), length(CFG$out_main) == 1, nzchar(CFG$out_main))
stopifnot(is.character(CFG$out_supp), length(CFG$out_supp) == 1, nzchar(CFG$out_supp))

dir.create(CFG$out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_main, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_supp, showWarnings = FALSE, recursive = TRUE)


# --- sanity: ensure scalar cutoff ---
if (is.null(CFG$loeuf_cutoff)) CFG$loeuf_cutoff <- 0.20
CFG$loeuf_cutoff <- as.numeric(CFG$loeuf_cutoff)[1]
stopifnot(is.finite(CFG$loeuf_cutoff))

dir.create(CFG$out_supp, showWarnings = FALSE, recursive = TRUE)



dir.create(CFG$OUT_MAIN_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(CFG$OUT_SUPP_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(CFG$CACHE_DIR, recursive=TRUE, showWarnings=FALSE)

# ----------------------------
# Figure output paths
# ----------------------------
FIG_OUT <- list(
  main = list(
    Figure_LOEUF_vs_truncating_clean_Table1Top15label =
      file.path(CFG$out_main, "Figure_LOEUF_vs_truncating_clean_Table1Top15label.png"),
    Figure_LOEUF_vs_trunc_AC_discordant_highlight =
      file.path(CFG$out_main, "Figure_LOEUF_vs_trunc_AC_discordant_highlight.png"),
    Figure1_inset_trunc_per_gene =
      file.path(CFG$out_main, "Figure1_inset_trunc_per_gene.png"),
    SuppFig_LOEUF_vs_observed_truncations_top15 =
      file.path(CFG$out_supp, "SuppFig_LOEUF_vs_observed_truncations_top15.png"),
    Figure2_discordant_vs_rest =
      file.path(CFG$out_main, "Figure2_discordant_vs_rest.png")
  ),
  supp = list(
    Supplementary_Figure_NMD_escape_per_gene =
      file.path(CFG$out_supp, "Supplementary_Figure_NMD_escape_per_gene.png"),
    SuppFig_length_adjusted_2panel =
      file.path(CFG$out_supp, "SuppFig_length_adjusted_2panel.png"),
    SuppFig_perm_null_2panel =
      file.path(CFG$out_supp, "SuppFig_perm_null_2panel.png"),
    SuppFig_negative_control_syn_mis =
      file.path(CFG$out_supp, "SuppFig_negative_control_syn_mis.png"),
    SuppFig_canonical_vs_min_LOEUF =
      file.path(CFG$out_supp, "SuppFig_canonical_vs_min_LOEUF.png"),
    SuppFig_trunc_burden_adjusted_by_syn_obs =
      file.path(CFG$out_supp, "SuppFig_trunc_burden_adjusted_by_syn_obs.png")
  )
)


# ----------------------------
# Helpers
# ----------------------------
norm_sym <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\s+", "")
  x <- toupper(x)
  x
}

add_gene_labels_like_length_panel <- function(p, df_labels,
                                              label_col = "Gene",
                                              size = 4.5,
                                              seed = 1) {
  
  p +
    ggrepel::geom_label_repel(
      data = df_labels,
      aes(label = .data[[label_col]]),
      
      size = size,
      fill = "white",
      color = "black",
      
      # MATCH SuppFig_length_adjusted_2panel
      label.size    = 0.25,
      box.padding   = 0.6,
      point.padding = 0.35,
      
      min.segment.length = 0,
      segment.color = "grey40",
      
      max.overlaps = Inf,
      seed = seed
    )
}


move_subtitle_to_bottom <- function(p, size = 3.0, pad_frac = 0.03) {
  sub <- p$labels$subtitle
  if (is.null(sub) || !nzchar(sub)) return(p)
  
  gb <- ggplot2::ggplot_build(p)
  pp <- gb$layout$panel_params[[1]]
  
  # panel y-range (often in *transformed* space)
  yr <- pp$y.range
  y0 <- yr[1]
  y1 <- yr[2]
  
  # get the y scale transform (newer ggplot2 stores it here)
  ytrans <- gb$layout$panel_scales_y[[1]]$trans
  has_inverse <- !is.null(ytrans) && !is.null(ytrans$inverse)
  
  # decide x position (left edge)
  is_x_discrete <- !is.null(pp$x$breaks)
  x_pos <- if (is_x_discrete) 1 else -Inf
  
  # compute y position safely *in panel space*
  if (has_inverse) {
    # convert panel-space bounds -> data-space, nudge in data-space, then re-transform
    y0_data <- ytrans$inverse(y0)
    y1_data <- ytrans$inverse(y1)
    
    # nudge upward from the bottom in data space
    y_pos_data <- y0_data + pad_frac * (y1_data - y0_data)
    
    # transform back to panel space for annotate()
    y_pos <- ytrans$transform(y_pos_data)
  } else {
    # no transform: just nudge in the same space
    y_pos <- y0 + pad_frac * (y1 - y0)
  }
  
  p +
    ggplot2::labs(subtitle = NULL) +
    ggplot2::annotate(
      "text",
      x = x_pos, y = y_pos,
      label = sub,
      hjust = if (is_x_discrete) 0 else -0.02,
      vjust = 0,
      size = size
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(t = 8, r = 10, b = 30, l = 10))
}


msg <- function(...) cat(sprintf(...), "\n")

read_rds_safe <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  readRDS(path)
}

detect_gene_col <- function(df) {
  candidates <- c("Gene","gene","gene_symbol","symbol","hgnc_symbol","hgnc","GENE")
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) stop("Could not detect gene column in: ", paste(names(df), collapse=", "))
  hit[1]
}

move_stats_to_bottom <- function(p, size = 9) {
  p +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(
        size = size,
        hjust = 0,
        vjust = 0,
        margin = ggplot2::margin(t = 0, b = 0)
      )
    ) +
    ggplot2::labs(
      subtitle = paste0("\n\n", p$labels$subtitle)
    )
}


# ----------------------------
# Robust AC detection / construction
# ----------------------------
detect_ac_cols <- function(v) {
  nms <- names(v)
  
  # helper: first hit by exact match (case-insensitive)
  pick <- function(cands) {
    hit <- cands[tolower(cands) %in% tolower(nms)][1]
    if (is.na(hit) || is.null(hit)) return(NULL)
    # return the actual column name as present in v
    nms[tolower(nms) == tolower(hit)][1]
  }
  
  exome <- pick(c("exome_ac", "exome_AC", "exomeAc", "exome.ac"))
  genome <- pick(c("genome_ac", "genome_AC", "genomeAc", "genome.ac"))
  
  # common single-column fallbacks
  single <- pick(c(
    "ac", "AC", "allele_count", "allelecount", "alleleCount",
    "total_ac", "Total_AC", "AN_AC", "obs_ac"
  ))
  
  list(exome = exome, genome = genome, single = single)
}

compute_total_ac <- function(v) {
  cols <- detect_ac_cols(v)
  
  # if exome/genome exist, sum whichever exist
  if (!is.null(cols$exome) || !is.null(cols$genome)) {
    ex <- if (!is.null(cols$exome)) suppressWarnings(as.numeric(v[[cols$exome]])) else 0
    gn <- if (!is.null(cols$genome)) suppressWarnings(as.numeric(v[[cols$genome]])) else 0
    tot <- ex + gn
    tot[is.na(tot)] <- 0
    return(list(total_ac = tot, ac_source = paste(na.omit(c(cols$exome, cols$genome)), collapse = " + ")))
  }
  
  # else try single AC column
  if (!is.null(cols$single)) {
    tot <- suppressWarnings(as.numeric(v[[cols$single]]))
    tot[is.na(tot)] <- 0
    return(list(total_ac = tot, ac_source = cols$single))
  }
  
  stop(
    "Variants CSV: could not detect allele count columns.\n",
    "Looked for exome_ac/genome_ac first, then common single columns (AC/ac/allele_count/etc).\n",
    "Available columns include:\n  ",
    paste(head(names(v), 80), collapse = ", ")
  )
}


detect_consequence_col <- function(df) {
  candidates <- c("Consequence","consequence","csq","most_severe_consequence","VEP_consequence")
  hit <- intersect(candidates, names(df))
  if (length(hit)==0) stop("Variants CSV: could not detect consequence column.")
  hit[1]
}

is_trunc_csq <- function(csq) {
  csq <- tolower(as.character(csq))
  trunc_terms <- c("stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant",
                   "start_lost","transcript_ablation","exon_loss_variant")
  str_detect(csq, paste(trunc_terms, collapse="|"))
}

theme_manuscript <- function(base_size = 18) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face="bold", size=base_size*1.35, hjust=0),
      plot.subtitle = element_text(size=base_size*0.9, hjust=0),
      axis.title = element_text(face="bold"),
      axis.text = element_text(color="black"),
      legend.position = "none"
    )
}

# ============================================================
# Shrink aesthetics ONLY for super-figures
# ============================================================
shrink_for_superfig <- function(p,
                                point_scale = 0.65,
                                text_scale  = 0.85,
                                y_title_scale = 0.75) {
  
  # shrink points + ggrepel labels only for the super-figure
  for (i in seq_along(p$layers)) {
    lyr <- p$layers[[i]]
    
    if (inherits(lyr$geom, c("GeomPoint", "GeomJitter"))) {
      if (!is.null(lyr$aes_params$size)) {
        lyr$aes_params$size <- lyr$aes_params$size * point_scale
      }
    }
    
    if (inherits(lyr$geom, c("GeomTextRepel", "GeomLabelRepel"))) {
      if (!is.null(lyr$aes_params$size)) {
        lyr$aes_params$size <- lyr$aes_params$size * text_scale
      }
    }
    
    p$layers[[i]] <- lyr
  }
  
  # shrink axis text + titles; specifically make Y title smaller + not bold
  p + ggplot2::theme(
    axis.text     = ggplot2::element_text(size = ggplot2::rel(text_scale), face = "plain"),
    axis.title.x  = ggplot2::element_text(size = ggplot2::rel(text_scale), face = "plain"),
    axis.title.y  = ggplot2::element_text(
      size = ggplot2::rel(text_scale * y_title_scale),
      face = "plain"
    )
  )
}




# ----------------------------
# Cache builders
# ----------------------------
cache_path <- function(key) file.path(CFG$CACHE_DIR, paste0(key, ".rds"))

cache_get_or_build <- function(key, force=FALSE, build_fn) {
  p <- cache_path(key)
  if (!force && file.exists(p)) return(readRDS(p))
  msg("[CACHE BUILD] %s", key)
  obj <- build_fn()
  saveRDS(obj, p)
  obj
}

build_constraint_gene <- function() {
  stopifnot(file.exists(CFG$CONSTRAINT_TSV))
  cg <- readr::read_tsv(CFG$CONSTRAINT_TSV, show_col_types = FALSE)  # <-- cg, not df_raw
  
  nms <- names(cg)
  
  # gene col
  gcol <- if ("gene" %in% nms) "gene" else detect_gene_col(cg)
  
  # LOEUF upper CI col (gnomAD varies)
  loeuf_col <- if ("lof.oe_ci.upper" %in% nms) {
    "lof.oe_ci.upper"
  } else if ("lof_oe_ci_upper" %in% nms) {
    "lof_oe_ci_upper"
  } else if ("lof_oe_ci_upper_bound" %in% nms) {
    "lof_oe_ci_upper_bound"
  } else {
    stop("Constraint TSV: cannot find LOEUF upper CI column. Have: ", paste(nms, collapse = ", "))
  }
  
  # LoF observed col
  lofobs_col <- if ("lof_hc_lc.obs" %in% nms) {
    "lof_hc_lc.obs"
  } else if ("lof.obs" %in% nms) {
    "lof.obs"
  } else if ("lof_obs" %in% nms) {
    "lof_obs"
  } else {
    stop("Constraint TSV: cannot find LoF observed column. Have: ", paste(nms, collapse = ", "))
  }
  
  synobs_col <- if ("syn.obs" %in% nms) "syn.obs" else NA_character_
  misoe_col  <- if ("mis.oe" %in% nms) "mis.oe" else NA_character_
  
  cg2 <- cg %>%
    dplyr::transmute(
      Gene    = norm_sym(.data[[gcol]]),
      LOEUF   = suppressWarnings(as.numeric(.data[[loeuf_col]])),
      lof_obs = suppressWarnings(as.numeric(.data[[lofobs_col]])),
      syn_obs = if (!is.na(synobs_col)) suppressWarnings(as.numeric(.data[[synobs_col]])) else NA_real_,
      mis_oe  = if (!is.na(misoe_col))  suppressWarnings(as.numeric(.data[[misoe_col]]))  else NA_real_
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", is.finite(LOEUF))
  
  if (nrow(cg2) == 0) {
    stop("Constraint TSV read OK, but produced 0 usable rows after parsing LOEUF/Gene. Check column formats.")
  }
  
  cg2
}

build_discordant_genes <- function() {
  d <- read_rds_safe(CFG$DISCORDANT_RDS)
  if (is.vector(d)) d <- tibble(Gene = d)
  gcol <- detect_gene_col(d)
  d %>% transmute(Gene = norm_sym(.data[[gcol]])) %>% distinct()
}

build_table1 <- function() {
  stopifnot(file.exists(CFG$TABLE1_CSV))
  t1 <- readr::read_csv(CFG$TABLE1_CSV, show_col_types = FALSE)
  stopifnot("Gene" %in% names(t1))
  
  t1 %>%
    dplyr::transmute(Gene = norm_sym(Gene)) %>%
    dplyr::filter(!is.na(Gene), Gene != "") %>%
    dplyr::distinct() %>%
    dplyr::slice_head(n = 15)
}


build_gene_trunc_ac <- function(force = FALSE) {
  
  if (cache_exists("gene_trunc_ac") && !force) {
    return(read_cache("gene_trunc_ac"))
  }
  
  stopifnot(file.exists(CFG$variants_csv))
  v <- readr::read_csv(CFG$variants_csv, show_col_types = FALSE)
  
  # --- detect required columns robustly ---
  gcol <- detect_gene_col(v)  # tries Gene/gene/gene_symbol/symbol/etc
  ccol <- detect_consequence_col(v)
  
  # variant id column (robust-ish)
  vid_candidates <- c("variant_id","variantId","id","var_id","VID")
  vid <- intersect(vid_candidates, names(v))[1]
  if (is.na(vid) || is.null(vid)) {
    # fallback: if no ID column, make one from row number
    v$.__rowid <- seq_len(nrow(v))
    vid <- ".__rowid"
  }
  
  # compute total AC using your robust AC detector
  ac_info <- compute_total_ac(v)
  v$total_ac <- ac_info$total_ac
  message("[INFO] variants AC source: ", ac_info$ac_source)
  
  trunc_tbl <- v %>%
    dplyr::mutate(
      Gene = norm_sym(.data[[gcol]]),
      consequence_norm = tolower(as.character(.data[[ccol]]))
    ) %>%
    dplyr::filter(grepl("stop_gained|frameshift_variant", consequence_norm)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      trunc_n_variants = dplyr::n_distinct(.data[[vid]]),
      trunc_AC = sum(total_ac, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_cache("gene_trunc_ac", trunc_tbl)
  trunc_tbl
}



build_constraint_tx <- function() {
  
  # Preferred (matches manuscript): use your precomputed table if present
  if (file.exists("canonical_vs_min_LOEUF_table.csv")) {
    x <- readr::read_csv("canonical_vs_min_LOEUF_table.csv", show_col_types = FALSE)
    
    # Try to detect columns robustly
    nms <- names(x)
    
    pick1 <- function(cands) {
      hit <- cands[cands %in% nms]
      if (length(hit) == 0) return(NULL)
      hit[1]
    }
    
    gene_col <- pick1(c("Gene", "gene", "gene_symbol"))
    min_col  <- pick1(c("min_LOEUF", "min_loeuf", "min_lof_oe_ci_upper", "min_lof_oe_ci.upper"))
    can_col  <- pick1(c("canonical_LOEUF", "canonical_loeuf", "canonical_lof_oe_ci_upper", "canonical_lof_oe_ci.upper"))
    
    if (is.null(gene_col) || is.null(min_col) || is.null(can_col)) {
      stop(
        "canonical_vs_min_LOEUF_table.csv is missing expected columns.\n",
        "Have: ", paste(nms, collapse = ", "), "\n",
        "Need gene + min_LOEUF + canonical_LOEUF (or close equivalents)."
      )
    }
    
    out <- x %>%
      dplyr::transmute(
        Gene = norm_sym(.data[[gene_col]]),
        min_LOEUF = suppressWarnings(as.numeric(.data[[min_col]])),
        canonical_LOEUF = suppressWarnings(as.numeric(.data[[can_col]]))
      ) %>%
      dplyr::filter(is.finite(min_LOEUF), is.finite(canonical_LOEUF))
    
    if (nrow(out) == 0) stop("canonical_vs_min_LOEUF_table.csv parsed but produced 0 finite rows.")
    
    return(out)
  }
  
  # Fallback: attempt to derive from constraint TSV IF it has transcript + canonical fields
  stopifnot(file.exists(CFG$constraint_tsv))
  df_raw <- readr::read_tsv(CFG$constraint_tsv, show_col_types = FALSE)
  nms <- names(df_raw)
  
  if (!("gene" %in% nms) || !("lof.oe_ci.upper" %in% nms)) {
    stop("Constraint TSV missing gene or lof.oe_ci.upper. Cannot build canonical-vs-min plot.")
  }
  
  tx_id_col <- dplyr::case_when(
    "transcript" %in% nms ~ "transcript",
    "transcript_id" %in% nms ~ "transcript_id",
    "transcriptId" %in% nms ~ "transcriptId",
    TRUE ~ NA_character_
  )
  
  if (is.na(tx_id_col) || !("canonical" %in% nms)) {
    stop(
      "No canonical_vs_min_LOEUF_table.csv found, and constraint TSV does not appear transcript-level ",
      "(missing transcript id and/or canonical column)."
    )
  }
  
  tx_df <- df_raw %>%
    dplyr::transmute(
      Gene = norm_sym(gene),
      transcript_id = as.character(.data[[tx_id_col]]),
      canonical = as.logical(canonical),
      LOEUF = suppressWarnings(as.numeric(`lof.oe_ci.upper`))
    ) %>%
    dplyr::filter(is.finite(LOEUF))
  
  # summarise to the table required for the plot
  out <- tx_df %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      min_LOEUF = min(LOEUF, na.rm = TRUE),
      canonical_LOEUF = suppressWarnings(min(LOEUF[canonical %in% TRUE], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(min_LOEUF), is.finite(canonical_LOEUF))
  
  if (nrow(out) == 0) stop("Derived transcript summary but got 0 rows with finite min/canonical LOEUF.")
  
  out
}

# ----------------------------
# helper theme (must exist before any figure functions run)
# ----------------------------
theme_paper <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text  = ggplot2::element_text(size = base_size * 0.85)
    )
}

add_panel_label <- function(p, label,
                            x = -Inf, y = Inf,
                            hjust = -0.3, vjust = 1.4,
                            size = 6) {
  p + ggplot2::annotate(
    "text",
    x = x, y = y,
    label = label,
    hjust = hjust, vjust = vjust,
    size = size,
    fontface = "bold"
  )
}

# ============================================================
# Super-figure axis normaliser (readable, consistent, not bold)
# ============================================================
apply_superfig_axes <- function(p,
                                axis_title_size = 11,
                                axis_text_size  = 10,
                                subtitle_size   = 10) {
  p + ggplot2::theme(
    axis.title.x  = ggplot2::element_text(size = axis_title_size, face = "plain"),
    axis.title.y  = ggplot2::element_text(size = axis_title_size, face = "plain"),
    axis.text.x   = ggplot2::element_text(size = axis_text_size, face = "plain", color = "black"),
    axis.text.y   = ggplot2::element_text(size = axis_text_size, face = "plain", color = "black"),
    plot.subtitle = ggplot2::element_text(size = subtitle_size, face = "plain")
  )
}

# ----------------------------
# Figure builders (MATCH manuscript)
# ----------------------------

# (1) Figure_LOEUF_vs_truncating_clean_Table1Top15label.png  [MAIN]
make_Figure_LOEUF_vs_truncating_clean_Table1Top15label <- function(constraint_gene, table1) {
  
  # enforce 1 row per gene (protect against transcript-level constraint input)
  df <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF   = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      lof_obs = suppressWarnings(max(lof_obs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(LOEUF), is.finite(lof_obs))
  
  df <- df %>%
    dplyr::mutate(
      y = lof_obs + 1,  # +1 pseudocount
      is_table1 = Gene %in% table1$Gene
    ) %>%
    dplyr::filter(is.finite(y), y > 0)
  
  label_genes <- intersect(table1$Gene, df$Gene)
  
  # ---- LOESS in log10 space (matches manuscript look on log y-axis) ----
  df_fit <- df %>% dplyr::filter(is.finite(LOEUF), is.finite(y), y > 0)
  
  lo <- stats::loess(log10(y) ~ LOEUF, data = df_fit, span = 0.75, degree = 1)
  
  xgrid <- data.frame(LOEUF = seq(0, 0.32, length.out = 300))
  xgrid$y_hat <- 10^(stats::predict(lo, newdata = xgrid))
  xgrid <- xgrid %>% dplyr::filter(is.finite(y_hat), y_hat > 0)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = LOEUF, y = y)) +
    ggplot2::geom_point(color = "grey70", alpha = 0.25, size = 2) +
    ggplot2::geom_line(
      data = xgrid,
      ggplot2::aes(x = LOEUF, y = y_hat),
      color = "grey40",
      linewidth = 1
    ) +
    ggplot2::geom_vline(
      xintercept = CFG$LOEUF_CUTOFF,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = 4,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(is_table1),
      ggplot2::aes(x = LOEUF, y = y),
      fill = "red",
      color = "black",
      shape = 21,
      size = 4,
      stroke = 0.6
    ) +
    ggrepel::geom_text_repel(
      data = df %>% dplyr::filter(Gene %in% label_genes),
      ggplot2::aes(label = Gene),
      size = 5,
      box.padding = 0.35,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = "grey40",
      max.overlaps = Inf
    ) +
    ggplot2::scale_y_log10(
      breaks = c(1, 10, 100, 1000),
      labels = c("1", "10", "100", "1000")
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 0.32), ylim = c(1, 1200)) +
    ggplot2::labs(
      title = "Highly constrained genes can still tolerate structured partial truncation",
      x = "LOEUF (lower = stronger constraint vs true LoF)",
      y = "Observed truncating variants per gene (log10; +1 pseudocount)"
    ) +
    theme_manuscript(base_size = 18)
  
  p
}


# (2) Figure_LOEUF_vs_trunc_AC_discordant_highlight.png  [MAIN]
make_Figure_LOEUF_vs_trunc_AC_discordant_highlight <- function(constraint_gene,
                                                               discordant_35,
                                                               gene_trunc_ac,
                                                               table1_top15 = NULL) {
  
  # ---- FORCE a standard Gene column no matter what came in ----
  cg <- constraint_gene
  if (!("Gene" %in% names(cg))) {
    gcol <- detect_gene_col(cg)  # will find gene/gene_symbol/symbol/etc
    cg <- cg %>% dplyr::mutate(Gene = norm_sym(.data[[gcol]]))
  } else {
    cg <- cg %>% dplyr::mutate(Gene = norm_sym(Gene))
  }
  
  dg <- discordant_35
  if (!("Gene" %in% names(dg))) {
    gcol <- detect_gene_col(dg)
    dg <- dg %>% dplyr::mutate(Gene = norm_sym(.data[[gcol]]))
  } else {
    dg <- dg %>% dplyr::mutate(Gene = norm_sym(Gene))
  }
  
  t1 <- table1_top15
  if (!is.null(t1)) {
    if (!("Gene" %in% names(t1))) {
      gcol <- detect_gene_col(t1)
      t1 <- t1 %>% dplyr::mutate(Gene = norm_sym(.data[[gcol]]))
    } else {
      t1 <- t1 %>% dplyr::mutate(Gene = norm_sym(Gene))
    }
  }
  
  gt <- gene_trunc_ac
  if (!("Gene" %in% names(gt))) {
    gcol <- detect_gene_col(gt)
    gt <- gt %>% dplyr::mutate(Gene = norm_sym(.data[[gcol]]))
  } else {
    gt <- gt %>% dplyr::mutate(Gene = norm_sym(Gene))
  }
  
  # ---- now proceed using cg/dg/gt/t1 (NOT the raw inputs) ----
  df <- cg %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", is.finite(LOEUF))
  
  gt2 <- gt %>%
    dplyr::transmute(
      Gene = Gene,
      trunc_AC = suppressWarnings(as.numeric(trunc_AC))
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "")
  
  df <- df %>%
    dplyr::left_join(gt2, by = "Gene") %>%
    dplyr::mutate(
      trunc_AC = dplyr::coalesce(trunc_AC, 0),
      y = trunc_AC + 1
    ) %>%
    dplyr::filter(is.finite(y), y > 0)
  
  disc_genes <- norm_sym(dg$Gene)
  
  top15_genes <- if (!is.null(t1) && "Gene" %in% names(t1)) {
    norm_sym(t1$Gene)
  } else {
    head(disc_genes, 15)
  }
  
  df <- df %>%
    dplyr::mutate(
      is_disc  = Gene %in% disc_genes,
      is_top15 = Gene %in% top15_genes
    )
  
  label_genes <- intersect(top15_genes, df$Gene)
  df_labels <- df %>% dplyr::filter(Gene %in% label_genes)
  
  ylab <- "Truncating allele count per gene (+1 pseudocount; log10)"
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = LOEUF, y = y)) +
    ggplot2::geom_point(color = "grey65", alpha = 0.5, size = 2.8) +
    ggplot2::geom_point(data = df %>% dplyr::filter(is_disc), color = "grey55", size = 4) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(is_top15),
      shape = 21, fill = "red", color = "black", size = 4.8, stroke = 0.7
    ) +
    ggplot2::scale_y_log10() +
    ggplot2::coord_cartesian(xlim = c(0, 0.22)) +
    ggplot2::labs(x = "LOEUF", y = ylab) +
    theme_manuscript(base_size = 18)
  
  p <- add_gene_labels_like_length_panel(p, df_labels, label_col = "Gene", size = 4.5, seed = CFG$SEED)
  p
}






# (3) Figure1_inset_trunc_per_gene.png  [MAIN]  (boxplot + points; log10)
make_Figure1_inset_trunc_per_gene <- function(constraint_gene, discordant_35, gene_trunc_ac) {
  
  # LOEUF<0.2 universe, one row per gene
  base <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(LOEUF), LOEUF < CFG$LOEUF_CUTOFF)
  
    disc_set <- unique(trimws(as.character(discordant_35$Gene)))
    disc_set <- disc_set[!is.na(disc_set) & disc_set != ""]
    
    # --- Clean & build base universe (no NA/blank Gene) ---
    base <- constraint_gene %>%
      dplyr::mutate(Gene = trimws(as.character(Gene))) %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(
        LOEUF = suppressWarnings(min(LOEUF, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      dplyr::filter(is.finite(LOEUF), LOEUF < CFG$LOEUF_CUTOFF)
    
    # --- Clean join table key too ---
    gt <- gene_trunc_ac %>%
      dplyr::mutate(Gene = trimws(as.character(Gene))) %>%
      dplyr::filter(!is.na(Gene), Gene != "") %>%
      dplyr::select(Gene, trunc_n_variants, trunc_AC)
    
    # --- Join + compute y + compute discordance deterministically ---
    df <- base %>%
      dplyr::left_join(gt, by = "Gene") %>%
      dplyr::mutate(
        y = suppressWarnings(as.numeric(trunc_n_variants))
      ) %>%
      dplyr::filter(is.finite(y), y >= 10) %>%
      dplyr::mutate(
        is_discordant = Gene %in% disc_set,
        group = dplyr::if_else(is_discordant, "Discordant", "Non-discordant")
      )
    
    top15_genes <- head(disc_set, 15)
    
    df <- df %>%
      dplyr::mutate(
        is_top15 = Gene %in% top15_genes
      )
    
    df$group <- factor(df$group, levels = c("Discordant", "Non-discordant"))
    
    # --- HARD FAIL if anything can create an NA x label ---
    if (anyNA(df$group)) {
      bad <- df %>% dplyr::filter(is.na(group) | is.na(Gene) | is.na(is_discordant))
      readr::write_csv(bad, file.path(CFG$OUT_DIR, "DEBUG_Fig2_BAD_ROWS.csv"))
      stop("NA detected in df$group. Wrote DEBUG_Fig2_BAD_ROWS.csv")
    }
    
    # --- ALWAYS write a small forensic bundle ---
    # pick a real directory even if CFG$OUT_DIR is unset
    out_dir <- CFG$OUT_DIR
    if (length(out_dir) == 0 || is.null(out_dir) || is.na(out_dir) || out_dir == "") {
      out_dir <- getwd()  # or tempdir()
      msg("[FIG2 DEBUG] CFG$OUT_DIR missing; using out_dir=%s", out_dir)
    }
    
    dbg_dir <- file.path(out_dir, "DEBUG_Fig2_bundle")
    msg("[FIG2 DEBUG] Writing bundle to: %s", dbg_dir)
    dir.create(dbg_dir, showWarnings = FALSE, recursive = TRUE)
    
    saveRDS(list(
      df = df,
      disc_set = disc_set,
      discordant_35_head = utils::head(discordant_35, 50),
      gene_trunc_ac_head = utils::head(gene_trunc_ac, 50),
      base_head = utils::head(base, 50)
    ), file.path(dbg_dir, "bundle.rds"))
    
    # --- Diagnose x variable exactly as ggplot sees it ---
    msg("[FIG2 DEBUG] x var class: %s", paste(class(df$group), collapse = ", "))
    msg("[FIG2 DEBUG] x var NA count: %d", sum(is.na(df$group)))
    msg("[FIG2 DEBUG] x var table:")
    print(table(df$group, useNA = "ifany"))
    
    # Common gotcha: literal "NA" string
    msg("[FIG2 DEBUG] Genes that are literal 'NA': %d", sum(df$Gene == "NA", na.rm = TRUE))
    
    # If anything could create an NA tick, STOP and write offenders
    if (anyNA(df$group)) {
      bad <- df %>%
        dplyr::mutate(
          Gene_dbg = dplyr::if_else(is.na(Gene), "<NA>", Gene),
          is_discordant_dbg = dplyr::if_else(is.na(is_discordant), NA_character_, as.character(is_discordant))
        ) %>%
        dplyr::filter(is.na(group) | is.na(Gene) | is.na(is_discordant) | Gene_dbg == "NA")
      
      readr::write_csv(bad, file.path(dbg_dir, "BAD_ROWS.csv"), na = "")
      msg("[FIG2 DEBUG] Wrote bad rows: %s", file.path(dbg_dir, "BAD_ROWS.csv"))
      stop("NA detected in df$group (x-axis). Inspect DEBUG_Fig2_bundle/BAD_ROWS.csv and bundle.rds")
    }
    
    # ---- FIG2 FORENSICS (add near the end, before ggplot) ----
    # =======================
    # FIG2 FORENSICS — DO NOT REMOVE
    # =======================
    dbg_dir <- file.path(getwd(), "DEBUG_Fig2_bundle")
    dir.create(dbg_dir, showWarnings = FALSE, recursive = TRUE)
    
    message("[FIG2 DEBUG] Writing forensic bundle to: ", dbg_dir)
    
    # What ggplot will see on the x axis
    message("[FIG2 DEBUG] x variable summary:")
    print(table(df$group, useNA = "ifany"))
    
    # Save full data frame
    saveRDS(df, file.path(dbg_dir, "df.rds"))
    readr::write_csv(df, file.path(dbg_dir, "df.csv"), na = "")
    
    # Catch *any* NA or literal 'NA' before plotting
    if (anyNA(df$group) || any(as.character(df$group) == "NA", na.rm = TRUE)) {
      bad <- df %>%
        dplyr::mutate(group_dbg = dplyr::if_else(is.na(group), "<NA>", as.character(group))) %>%
        dplyr::filter(is.na(group) | group_dbg == "NA")
      
      readr::write_csv(bad, file.path(dbg_dir, "BAD_ROWS.csv"), na = "")
      stop("FIG2 ERROR: NA (or literal 'NA') detected in x-axis variable. See DEBUG_Fig2_bundle/")
    }
    
  
    p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = y)) +
      
      # violin (behind everything)
      ggplot2::geom_violin(
        fill = "white",
        color = "black",
        linewidth = 1.2,
        trim = TRUE
      ) +
      
      # non-discordant genes (grey)
      ggplot2::geom_jitter(
        data = df %>% dplyr::filter(group == "Non-discordant"),
        width = 0.16, height = 0,
        size = 3.4, alpha = 0.45,
        color = "grey50"
      ) +
      
      # discordant genes NOT in top15 (black)
      ggplot2::geom_jitter(
        data = df %>% dplyr::filter(group == "Discordant", !is_top15),
        width = 0.16, height = 0,
        size = 3.8, alpha = 0.60,
        color = "grey55"
      ) +
      
      # discordant Top15 genes (red, outlined)
      ggplot2::geom_jitter(
        data = df %>% dplyr::filter(group == "Discordant", is_top15),
        width = 0.16, height = 0,
        size = 4.2, alpha = 0.95,
        shape = 21,
        fill = "red",
        color = "black",
        stroke = 0.6
      ) +
      
      # optional: small boxplot on top for median/IQR readability
      ggplot2::geom_boxplot(
        width = 0.20,
        outlier.shape = NA,
        
        fill  = "white",   # white box
        alpha = 0.55,      # <-- transparent-ish
        color = "black",
        linewidth = 0.9
      ) +
      
      ggplot2::scale_x_discrete(
        labels = c(
          "Discordant"     = "Discordant LOEUF < 0.2",
          "Non-discordant" = "Other LOEUF < 0.2"
        )
      ) +
      
      ggplot2::scale_y_log10(
        breaks = c(10, 30, 100, 300),
        labels = c("10", "30", "100", "300")
      ) +
      
      ggplot2::labs(
        x = NULL,
        y = "Truncating variants per gene (log10)"
      ) +
      
      ggplot2::theme_classic(base_size = 20) +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(face = "bold"),
        axis.text.x = ggplot2::element_text(size = 16),
        axis.text.y = ggplot2::element_text(size = 16)
      )
    
    p
}


# (4) Figure2_discordant_vs_rest.png  [MAIN]  (3-panel boxplots)
make_Figure2_discordant_vs_rest <- function(constraint_gene, discordant_35, gene_trunc_ac = NULL) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required. install.packages('patchwork')")
  }
  
  stopifnot(file.exists(CFG$CONSTRAINT_TSV))
  raw <- readr::read_tsv(CFG$CONSTRAINT_TSV, show_col_types = FALSE)
  
  # --- gene-level collapse from raw constraint TSV for the metrics we need ---
  # NOTE: do NOT canonical-filter here (it was zeroing you out)
  gene_metrics <- raw %>%
    dplyr::transmute(
      Gene   = norm_sym(.data[["gene"]]),
      lof_oe = suppressWarnings(as.numeric(.data[["lof.oe"]])),
      lo_ci  = suppressWarnings(as.numeric(.data[["lof.oe_ci.lower"]])),
      up_ci  = suppressWarnings(as.numeric(.data[["lof.oe_ci.upper"]]))
    ) %>%
    dplyr::filter(Gene != "", is.finite(lof_oe), is.finite(lo_ci), is.finite(up_ci)) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      lof_oe = median(lof_oe, na.rm = TRUE),
      lo_ci  = median(lo_ci,  na.rm = TRUE),
      up_ci  = median(up_ci,  na.rm = TRUE),
      .groups = "drop"
    )
  
  # --- make sure constraint_gene is 1 row per gene ---
  cg_gene <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF   = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      lof_obs = suppressWarnings(max(lof_obs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(LOEUF), !is.na(Gene), Gene != "")
  
  
  # --- join + define universe ---
  df <- cg_gene %>%
    dplyr::left_join(gene_metrics, by = "Gene") %>%
    dplyr::mutate(
      ci_width = pmax(up_ci - lo_ci, 1e-12),
      group = dplyr::if_else(Gene %in% discordant_35$Gene, "Discordant", "Non-discordant")
    ) %>%
    dplyr::filter(LOEUF < CFG$LOEUF_CUTOFF) %>%
    dplyr::filter(is.finite(lof_oe), is.finite(up_ci), is.finite(lo_ci))
  
  # ---- debug if empty (shouldn't be now, but keep it) ----
  if (nrow(df) == 0) {
    msg("[FIG2 DEBUG] df is empty after join+filters.")
    msg(" - nrow(cg_gene) = %d", nrow(cg_gene))
    msg(" - nrow(gene_metrics) = %d", nrow(gene_metrics))
    msg(" - LOEUF cutoff = %s", CFG$LOEUF_CUTOFF)
    stop("Figure2: df empty after join; check gene symbol harmonisation or input TSV.")
  }
  
  df$group <- factor(df$group, levels = c("Non-discordant", "Discordant"))
  tab <- table(df$group)
  
  msg("[FIG2 DEBUG] After filtering: Non-discordant=%d | Discordant=%d",
      as.integer(tab[["Non-discordant"]]), as.integer(tab[["Discordant"]]))
  
  if (any(tab == 0)) {
    stop("wilcox.test requires both groups; one group has 0 rows after filtering.")
  }
  
  n_nd <- as.integer(tab[["Non-discordant"]])
  n_d  <- as.integer(tab[["Discordant"]])
  
  df$group_lbl <- factor(
    ifelse(df$group == "Discordant",
           "Discordant",
           "Other LOEUF < 0.2"),
    levels = c("Other LOEUF < 0.2", "Discordant")
  )
  
  
  msg("[FIG2 DEBUG] group_lbl table:")
  print(table(df$group_lbl, useNA = "ifany"))
  
  top15 <- unique(discordant_35$Gene)[1:15]  # OR better: pass in a real top15 vector from Table 1
  df <- df %>%
    dplyr::mutate(
      pt_class = dplyr::if_else(Gene %in% top15, "Top 15 discordant", "All other genes")
    )
  df$pt_class <- factor(df$pt_class, levels = c("All other genes", "Top 15 discordant"))
  
  
  fmt_p <- function(p) {
    if (!is.finite(p)) return("Wilcoxon p = NA")
    if (p < 1e-16) return("Wilcoxon p < 1e-16")
    paste0("Wilcoxon p = ", format(p, digits = 3, scientific = TRUE))
  }
  
  add_panel <- function(ycol, title, ylab) {
    y <- df[[ycol]]
    y_pos <- suppressWarnings(min(y, na.rm = TRUE))
    pval <- suppressWarnings(stats::wilcox.test(y ~ df$group)$p.value)
    ptxt <- fmt_p(pval)
    
    # split once for clean layering
    df_other <- df[df$pt_class == "All other genes", , drop = FALSE]
    df_top15 <- df[df$pt_class == "Top 15 discordant", , drop = FALSE]
    
    ggplot2::ggplot() +
      
      # 1) background points (behind)
      ggplot2::geom_jitter(
        data = df_other,
        ggplot2::aes(x = group_lbl, y = .data[[ycol]]),
        width = 0.12, height = 0,
        size = 0.9, alpha = 0.18, color = "grey30"
      ) +
      
      # 2) violin outline
      ggplot2::geom_violin(
        data = df,
        ggplot2::aes(x = group_lbl, y = .data[[ycol]], group = group_lbl),
        fill = NA, color = "black",
        width = 0.85, trim = TRUE, linewidth = 0.6
      ) +
      
      # 3) IQR box (force constant transparency)
      ggplot2::geom_boxplot(
        data = df,
        inherit.aes = FALSE,
        ggplot2::aes(x = group_lbl, y = .data[[ycol]], group = group_lbl),
        width = 0.20,
        outlier.shape = NA,
        fill = "white", alpha = 0.70,          # <- consistent for every panel
        color = "black", linewidth = 0.95
      ) +
      
      # 4) median tick on top of box
      ggplot2::stat_summary(
        data = df,
        inherit.aes = FALSE,
        ggplot2::aes(x = group_lbl, y = .data[[ycol]], group = group_lbl),
        fun = median, geom = "point",
        shape = 95, size = 7, color = "black"
      ) +
      
      # 5) TOP15 red points LAST (on top of box/median)
      ggplot2::geom_jitter(
        data = df_top15,
        ggplot2::aes(x = group_lbl, y = .data[[ycol]]),
        width = 0.12, height = 0,
        size = 1.3, alpha = 1,
        shape = 21, fill = "red3", color = "black", stroke = 0.25
      ) +
      
      ggplot2::annotate(
        "text",
        x = Inf, y = y_pos,
        label = ptxt,
        hjust = 1.03, vjust = -0.6,
        size = 3.2,
        color = "black"
      )+
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(plot.margin = ggplot2::margin(5.5, 10, 5.5, 5.5))+
      ggplot2::labs(x = NULL, y = ylab) +
      ggplot2::theme_classic(base_size = 12) +
      ggplot2::theme(
        plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5, size = 18),
        plot.subtitle = ggplot2::element_blank(),
        axis.text.x   = ggplot2::element_text(size = 8),
        axis.text.y   = ggplot2::element_text(size = 10)
      )
    
  }
  
  
  pA <- add_panel("lof_oe",   "Fig 2A", "LoF o/e (point estimate)")
  pB <- add_panel("LOEUF", "Fig 2B", "LOEUF (upper CI bound)")
  
  pC <- add_panel("ci_width", "Fig 2C", "LoF o/e CI width") +
    ggplot2::scale_y_log10()
  
  pA <- pA + ggplot2::annotate("text", x = -Inf, y = Inf,
                               label = "A", hjust = -0.3, vjust = 1.4,
                               size = 6, fontface = "bold")
  
  pB <- pB + ggplot2::annotate("text", x = -Inf, y = Inf,
                               label = "B", hjust = -0.3, vjust = 1.4,
                               size = 6, fontface = "bold")
  
  pC <- pC + ggplot2::annotate("text", x = -Inf, y = Inf,
                               label = "C", hjust = -0.3, vjust = 1.4,
                               size = 6, fontface = "bold")
  
  
  pA + pB + pC + patchwork::plot_layout(ncol = 3)
}


# (5) Supplementary_Figure_NMD_escape_per_gene.png  [SUPP] (boxplot)
make_Supplementary_Figure_NMD_escape_per_gene <- function(nmd_df) {
  
  top15_genes <- c(
    "ARID1B","BIRC6","BTAF1","CACNA1C","CHD4","CIT","COL2A1",
    "KMT2C","MGA","PCLO","PRKDC","TRIO","TRIP12","WDFY3","ZFHX3"
  )
  
  # ----------------------------
  # Clean + normalise gene names
  # ----------------------------
  nmd_df <- nmd_df %>%
    dplyr::mutate(
      Gene  = toupper(trimws(as.character(Gene))),
      group = trimws(as.character(group))
    )
  
  top15_genes <- toupper(trimws(top15_genes))
  
  # ----------------------------
  # Define Top15 flag
  # ----------------------------
  nmd_df <- nmd_df %>%
    dplyr::mutate(
      is_top15 = (group == "Discordant") & (Gene %in% top15_genes),
      pt_class = dplyr::case_when(
        is_top15 ~ "Top15",
        group == "Discordant" ~ "Other_discordant",
        TRUE ~ "Background"
      )
    )
  
  # ----------------------------
  # Sanity check
  # ----------------------------
  n_top15_found <- sum(nmd_df$is_top15)
  message("[NMD FIG] Top15 genes matched in discordant group: ", n_top15_found)
  if (n_top15_found == 0) {
    warning("No Top15 genes were matched — check Gene symbols in nmd_df!")
    message("Top genes in nmd_df Discordant group are:")
    print(head(nmd_df %>% dplyr::filter(group == "Discordant") %>% dplyr::pull(Gene), 20))
  }
  
  # ----------------------------
  # Plot
  # ----------------------------
  p <- ggplot(nmd_df, aes(x = group, y = prop_escape)) +
    
    geom_violin(fill = "grey95", color = NA, trim = FALSE) +
    
    # Background jitter FIRST
    ggbeeswarm::geom_quasirandom(
      data = nmd_df %>% dplyr::filter(pt_class == "Background"),
      alpha = 0.25,
      size = 1.4,
      width = 0.18,
      color = "grey60"
    ) +
    
    # Other discordant (black)
    ggbeeswarm::geom_quasirandom(
      data = nmd_df %>% dplyr::filter(pt_class == "Other_discordant"),
      alpha = 0.35,
      size = 1.6,
      width = 0.18,
      color = "grey55"
    ) +
    
    # Top15 discordant (red)
    ggbeeswarm::geom_quasirandom(
      data = nmd_df %>% dplyr::filter(pt_class == "Top15"),
      alpha = 0.95,
      size = 2,
      width = 0.18,
      shape = 21,
      fill = "red",
      color = "black",
      stroke = 0.4
    ) +
    
    # Global boxplot (slightly transparent)
    geom_boxplot(
      width = 0.22,
      outlier.shape = NA,
      fill = "white",
      alpha = 0.65,
      color = "black",
      linewidth = 0.7
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)
    ) +
    
    labs(
      x = NULL,
      y = "Proportion of truncating variants predicted to escape NMD"
    ) +
    scale_x_discrete(
      labels = c(
        "Discordant"     = "Discordant LOEUF < 0.2",
        "Non-discordant" = "Other LOEUF < 0.2"
      )
    )+
    
    theme_classic(base_size = 18) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      axis.line  = element_line(linewidth = 0.8),
      axis.title.y = element_text(
        size = 12,        # smaller font
        face = "plain"    # not bold
    ),
    axis.text.x = element_text(size = 10))
  
  p
}

# (6) SuppFig_canonical_vs_min_LOEUF.png  [SUPP]
# Canonical vs minimum-transcript LOEUF per gene (screenshot-style)
make_SuppFig_canonical_vs_min_LOEUF <- function(constraint_tx = NULL,
                                                discordant_35 = NULL,
                                                constraint_gene = NULL,
                                                table1_top15 = NULL,
                                                ...) {
  # We deliberately read the TSV because it contains canonical + cds fields reliably
  stopifnot(file.exists(CFG$CONSTRAINT_TSV))
  raw <- readr::read_tsv(CFG$CONSTRAINT_TSV, show_col_types = FALSE)
  
  # ---- column detection ----
  if (!("gene" %in% names(raw))) stop("Constraint TSV missing 'gene' column.")
  if (!("canonical" %in% names(raw))) stop("Constraint TSV missing 'canonical' column.")
  
  # gnomAD v4 uses lof.oe_ci.upper for LOEUF (upper CI)
  loeuf_col <- dplyr::case_when(
    "lof.oe_ci.upper" %in% names(raw) ~ "lof.oe_ci.upper",
    "lof_oe_ci_upper" %in% names(raw) ~ "lof_oe_ci_upper",
    "lof_oe_ci_upper_bound" %in% names(raw) ~ "lof_oe_ci_upper_bound",
    TRUE ~ NA_character_
  )
  if (is.na(loeuf_col)) stop("Constraint TSV: cannot find LOEUF upper CI column (lof.oe_ci.upper).")
  
  tx <- raw %>%
    dplyr::transmute(
      Gene      = norm_sym(.data$gene),
      canonical = as.logical(.data$canonical),
      LOEUF     = suppressWarnings(as.numeric(.data[[loeuf_col]]))
    ) %>%
    dplyr::filter(Gene != "", is.finite(LOEUF))
  
  # ---- gene-level summary: canonical vs minimum across transcripts ----
  df <- tx %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      canonical_LOEUF = suppressWarnings(min(LOEUF[canonical %in% TRUE], na.rm = TRUE)),
      min_LOEUF       = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(canonical_LOEUF), is.finite(min_LOEUF))
  
  # ---- discordant gene set (accept several possible objects) ----
  disc_genes <- character(0)
  if (!is.null(discordant_35) && "Gene" %in% names(discordant_35)) disc_genes <- discordant_35$Gene
  if (!is.null(table1_top15)  && "Gene" %in% names(table1_top15))  disc_genes <- unique(c(disc_genes, table1_top15$Gene))
  disc_genes <- norm_sym(disc_genes)
  
  df <- df %>%
    dplyr::mutate(is_discordant = Gene %in% disc_genes)
  
  # ---- plot params to match screenshot ----
  cutoff <- CFG$LOEUF_CUTOFF %||% 0.2
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = canonical_LOEUF, y = min_LOEUF)) +
    ggplot2::geom_point(color = "grey70", alpha = 0.35, size = 2) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(is_discordant),
      ggplot2::aes(x = canonical_LOEUF, y = min_LOEUF),
      shape = 21, fill = "red", color = "black", size = 3.2, stroke = 0.6
    ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1.0, color = "black") +
    ggplot2::geom_vline(xintercept = cutoff, linetype = "dotted", linewidth = 1.0, color = "grey35") +
    ggplot2::annotate(
      "text",
      x = cutoff,
      y = 2.05,                    # near top of plotting range
      label = paste0("LOEUF = ", cutoff),
      angle = 90,                  # rotate text
      vjust = -0.4,                # nudge left/right relative to line
      hjust = 1,
      size = 3,
      color = "grey35"
    )+
    ggplot2::geom_hline(yintercept = cutoff, linetype = "dotted", linewidth = 1.0, color = "grey35") +
    ggplot2::annotate(
      "text",
      x = 2.05,                    # far right of plotting range
      y = cutoff,
      label = paste0("LOEUF = ", cutoff),
      hjust = 1,
      vjust = -0.6,                # sit just above the line
      size = 3,
      color = "grey35"
    )+
    ggplot2::coord_cartesian(xlim = c(0, 2.1), ylim = c(0, 2.1)) +
    ggplot2::labs(
      x = "Canonical transcript LOEUF",
      y = "Minimum LOEUF across transcripts"
    ) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "plain", size = 9),
      axis.text  = ggplot2::element_text(color = "black", size = 9),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
  
  p
}

make_SuppFig_LOEUF_vs_observed_truncations_top15 <- function(constraint_gene) {
  
  top15 <- c(
    "ARID1B","BIRC6","BTAF1","CACNA1C","CHD4","CIT","COL2A1","KMT2C",
    "MGA","PCLO","PRKDC","TRIO","TRIP12","WDFY3","ZFHX3"
  ) %>% norm_sym()
  
  cg_ultra_gene <- constraint_gene %>%
    dplyr::filter(LOEUF < CFG$LOEUF_CUTOFF) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF   = min(LOEUF, na.rm = TRUE),
      lof_obs = max(lof_obs, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(LOEUF), is.finite(lof_obs)) %>%
    dplyr::mutate(is_top15 = Gene %in% top15)
  
  # Sanity: all 15 present
  stopifnot(sum(cg_ultra_gene$is_top15) == 15)
  
  p <- ggplot2::ggplot(cg_ultra_gene, ggplot2::aes(x = LOEUF, y = lof_obs)) +
    ggplot2::geom_point(alpha = 0.4) +
    ggplot2::geom_point(
      data = cg_ultra_gene %>% dplyr::filter(is_top15),
      size = 2
    ) +
    ggrepel::geom_text_repel(
      data = cg_ultra_gene %>% dplyr::filter(is_top15),
      ggplot2::aes(label = Gene),
      max.overlaps = Inf
    ) +
    ggplot2::labs(
      x = "LOEUF (min per gene within LOEUF < 0.2 set)",
      y = "Observed LoF variants (lof_obs)"
    ) +
    ggplot2::theme_classic()
  
  p
}


# small helper if your script doesn't already have %||%
`%||%` <- function(a, b) if (!is.null(a)) a else b


# (7) SuppFig_length_adjusted_2panel.png  [SUPP]  -- screenshot-style
make_SuppFig_length_adjusted_2panel <- function(constraint_gene,
                                                discordant_35,
                                                table1_top15,
                                                return = c("combined","panels")) {
  return <- match.arg(return)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required. install.packages('patchwork')")
  }
  
  # ----------------------------
  # 1) Read constraint TSV for canonical CDS lengths
  # ----------------------------
  raw <- readr::read_tsv(CFG$CONSTRAINT_TSV, show_col_types = FALSE)
  
  if (!("gene" %in% names(raw))) stop("Constraint TSV missing 'gene' column.")
  if (!("cds_length" %in% names(raw))) stop("Constraint TSV missing 'cds_length' column.")
  if (!("canonical" %in% names(raw))) stop("Constraint TSV missing 'canonical' column.")
  
  len_df <- raw %>%
    dplyr::filter(.data$canonical %in% TRUE) %>%
    dplyr::transmute(
      Gene   = norm_sym(.data$gene),
      cds_bp = suppressWarnings(as.numeric(.data$cds_length))
    ) %>%
    dplyr::filter(is.finite(cds_bp), cds_bp > 0) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(cds_bp = max(cds_bp, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(cds_kb = cds_bp / 1000)
  
  # ----------------------------
  # 2) One row per gene from constraint_gene
  # ----------------------------
  cg <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF   = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      lof_obs = suppressWarnings(max(lof_obs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(LOEUF), is.finite(lof_obs))
  
  # ----------------------------
  # 3) Define discordant(35) vs background; then Top15 within discordant
  # ----------------------------
  disc35 <- norm_sym(discordant_35$Gene)
  top15  <- norm_sym(table1_top15$Gene)
  
  df <- cg %>%
    dplyr::inner_join(len_df, by = "Gene") %>%
    dplyr::filter(
      LOEUF < CFG$LOEUF_CUTOFF,
      lof_obs > 0,
      cds_kb > 0
    ) %>%
    dplyr::mutate(
      group = dplyr::if_else(Gene %in% disc35, "Discordant", "Other LOEUF<0.2"),
      trunc_per_kb = lof_obs / cds_kb,
      pt_class = dplyr::case_when(
        group != "Discordant" ~ "Background",
        Gene %in% top15       ~ "Top15",
        TRUE                  ~ "Other_discordant"
      )
    ) %>%
    dplyr::filter(is.finite(trunc_per_kb), trunc_per_kb > 0)
  
  # guard
  tab <- table(df$group)
  n_disc <- if ("Discordant" %in% names(tab)) unname(tab["Discordant"]) else 0
  n_rest <- if ("Other LOEUF<0.2" %in% names(tab)) unname(tab["Other LOEUF<0.2"]) else 0
  if (n_disc == 0 || n_rest == 0) {
    stop(
      "SuppFig_length_adjusted_2panel: need both groups after filtering.\n",
      "Counts: Discordant=", n_disc, " | Other LOEUF<0.2=", n_rest
    )
  }
  
  wt <- stats::wilcox.test(trunc_per_kb ~ group, data = df, exact = FALSE)
  p_txt <- formatC(wt$p.value, format = "e", digits = 2)
  
  # ✅ label only Top15
  label_genes <- top15
  
  # ----------------------------
  # Panel A
  # ----------------------------
  pA <- ggplot2::ggplot(df, ggplot2::aes(x = cds_kb, y = lof_obs)) +
    
    ggplot2::geom_point(
      data = df %>% dplyr::filter(pt_class == "Background"),
      color = "grey75", alpha = 0.45, size = 2.8
    ) +
    
    ggplot2::geom_point(
      data = df %>% dplyr::filter(pt_class == "Other_discordant"),
      shape = 21, fill = "grey55", color = "black",
      size = 4.0, stroke = 0.55
    ) +
    
    ggplot2::geom_point(
      data = df %>% dplyr::filter(pt_class == "Top15"),
      shape = 21, fill = "red", color = "black",
      size = 4.2, stroke = 0.6
    ) +
    
    ggrepel::geom_label_repel(
      data = df %>% dplyr::filter(Gene %in% label_genes),
      ggplot2::aes(label = Gene),
      size = 4.5,
      label.size = 0.25,
      fill = "white",
      color = "black",
      box.padding = 0.6,
      point.padding = 0.35,
      min.segment.length = 0,
      segment.color = "grey40",
      max.overlaps = Inf
    ) +
    
    ggplot2::scale_y_log10(
      breaks = c(1, 3, 10, 30, 100),
      labels = c("1", "3", "10", "30", "100")
    ) +
    
    ggplot2::labs(
      title = "A",
      x = "CDS length (kb)",
      y = "Observed truncating variants (log10)"
    ) +
    
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 22, hjust = 0),
      axis.title = ggplot2::element_text(face = "bold", size = 18),
      axis.text  = ggplot2::element_text(color = "black", size = 14)
    )
  
  # ----------------------------
  # Panel B
  # ----------------------------
  dfB <- df %>%
    dplyr::mutate(group = factor(group, levels = c("Discordant", "Other LOEUF<0.2")))
  
  pB <- ggplot2::ggplot(dfB, ggplot2::aes(x = group, y = trunc_per_kb)) +
    ggplot2::geom_violin(fill = "white", color = "black", linewidth = 0.9, trim = TRUE) +
    ggplot2::scale_x_discrete(
      labels = c(
        "Discordant"     = "Discordant LOEUF < 0.2",
        "Other LOEUF<0.2"= "Other LOEUF < 0.2"
      )
    ) +
    
    ggplot2::geom_jitter(
      data = dfB %>% dplyr::filter(pt_class == "Background"),
      color = "grey65", alpha = 0.35, width = 0.12, size = 2
    ) +
    
    ggplot2::geom_jitter(
      data = dfB %>% dplyr::filter(pt_class == "Other_discordant"),
      shape = 21, fill = "grey55", color = "black",
      stroke = 0.30, alpha = 0.85, width = 0.10, size = 3.0
    ) +
    
    ggplot2::geom_jitter(
      data = dfB %>% dplyr::filter(pt_class == "Top15"),
      shape = 21, fill = "red", color = "black",
      stroke = 0.40, alpha = 0.95, width = 0.10, size = 3.3
    ) +
    
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      alpha = 0.55,
      linewidth = 0.9,
      fatten = 2.5
    ) +
    
    ggplot2::scale_y_log10(
      breaks = c(0.3, 1, 3, 10),
      labels = c("0.3", "1.0", "3.0", "10.0")
    ) +
    
    ggplot2::labs(
      title = "B",
      subtitle = sprintf("Wilcoxon p = %s", p_txt),
      x = NULL,
      y = "Truncating variants per kb CDS (log10)"
    ) +
    
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 22, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0),
      axis.title    = ggplot2::element_text(face = "bold", size = 18),
      axis.text.x   = ggplot2::element_text(color = "black", size = 14),
      axis.text.y   = ggplot2::element_text(color = "black", size = 14)
    )
  
  combined <- pA + pB + patchwork::plot_layout(ncol = 2, widths = c(1, 1))
  
  if (return == "panels") return(list(A = pA, B = pB))
  combined
}


# (8) SuppFig_negative_control_syn_mis.png  [SUPP]
make_SuppFig_negative_control_syn_mis <- function(constraint_gene, discordant_35,
                                                  return = c("combined","panels")) {
  return <- match.arg(return)
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork required")
  
  # ---- helper: pick first existing column name from a set ----
  pick_col <- function(df, candidates) {
    hit <- candidates[candidates %in% names(df)]
    if (length(hit) == 0) return(NA_character_)
    hit[[1]]
  }
  
  # ---- detect columns (supports both dot + underscore styles) ----
  col_gene  <- pick_col(constraint_gene, c("Gene", "gene"))
  col_loeuf <- pick_col(constraint_gene, c("LOEUF", "loeuf", "lof.oe_ci.upper", "lof_oe_ci_upper"))
  col_syn   <- pick_col(constraint_gene, c("syn_obs", "syn.obs", "synonymous_obs", "synonymous.obs"))
  col_lof   <- pick_col(constraint_gene, c("lof_obs", "lof.obs", "trunc_obs", "trunc.obs"))
  col_mis   <- pick_col(constraint_gene, c("mis_oe", "mis.oe", "mis.oe_ci.upper", "mis_oe_ci_upper"))
  
  if (is.na(col_gene)  || is.na(col_loeuf)) stop("Need Gene + LOEUF columns in constraint_gene.")
  if (is.na(col_syn)) stop("No synonymous observed column found (syn_obs / syn.obs).")
  if (is.na(col_lof)) stop("No LoF/trunc observed column found (lof_obs / lof.obs).")
  if (is.na(col_mis)) stop("No missense o/e column found (mis_oe / mis.oe).")
  
  disc_set <- unique(trimws(as.character(discordant_35$Gene)))
  disc_set <- disc_set[!is.na(disc_set) & disc_set != ""]
  
  top15 <- unique(discordant_35$Gene)[1:15]
  top15 <- trimws(as.character(top15))
  top15 <- top15[!is.na(top15) & top15 != ""]
  
  # ---- one row per gene ----
  # NOTE: since constraint_gene is already gene-level in your current pipeline,
  # we keep first row per gene after filtering rather than using max() across duplicates.
  df <- constraint_gene %>%
    dplyr::transmute(
      Gene    = trimws(as.character(.data[[col_gene]])),
      LOEUF   = suppressWarnings(as.numeric(.data[[col_loeuf]])),
      syn_obs = suppressWarnings(as.numeric(.data[[col_syn]])),
      lof_obs = suppressWarnings(as.numeric(.data[[col_lof]])),
      mis_oe  = suppressWarnings(as.numeric(.data[[col_mis]]))
    ) %>%
    dplyr::filter(
      !is.na(Gene), Gene != "",
      is.finite(LOEUF), LOEUF < CFG$LOEUF_CUTOFF,
      is.finite(syn_obs), syn_obs >= 0,
      is.finite(lof_obs), lof_obs >= 0,
      is.finite(mis_oe)
    ) %>%
    dplyr::arrange(Gene, LOEUF) %>%
    dplyr::group_by(Gene) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      is_disc = Gene %in% disc_set,
      plot_class = dplyr::case_when(
        Gene %in% top15 ~ "Top15",
        is_disc         ~ "Other_discordant",
        TRUE            ~ "Background"
      ),
      # ---- key negative-control quantity (synonymous-adjusted LoF excess) ----
      # residual of log10(lof_obs+1) ~ log10(syn_obs+1)
      lof_adj_by_syn = as.numeric(stats::residuals(
        stats::lm(log10(lof_obs + 1) ~ log10(syn_obs + 1))
      ))
    )
  
  # ---- HARD sanity checks ----
  msg("[SuppFig_neg_ctrl] rows plotted (should == #genes) = %d", nrow(df))
  msg("[SuppFig_neg_ctrl] unique genes = %d", dplyr::n_distinct(df$Gene))
  msg("[SuppFig_neg_ctrl] discordant genes in universe = %d", dplyr::n_distinct(df$Gene[df$is_disc]))
  
  # write debug so you can inspect exactly what got plotted
  dbg <- file.path(getwd(), "DEBUG_SuppFig_negative_control_syn_mis_genelevel.csv")
  readr::write_csv(df, dbg, na = "")
  msg("[SuppFig_neg_ctrl] wrote plotted table: %s", dbg)
  
  # NOTE: you said you've removed printing the p-value on the plot.
  # We also no longer compute a Wilcoxon p-value on raw syn_obs, because it is confounded by gene size/callability.
  # If you still want a number for QC logs only, uncomment:
  # wtA <- suppressWarnings(stats::wilcox.test(lof_adj_by_syn ~ is_disc, data = df, exact = FALSE))
  # msg("[SuppFig_neg_ctrl] Wilcoxon p (LoF residual | syn) = %s", formatC(wtA$p.value, format="e", digits=2))
  
  # ----------------------------
  # Panel A: Synonymous-adjusted LoF excess vs LOEUF
  # ----------------------------
  pA <- ggplot2::ggplot(df, ggplot2::aes(x = LOEUF, y = lof_adj_by_syn)) +
    
    ggplot2::geom_point(
      data = df %>% dplyr::filter(plot_class == "Background"),
      color = "grey75", alpha = 0.45, size = 2.3
    ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(plot_class == "Other_discordant"),
      color = "grey55", alpha = 0.80, size = 2.6
    ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(plot_class == "Top15"),
      shape = 21, fill = "red", color = "black",
      size = 3.4, stroke = 0.5, alpha = 0.95
    ) +
    
    ggplot2::labs(
      title = "A",
      x = "LOEUF",
      y = "LoF burden residual after adjusting for synonymous burden"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title    = ggplot2::element_text(size = 11),
      axis.text     = ggplot2::element_text(size = 10)
    )
  
  # ----------------------------
  # Panel B: Missense tolerance vs truncating burden
  # ----------------------------
  pB <- ggplot2::ggplot(df, ggplot2::aes(x = mis_oe, y = lof_obs + 1)) +
    
    ggplot2::geom_point(
      data = df %>% dplyr::filter(plot_class == "Background"),
      color = "grey75", alpha = 0.45, size = 2.3
    ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(plot_class == "Other_discordant"),
      color = "grey55", alpha = 0.80, size = 2.6
    ) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(plot_class == "Top15"),
      shape = 21, fill = "red", color = "black",
      size = 3.4, stroke = 0.5, alpha = 0.95
    ) +
    
    ggplot2::scale_y_log10(breaks = c(1, 3, 10, 30), labels = c("1", "3", "10", "30")) +
    ggplot2::labs(
      title = "B",
      x = "Missense observed/expected",
      y = "Observed truncating variants per gene (log10)"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title    = ggplot2::element_text(size = 11),
      axis.text     = ggplot2::element_text(size = 10)
    )
  
  pA + pB + patchwork::plot_layout(ncol = 2)
  
  if (return == "panels") return(list(A = pA, B = pB))
  pA + pB + patchwork::plot_layout(ncol = 2)
}



# ============================================================
# ADD: new figure builder
# Place this with the other "Figure builders" (e.g., after make_SuppFig_negative_control_syn_mis)
# ============================================================
make_SuppFig_trunc_burden_adjusted_by_syn_obs <- function(constraint_gene, discordant_35) {
  
  disc_set <- norm_sym(discordant_35$Gene)
  top15 <- head(disc_set, 15)
  
  # one row per gene (protect against transcript-level inputs)
  df <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF   = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      lof_obs = suppressWarnings(max(lof_obs, na.rm = TRUE)),
      syn_obs = suppressWarnings(max(syn_obs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      !is.na(Gene), Gene != "",
      is.finite(LOEUF), LOEUF < CFG$LOEUF_CUTOFF,
      is.finite(lof_obs), is.finite(syn_obs)
    ) %>%
    dplyr::mutate(
      is_discordant = Gene %in% disc_set,
      group = dplyr::if_else(is_discordant, "Discordant", "Other LOEUF < 0.2")
    )
  
  # guard: need both groups
  tab <- table(df$group)
  if (length(tab) != 2 || any(tab == 0)) {
    stop("Need both groups with >=1 row. Got: ", paste(names(tab), tab, sep = "=", collapse = " | "))
  }
  
  # regress out synonymous burden (log1p handles zeros)
  fit <- stats::lm(log1p(lof_obs) ~ log1p(syn_obs), data = df)
  df$lof_resid <- as.numeric(stats::resid(fit))
  
  # drop non-finite residuals (prevents violin warnings)
  df <- df %>% dplyr::filter(is.finite(lof_resid))
  
  # Wilcoxon on residuals
  wt <- stats::wilcox.test(lof_resid ~ group, data = df, exact = FALSE)
  p_txt <- formatC(wt$p.value, format = "e", digits = 2)
  
  # factor order
  df <- df %>%
    dplyr::mutate(
      is_discordant = Gene %in% disc_set,
      is_top15      = Gene %in% top15,
      
      group = dplyr::if_else(
        is_discordant,
        "Discordant",
        "Other LOEUF < 0.2"
      )
    )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = lof_resid)) +
    
    ggplot2::geom_violin(fill = "white", color = "black", linewidth = 0.9, trim = TRUE) +
    
    ggplot2::geom_jitter(
      data = df %>% dplyr::filter(group == "Other LOEUF < 0.2"),
      color = "grey40",
      alpha = 0.35,
      width = 0.12,
      size = 2.0
    ) +
    # --- Discordant but NOT Top15 (grey55)
    ggplot2::geom_jitter(
      data = df %>% dplyr::filter(group == "Discordant", !is_top15),
      shape  = 21,
      fill   = "grey55",
      color  = "black",
      stroke = 0.30,
      alpha  = 0.85,
      width  = 0.10,
      size   = 3.0
    ) +
    
    # --- Top15 discordant (red)
    ggplot2::geom_jitter(
      data = df %>% dplyr::filter(group == "Discordant", is_top15),
      shape  = 21,
      fill   = "red",
      color  = "black",
      stroke = 0.40,
      alpha  = 0.95,
      width  = 0.10,
      size   = 3.6
    ) +
    
    
    ggplot2::geom_boxplot(
      width = 0.18,
      outlier.shape = NA,
      fill = "white",
      alpha = 0.55,
      color = "black",
      linewidth = 0.9,
      fatten = 2.5
    ) +
    
    ggplot2::scale_x_discrete(
      labels = c(
        "Discordant"        = "Discordant LOEUF < 0.2",
        "Other LOEUF < 0.2" = "Other LOEUF < 0.2"
      )
    ) +
    
    ggplot2::labs(
      subtitle = sprintf("Wilcoxon p = %s", p_txt),
      x = NULL,
      y = "LoF burden adjusted for synonymous burden"
    ) +
    ggplot2::theme_classic(base_size = 18) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0),
      axis.title    = ggplot2::element_text(face = "bold", size = 18),
      axis.text.x   = ggplot2::element_text(color = "black", size = 14),
      axis.text.y   = ggplot2::element_text(color = "black", size = 14)
    )
  
  p
}



# (9) SuppFig_perm_null_2panel.png  [SUPP]

make_SuppFig_perm_null_2panel <- function(constraint_gene,
                                          discordant_table1,
                                          constraint_tx = NULL,
                                          reps = 1000,
                                          seed = 1,
                                          return = c("combined", "panels"),
                                          ...) {
  
  return <- match.arg(return)
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required. install.packages('patchwork')")
  }
  
  stopifnot(file.exists(CFG$CONSTRAINT_TSV))
  raw <- readr::read_tsv(CFG$CONSTRAINT_TSV, show_col_types = FALSE)
  
  # ---- canonical CDS length per gene (bp -> kb) ----
  if (!("gene" %in% names(raw))) stop("Constraint TSV missing 'gene' column.")
  if (!("cds_length" %in% names(raw))) stop("Constraint TSV missing 'cds_length' column.")
  if (!("canonical" %in% names(raw))) stop("Constraint TSV missing 'canonical' column.")
  
  len_df <- raw %>%
    dplyr::filter(.data$canonical %in% TRUE) %>%
    dplyr::transmute(
      Gene   = norm_sym(.data$gene),
      cds_bp = suppressWarnings(as.numeric(.data$cds_length))
    ) %>%
    dplyr::filter(Gene != "", is.finite(cds_bp), cds_bp > 0) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(cds_bp = max(cds_bp, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(cds_kb = cds_bp / 1000)
  
  # ---- gene-level LOEUF + truncating burden ----
  cg_gene <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF   = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      lof_obs = suppressWarnings(max(lof_obs, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(Gene != "", is.finite(LOEUF), is.finite(lof_obs), lof_obs > 0)
  
  df <- cg_gene %>%
    dplyr::inner_join(len_df, by = "Gene") %>%
    dplyr::filter(is.finite(cds_kb), cds_kb > 0) %>%
    dplyr::mutate(trunc_per_kb = lof_obs / cds_kb) %>%
    dplyr::filter(is.finite(trunc_per_kb))
  
  # ---- Universe: LOEUF < cutoff ----
  universe <- df %>%
    dplyr::filter(LOEUF < CFG$LOEUF_CUTOFF)
  
  if (nrow(universe) < 10) {
    stop("Universe too small after LOEUF cutoff; check CFG$LOEUF_CUTOFF and inputs.")
  }
  
  # ---- Observed set: ALL discordant genes (N should be 35) ----
  obs <- universe %>%
    dplyr::filter(Gene %in% norm_sym(discordant_table1$Gene))
  
  N <- nrow(obs)
  if (N < 2) {
    stop("Observed discordant set too small after filtering (N<2). Check discordant_table1$Gene naming/normalisation.")
  }
  
  obs_median <- stats::median(obs$trunc_per_kb, na.rm = TRUE)
  
  # ---- Null: sample N genes from universe, compute median trunc_per_kb ----
  set.seed(seed)
  null_medians <- replicate(reps, {
    samp <- universe %>% dplyr::slice_sample(n = N, replace = FALSE)
    stats::median(samp$trunc_per_kb, na.rm = TRUE)
  })
  
  # Empirical one-sided p (>= observed)
  p_emp <- (1 + sum(null_medians >= obs_median, na.rm = TRUE)) / (reps + 1)
  p_txt <- formatC(p_emp, format = "f", digits = 6)
  
  subtitle <- sprintf("N = %d | reps = %d | empirical one-sided p = %s", N, reps, p_txt)
  
  # ----------------------------
  # Text sizing (IMPORTANT)
  # ----------------------------
  axis_text_pt  <- 9
  axis_title_pt <- 10
  subtitle_pt   <- 10
  
  # ggplot annotation size is NOT points; convert points -> ggplot "size"
  ann_size <- axis_text_pt / ggplot2::.pt
  
  # line widths (thin dashed)
  vline_lw <- 0.7
  curve_lw <- 1.0
  hist_lw  <- 0.3
  
  # ---- Panel A ----
  df_null <- data.frame(null_median = null_medians)
  x_min <- min(df_null$null_median, obs_median, na.rm = TRUE)
  x_max <- max(df_null$null_median, obs_median, na.rm = TRUE)
  x_range <- x_max - x_min
  x_nudge <- 0.02 * x_range
  
  pA <- ggplot2::ggplot(df_null, ggplot2::aes(x = null_median)) +
    ggplot2::geom_histogram(
      bins = 60,
      fill = "grey90",
      color = "grey40",
      linewidth = hist_lw
    ) +
    ggplot2::geom_vline(
      xintercept = obs_median,
      linewidth = vline_lw,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = obs_median + x_nudge,
      y = Inf,
      label = "Observed median",
      angle = 90,
      vjust = 1.1,
      hjust = 1,
      size = ann_size,
      fontface = "plain"
    ) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    ggplot2::labs(
      title = "A",
      x = "Median truncating variants per kb CDS",
      y = "Number of resamples"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
      axis.title = ggplot2::element_text(face = "plain", size = axis_title_pt),
      axis.text  = ggplot2::element_text(size = axis_text_pt)
    )
  
  # ---- Panel B ----
  ecdf_fun <- stats::ecdf(null_medians)
  ecdf_at_obs <- as.numeric(ecdf_fun(obs_median))
  right_tail  <- 1 - ecdf_at_obs
  
  dens  <- stats::density(null_medians, n = 2048)
  cdf_y <- cumsum(dens$y) / sum(dens$y)
  df_cdf <- data.frame(x = dens$x, cdf = cdf_y)
  cdf_at_obs_smooth <- stats::approx(df_cdf$x, df_cdf$cdf, xout = obs_median, rule = 2)$y
  
  ecdf_label <- sprintf("ECDF = %.3f | Right tail = %.3f", ecdf_at_obs, right_tail)
  
  pB <- ggplot2::ggplot(df_cdf, ggplot2::aes(x = x, y = cdf)) +
    ggplot2::geom_line(linewidth = curve_lw, color = "black") +
    ggplot2::geom_vline(
      xintercept = obs_median,
      linewidth = vline_lw,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_point(
      data = data.frame(x = obs_median, y = cdf_at_obs_smooth),
      ggplot2::aes(x = x, y = y),
      inherit.aes = FALSE,
      size = 2.6
    ) +
    ggplot2::annotate(
      "text",
      x = obs_median + (2 * x_nudge),
      y = 0.68,
      label = ecdf_label,
      angle = 90,
      size = ann_size,
      fontface = "plain"
    ) +
    ggplot2::labs(
      title = "B",
      x = "Median truncating variants per kb CDS",
      y = "Cumulative probability"
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1.05), clip = "off") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title  = ggplot2::element_text(face = "bold", size = 16, hjust = 0),
      axis.title  = ggplot2::element_text(face = "plain", size = axis_title_pt),
      axis.text   = ggplot2::element_text(size = axis_text_pt),
      plot.margin = ggplot2::margin(t = 10, r = 40, b = 10, l = 10)
    )
  
  # If the caller wants the panels (for super-fig labelling), return them here
  if (return == "panels") return(list(A = pA, B = pB))
  
  # Default: return combined (exactly as before)
  (pA + pB + patchwork::plot_layout(ncol = 2)) +
    patchwork::plot_annotation(
      title = NULL,
      subtitle = subtitle,
      theme = ggplot2::theme(
        plot.subtitle = ggplot2::element_text(size = subtitle_pt, hjust = 0),
        plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10)
      )
    )
}


#insert updated superfig
make_SuperFig_all_key_supp <- function(
    nmd_df,
    constraint_gene,
    discordant_35,
    reps = 1000,
    seed = 1
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork required")
  }
  
  # ------------------------------------------------------------
  # Local panel tagger (ggplot-safe)
  # ------------------------------------------------------------
  tag_panel <- function(p, label) {
    p + ggplot2::annotate(
      "text",
      x = -Inf, y = Inf,
      label = label,
      hjust = -0.05,
      vjust = 1.25,
      fontface = "bold",
      size = 3.9
    )
  }
  
  # ------------------------------------------------------------
  # Build plots (CALL existing functions so edits propagate)
  # ------------------------------------------------------------
  pA <- make_Supplementary_Figure_NMD_escape_per_gene(nmd_df)
  
  neg_pan <- make_SuppFig_negative_control_syn_mis(
    constraint_gene, discordant_35,
    return = "panels"
  )
  pB <- neg_pan$A
  pC <- neg_pan$B
  
  pD <- make_SuppFig_canonical_vs_min_LOEUF(
    discordant_35 = discordant_35
  )
  
  perm_pan <- make_SuppFig_perm_null_2panel(
    constraint_gene,
    discordant_table1 = discordant_35,
    reps = reps,
    seed = seed,
    return = "panels"
  )
  pE <- perm_pan$A
  pF <- perm_pan$B
  
  # ------------------------------------------------------------
  # Shrink blobs slightly (superfig only)
  # ------------------------------------------------------------
  if (exists("shrink_for_superfig", mode = "function")) {
    pA <- shrink_for_superfig(pA, point_scale = 0.55, text_scale = 0.88)
    pB <- shrink_for_superfig(pB, point_scale = 0.60, text_scale = 0.88)
    pC <- shrink_for_superfig(pC, point_scale = 0.60, text_scale = 0.88)
    pD <- shrink_for_superfig(pD, point_scale = 0.60, text_scale = 0.88)
    pE <- shrink_for_superfig(pE, point_scale = 0.60, text_scale = 0.88)
    pF <- shrink_for_superfig(pF, point_scale = 0.60, text_scale = 0.88)
  }
  
  # ------------------------------------------------------------
  # Uniform axis sizing (SAFE: applied per-panel, not via patchwork merge)
  # ------------------------------------------------------------
  uniform_axis_theme <- ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = 10, face = "plain"),
    axis.title.y = ggplot2::element_text(size = 10, face = "plain"),
    axis.text.x  = ggplot2::element_text(size = 9),
    axis.text.y  = ggplot2::element_text(size = 9),
    plot.subtitle = ggplot2::element_text(size = 9)
  )
  
  pA <- pA + uniform_axis_theme
  pB <- pB + uniform_axis_theme
  pC <- pC + uniform_axis_theme
  pD <- pD + uniform_axis_theme
  pE <- pE + uniform_axis_theme
  pF <- pF + uniform_axis_theme
  
  # Remove internal titles (superfig handles labelling)
  pA <- pA + ggplot2::labs(title = NULL)
  pB <- pB + ggplot2::labs(title = NULL)
  pC <- pC + ggplot2::labs(title = NULL)
  pD <- pD + ggplot2::labs(title = NULL)
  pE <- pE + ggplot2::labs(title = NULL)
  pF <- pF + ggplot2::labs(title = NULL)
  
  # ------------------------------------------------------------
  # Apply A–F labels in Results order:
  # A syn ctrl, B missense, C canonical/min, D perm hist, E perm ECDF, F NMD
  # ------------------------------------------------------------
  
  pA_new <- tag_panel(pB, "A")
  pB_new <- tag_panel(pC, "B")
  pC_new <- tag_panel(pD, "C")
  pD_new <- tag_panel(pE, "D")
  pE_new <- tag_panel(pF, "E")
  pF_new <- tag_panel(pA, "F")
  
  # ------------------------------------------------------------
  # Layout: same geometry, reordered panels
  # ------------------------------------------------------------
  # ------------------------------------------------------------
  # Layout: A B C on top; D E F on bottom
  # ------------------------------------------------------------
  super <- patchwork::wrap_plots(
    pA_new, pB_new, pC_new,
    pD_new, pE_new, pF_new,
    ncol = 3,
    byrow = TRUE
  ) +
    patchwork::plot_layout(nrow = 2)
  
  
  super <- super & ggplot2::theme(
    plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 12)
  )
  
  return(super)
}





make_SuperFig_discordant_summary <- function(
    constraint_gene,
    discordant_35,
    gene_trunc_ac,
    table1_top15
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork required")
  
  # ----------------------------
  # 1) Build plots FIRST
  # ----------------------------
  pTop <- make_Figure_LOEUF_vs_trunc_AC_discordant_highlight(
    constraint_gene, discordant_35, gene_trunc_ac, table1_top15
  )
  
  pB <- make_Figure1_inset_trunc_per_gene(
    constraint_gene, discordant_35, gene_trunc_ac
  )
  
  pC <- make_SuppFig_trunc_burden_adjusted_by_syn_obs(
    constraint_gene, discordant_35
  )
  
  # pD: we need the TWO panels (A and B) separately for D/E tagging
  if ("return" %in% names(formals(make_SuppFig_length_adjusted_2panel))) {
    pD_panels <- make_SuppFig_length_adjusted_2panel(
      constraint_gene, discordant_35, table1_top15, return = "panels"
    )
    pD_left  <- pD_panels$A
    pD_right <- pD_panels$B
    
    # after shrink_for_superfig(...)
    
    pC       <- move_subtitle_to_bottom(pC,       size = 3.0)
    pD_left  <- move_subtitle_to_bottom(pD_left,  size = 3.0)
    pD_right <- move_subtitle_to_bottom(pD_right, size = 3.0)
    
  } else {
    stop(
      "Your make_SuppFig_length_adjusted_2panel() does not support return='panels'.\n",
      "Add a 'return' argument to that function so it can return list(A=pA, B=pB) when return=='panels'."
    )
  }
  
  # ----------------------------
  # 2) Shrink ONLY for super-fig (does not change individual figures)
  # ----------------------------
  pTop <- shrink_for_superfig(pTop, point_scale = 0.65, text_scale = 0.85, y_title_scale = 0.75)
  pB   <- shrink_for_superfig(pB,   point_scale = 0.55, text_scale = 0.85, y_title_scale = 0.75) # smaller points just for B
  pC   <- shrink_for_superfig(pC,   point_scale = 0.65, text_scale = 0.85, y_title_scale = 0.75)
  pD_left  <- shrink_for_superfig(pD_left,  point_scale = 0.65, text_scale = 0.85, y_title_scale = 0.75)
  pD_right <- shrink_for_superfig(pD_right, point_scale = 0.65, text_scale = 0.85, y_title_scale = 0.75)
  
  # ----------------------------
  # 3) Remove plot titles (keep subtitles like Wilcoxon p)
  # ----------------------------
  pTop     <- pTop     + ggplot2::labs(title = NULL)
  pB       <- pB       + ggplot2::labs(title = NULL)
  pC       <- pC       + ggplot2::labs(title = NULL)
  pD_left  <- pD_left  + ggplot2::labs(title = NULL)
  pD_right <- pD_right + ggplot2::labs(title = NULL)
  
  # ----------------------------
  # 3b) Uniform, readable axes for super-figure only (slightly larger; not bold)
  # ----------------------------
  pTop     <- apply_superfig_axes(pTop,     axis_title_size = 11, axis_text_size = 10)
  pB       <- apply_superfig_axes(pB,       axis_title_size = 11, axis_text_size = 10)
  pC       <- apply_superfig_axes(pC,       axis_title_size = 11, axis_text_size = 10)
  pD_left  <- apply_superfig_axes(pD_left,  axis_title_size = 11, axis_text_size = 10)
  pD_right <- apply_superfig_axes(pD_right, axis_title_size = 11, axis_text_size = 10)
  
  # ----------------------------
  # 4) Define highlighting sets
  # ----------------------------
  
  
  
  # ----------------------------
  # 4) Add panel labels (A–E) with original names in parentheses
  # ----------------------------
  pTop <- add_panel_label(pTop, "A")
  pB   <- add_panel_label(pB,   "B")
  pC   <- add_panel_label(pC,   "C")
  pD_left  <- add_panel_label(pD_left,  "D")
  pD_right <- add_panel_label(pD_right, "E")
  
  # ----------------------------
  # 5) Combine layout (same as before)
  # ----------------------------
  pD <- pD_left | pD_right
  
  super_fig <- (pTop / (pB | pC) / pD) +
    patchwork::plot_layout(heights = c(1.15, 1, 1))
  
  # ----------------------------
  # 6) Global styling for the super-fig only
  # ----------------------------
  super_fig <- super_fig &
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 18, r = 22, b = 18, l = 22),
      plot.subtitle = ggplot2::element_text(margin = ggplot2::margin(t = 12))
    )
  df <- constraint_gene %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      LOEUF = suppressWarnings(min(LOEUF, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(Gene), Gene != "", is.finite(LOEUF))
  return(super_fig)
}






# ----------------------------
# Master build + dispatcher
# ----------------------------
build_all_figures <- function(mode=c("all","plots-only","one"), force=FALSE, fig=NULL) {
  mode <- match.arg(mode)
  msg("[INFO] Working directory: %s", getwd())
  msg("[INFO] Mode: %s | force=%s", mode, force)

  constraint_gene <- cache_get_or_build("constraint_gene", force=force, build_fn=build_constraint_gene)
  discordant_35   <- cache_get_or_build("discordant_35", force=force, build_fn=build_discordant_genes)
  table1_top15    <- cache_get_or_build("table1_top15", force=force, build_fn=build_table1)
  gene_trunc_ac   <- cache_get_or_build("gene_trunc_ac", force=force, build_fn=build_gene_trunc_ac)
  constraint_tx   <- cache_get_or_build("constraint_tx", force=force, build_fn=build_constraint_tx)

  # NMD cache: if you have a real NMD table, replace here. For now, derive placeholder from lof_obs distribution.
  nmd_df <- cache_get_or_build("nmd_per_gene", force=force, build_fn=function() {
    
    # Prefer the true per-gene table if it exists
    if (file.exists("Supplementary_Figure_NMD_escape_per_gene.png") &&
        file.exists("nmd_heuristic_truncations.csv")) {
      message("[NMD] Using nmd_heuristic_truncations.csv to build per-gene prop_escape")
    }
    
    if (file.exists("nmd_heuristic_truncations.csv")) {
      v <- readr::read_csv("nmd_heuristic_truncations.csv", show_col_types = FALSE)
      
      # gene_symbol + group + nmd_escape are present per your grep output
      stopifnot(all(c("gene_symbol","group","nmd_escape") %in% names(v)))
      
      gene_summ <- v %>%
        dplyr::filter(!is.na(gene_symbol), gene_symbol != "", !is.na(group)) %>%
        dplyr::mutate(
          Gene = norm_sym(gene_symbol),
          group = dplyr::if_else(group == "Discordant", "Discordant", "Non-discordant")
        ) %>%
        dplyr::group_by(Gene, group) %>%
        dplyr::summarise(
          n_trunc = dplyr::n(),
          prop_escape = mean(nmd_escape %in% TRUE, na.rm = TRUE),
          .groups = "drop"
        )
      
      # HARD GUARD: should be <= 35 discordant genes
      message("[NMD] Genes per group:")
      print(dplyr::count(gene_summ, group))
      
      message("[NMD] Max rows per gene check (should all be 1):")
      print(gene_summ %>% dplyr::count(group, Gene) %>% dplyr::arrange(desc(n)) %>% head(10))
      
      return(gene_summ %>% dplyr::select(Gene, group, prop_escape))
    }
    
    stop("Cannot build NMD per-gene summary: missing nmd_heuristic_truncations.csv")
  })
  
  

  figs <- list(
    Figure_LOEUF_vs_truncating_clean_Table1Top15label = function() {
      make_Figure_LOEUF_vs_truncating_clean_Table1Top15label(constraint_gene, table1_top15)
    },
    
    Figure_LOEUF_vs_trunc_AC_discordant_highlight = function() {
      make_Figure_LOEUF_vs_trunc_AC_discordant_highlight(
        constraint_gene, discordant_35, gene_trunc_ac, table1_top15
      )
    },
    
    Figure1_inset_trunc_per_gene = function() {
      make_Figure1_inset_trunc_per_gene(constraint_gene, discordant_35, gene_trunc_ac)
    },
    
    Figure2_discordant_vs_rest = function() {
      make_Figure2_discordant_vs_rest(constraint_gene, discordant_35, gene_trunc_ac)
    },
    
    Supplementary_Figure_NMD_escape_per_gene = function() {
      make_Supplementary_Figure_NMD_escape_per_gene(nmd_df)
    },
    
    SuppFig_canonical_vs_min_LOEUF = function() {
      make_SuppFig_canonical_vs_min_LOEUF(discordant_35 = discordant_35, table1_top15 = table1_top15)
    },
    
    SuppFig_LOEUF_vs_observed_truncations_top15 = function() {
      make_SuppFig_LOEUF_vs_observed_truncations_top15(constraint_gene)
    },
    
    SuppFig_length_adjusted_2panel = function() {
      make_SuppFig_length_adjusted_2panel(constraint_gene, discordant_35, table1_top15)
    },
    
    SuppFig_negative_control_syn_mis = function() {
      make_SuppFig_negative_control_syn_mis(constraint_gene, discordant_35)
    },
    
    SuppFig_perm_null_2panel = function() {
      make_SuppFig_perm_null_2panel(constraint_gene, discordant_35, reps = 1000, seed = CFG$SEED)
    },
    
    # NEW: trunc burden adjusted by synonymous burden
    SuppFig_trunc_burden_adjusted_by_syn_obs = function() {
      make_SuppFig_trunc_burden_adjusted_by_syn_obs(constraint_gene, discordant_35)
    },
    
    SuperFig_discordant_summary = function() {
      make_SuperFig_discordant_summary(
        constraint_gene,
        discordant_35,
        gene_trunc_ac,
        table1_top15
      )
    },
    
    SuperFig_all_key_supp = function() {
      make_SuperFig_all_key_supp(
        nmd_df         = nmd_df,
        constraint_gene = constraint_gene,
        discordant_35   = discordant_35,
        reps            = 1000,
        seed            = CFG$SEED
      )
    }
    
  )
  
  out_path_for <- function(name) {
    main_names <- c(
      "Figure_LOEUF_vs_truncating_clean_Table1Top15label",
      "Figure_LOEUF_vs_trunc_AC_discordant_highlight",
      "Figure1_inset_trunc_per_gene",
      "Figure2_discordant_vs_rest",
      "SuperFig_discordant_summary"
    )
    
    if (name %in% main_names) {
      file.path(CFG$OUT_MAIN_DIR, paste0(name, ".png"))
    } else {
      file.path(CFG$OUT_SUPP_DIR, paste0(name, ".png"))
    }
  }
  
  build_one <- function(name) {
    if (!name %in% names(figs)) stop("Unknown fig target: ", name)
    p <- figs[[name]]()
    out <- out_path_for(name)
    
    # Give the super-figure more space
if (name == "SuperFig_discordant_summary") {
  w <- 12; h <- 14
} else if (name == "Figure2_discordant_vs_rest") {
  w <- 12; h <- 5
} else {
  w <- 7;  h <- 6
}
    
    if (name == "SuperFig_all_key_supp") {
      w <- 20; h <- 12   # try 12–14
    } else if (name == "SuperFig_discordant_summary") {
      w <- 12; h <- 14
    } else if (name == "Figure2_discordant_vs_rest") {
      w <- 12; h <- 5
    } else {
      w <- 7;  h <- 6
    }
    
    
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(out, plot = p, device = ragg::agg_png,
                      width = w, height = h, units = "in", dpi = 300, bg = "white")
    } else {
      ggplot2::ggsave(out, plot = p, device = "png",
                      width = w, height = h, units = "in", dpi = 300, bg = "white")
    }
    
    msg("[FIG] wrote %s", out)
  }
  

  if (mode == "one") {
    stopifnot(!is.null(fig))
    build_one(fig)
  } else {
    for (nm in names(figs)) build_one(nm)
  }

  invisible(TRUE)
}

run_fig_build <- function(mode=c("all","plots-only","one"), fig=NULL, force=FALSE) {
  build_all_figures(mode=mode, fig=fig, force=force)
}

# ----------------------------
# CLI
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (!interactive() && length(args) > 0) {
  mode <- "all"; force <- FALSE; fig <- NULL
  for (a in args) {
    if (a %in% c("--all")) mode <- "all"
    if (a %in% c("--plots-only")) mode <- "plots-only"
    if (a %in% c("--force")) force <- TRUE
    if (str_detect(a, "^--one")) mode <- "one"
    if (str_detect(a, "^--fig=")) fig <- str_replace(a, "^--fig=", "")
  }
  run_fig_build(mode=mode, fig=fig, force=force)
}
