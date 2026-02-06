#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(httr)
  library(jsonlite)
  library(tibble)
  library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

INFILE  <- "gnomad_variants_all.csv"
OUT_CSV <- "nmd_heuristic_truncations.csv"
OUT_TXT <- "nmd_heuristic_results.txt"

stopifnot(file.exists(INFILE))
#KEEP FOR PAPER

# -------------------------
# Load + keep truncating
# -------------------------
v <- readr::read_csv(INFILE, show_col_types = FALSE)

if (!all(c("gene_symbol","set","chrom","pos","ref","alt","consequence") %in% names(v))) {
  stop("gnomad_variants_all.csv is missing required columns. Need: gene_symbol,set,chrom,pos,ref,alt,consequence")
}

x <- v %>%
  mutate(
    consequence_l = tolower(consequence %||% ""),
    is_trunc = str_detect(consequence_l, "stop_gained|frameshift_variant")
  ) %>%
  filter(is_trunc) %>%
  mutate(
    group = case_when(
      set == "Discordant" ~ "Discordant",
      set == "Non-discordant" ~ "Non-discordant",
      TRUE ~ as.character(set)
    )
  ) %>%
  filter(group %in% c("Discordant","Non-discordant")) %>%
  mutate(
    chrom = as.character(chrom),
    chrom = str_replace(chrom, "^chr", ""),
    pos = as.integer(pos)
  ) %>%
  filter(!is.na(pos), pos > 0, nchar(ref) >= 1, nchar(alt) >= 1)

if (nrow(x) == 0) stop("No truncating variants found after filtering.")

# -------------------------
# Ensembl GRCh38 REST helpers
# -------------------------
ENSEMBL_BASE <- "https://rest.ensembl.org"
UA <- "constraint-paper/1.0 (ellie)"

request_json <- function(url, max_tries = 4, backoff = 0.7) {
  for (i in seq_len(max_tries)) {
    res <- tryCatch(
      httr::GET(url, httr::accept("application/json"), httr::user_agent(UA), httr::timeout(30)),
      error = function(e) NULL
    )
    if (is.null(res)) {
      Sys.sleep(backoff * (2^(i-1))); next
    }
    code <- httr::status_code(res)
    txt  <- httr::content(res, "text", encoding = "UTF-8")
    
    if (code == 200) {
      return(tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL))
    }
    
    # transient
    if (code %in% c(408, 429, 500, 502, 503, 504)) {
      Sys.sleep(backoff * (2^(i-1)) + runif(1, 0, 0.4))
      next
    }
    
    # non-transient
    return(structure(list(), class = "ensembl_null", http_code = code, text = txt))
  }
  
  structure(list(), class = "ensembl_null", http_code = NA_integer_, text = "")
}

# canonical transcript (protein-coding) for gene symbol
get_canonical_tx <- function(gene_symbol) {
  g <- toupper(trimws(gene_symbol))
  url <- paste0(ENSEMBL_BASE, "/lookup/symbol/homo_sapiens/", g, "?expand=1")
  js <- request_json(url)
  if (is.null(js) || inherits(js, "ensembl_null")) return(NULL)
  txs <- js$Transcript
  if (is.null(txs) || length(txs) == 0) return(NULL)
  
  # choose canonical if flagged; else longest translation
  canon_idx <- which(vapply(txs, function(t) isTRUE(t$is_canonical %||% FALSE), logical(1)))
  if (length(canon_idx) > 0) {
    tx <- txs[[canon_idx[1]]]
  } else {
    lens <- vapply(txs, function(t) {
      if (!is.null(t$Translation$length)) suppressWarnings(as.numeric(t$Translation$length)) else NA_real_
    }, numeric(1))
    if (all(is.na(lens))) return(NULL)
    tx <- txs[[which.max(lens)]]
  }
  
  if (is.null(tx$Exon) || length(tx$Exon) == 0) return(NULL)
  if (is.null(tx$Translation$start) || is.null(tx$Translation$end)) return(NULL)
  
  list(
    gene = g,
    tx_id = tx$id %||% NA_character_,
    strand = as.integer(tx$strand %||% 1L),
    tr_start = as.integer(tx$Translation$start),
    tr_end   = as.integer(tx$Translation$end),
    exons = tx$Exon
  )
}

# compute: last coding exon genomic interval + last junction genomic position
compute_nmd_heuristic_meta <- function(tx) {
  if (is.null(tx)) return(NULL)
  
  tr_lo <- min(tx$tr_start, tx$tr_end)
  tr_hi <- max(tx$tr_start, tx$tr_end)
  
  ex_tbl <- tibble(
    start = as.integer(vapply(tx$exons, function(e) e$start %||% NA, numeric(1))),
    end   = as.integer(vapply(tx$exons, function(e) e$end   %||% NA, numeric(1)))
  ) %>% filter(!is.na(start), !is.na(end))
  
  if (nrow(ex_tbl) == 0) return(NULL)
  
  # order exons in transcription order
  if (tx$strand == 1) {
    ex_tbl <- ex_tbl %>% arrange(start, end)
  } else {
    ex_tbl <- ex_tbl %>% arrange(desc(start), desc(end))
  }
  
  # mark coding overlap with CDS span (translation start/end are genomic coords in transcript)
  overlap_len <- function(a1,a2,b1,b2){
    lo <- max(a1,b1); hi <- min(a2,b2)
    if (hi < lo) 0L else as.integer(hi - lo + 1L)
  }
  ex_tbl$coding_bp <- mapply(function(s,e) overlap_len(s,e,tr_lo,tr_hi), ex_tbl$start, ex_tbl$end)
  
  coding_ex <- ex_tbl %>% filter(coding_bp > 0)
  if (nrow(coding_ex) == 0) return(NULL)
  
  last_coding_exon <- coding_ex[nrow(coding_ex), , drop = FALSE]
  penult_coding_exon <- if (nrow(coding_ex) >= 2) coding_ex[nrow(coding_ex)-1, , drop = FALSE] else NULL
  
  # last exon-exon junction (between penultimate and last coding exon) in genomic coords:
  # For + strand: junction at end of penultimate exon; for - strand: junction at start of penultimate exon
  last_junction_genomic <- NA_integer_
  if (!is.null(penult_coding_exon)) {
    last_junction_genomic <- if (tx$strand == 1) {
      as.integer(penult_coding_exon$end)
    } else {
      as.integer(penult_coding_exon$start)
    }
  }
  
  list(
    strand = tx$strand,
    last_exon_start = as.integer(last_coding_exon$start),
    last_exon_end   = as.integer(last_coding_exon$end),
    last_junction_genomic = last_junction_genomic
  )
}

# per gene meta cache
meta_cache <- new.env(parent = emptyenv())

get_gene_meta <- function(gene) {
  g <- toupper(trimws(gene))
  if (exists(g, envir = meta_cache, inherits = FALSE)) return(get(g, envir = meta_cache))
  tx <- get_canonical_tx(g)
  meta <- compute_nmd_heuristic_meta(tx)
  assign(g, meta, envir = meta_cache)
  meta
}

# flag NMD escape:
# - in last coding exon (genomic interval)
# - OR within 50 bp of last junction (towards the 3' end of the penultimate exon)
is_nmd_escape <- function(pos, meta) {
  if (is.null(meta)) return(NA)
  pos <- as.integer(pos)
  if (is.na(pos)) return(NA)
  
  in_last_exon <- pos >= meta$last_exon_start && pos <= meta$last_exon_end
  
  near_last_junction <- NA
  if (!is.na(meta$last_junction_genomic)) {
    near_last_junction <- if (meta$strand == 1) {
      # within 50 bp upstream of junction => [junction-50, junction]
      pos >= (meta$last_junction_genomic - 50L) && pos <= meta$last_junction_genomic
    } else {
      # on - strand, upstream of junction is increasing coords away from start
      pos >= meta$last_junction_genomic && pos <= (meta$last_junction_genomic + 50L)
    }
  }
  
  isTRUE(in_last_exon) || isTRUE(near_last_junction)
}

# -------------------------
# Apply gene meta + compute escape flags
# -------------------------
message("[RUN] computing NMD heuristic per trunc variant... (genes=", n_distinct(x$gene_symbol), ")")

# compute once per gene
genes <- sort(unique(toupper(trimws(x$gene_symbol))))
invisible(lapply(genes, get_gene_meta))

x2 <- x %>%
  mutate(
    gene_u = toupper(trimws(gene_symbol)),
    nmd_escape = pmap_lgl(list(pos, gene_u), function(p, g) {
      meta <- get(g, envir = meta_cache, inherits = FALSE)
      val <- is_nmd_escape(p, meta)
      if (is.na(val)) FALSE else val
    })
  )

write_csv(x2, OUT_CSV, na = "")
message("[DONE] wrote ", OUT_CSV, " (rows=", nrow(x2), ")")

# -------------------------
# Compare proportions (Discordant vs Non-discordant)
# -------------------------
tab <- with(x2, table(group, nmd_escape))
tab <- tab[c("Discordant","Non-discordant"), , drop = FALSE]

ft <- NULL
if (all(dim(tab) == c(2,2))) {
  ft <- fisher.test(tab)
}

summ <- x2 %>%
  group_by(group) %>%
  summarise(
    n_trunc = n(),
    n_escape = sum(nmd_escape, na.rm = TRUE),
    prop_escape = n_escape / pmax(1, n_trunc),
    .groups = "drop"
  )

lines <- c(
  "---- NMD heuristic: last coding exon OR within 50bp of last junction ----",
  paste0("Input: ", INFILE),
  paste0("Output variants: ", OUT_CSV),
  "",
  "Counts by group:",
  capture.output(print(summ)),
  "",
  "2x2 table (group x nmd_escape):",
  capture.output(print(tab))
)

if (!is.null(ft)) {
  lines <- c(lines, "",
             paste0("Fisher's exact test p = ", signif(ft$p.value, 4)),
             paste0("OR = ", signif(unname(ft$estimate), 4),
                    " (95% CI ", signif(ft$conf.int[1], 4), "–", signif(ft$conf.int[2], 4), ")")
  )
} else {
  lines <- c(lines, "", "[SKIP] Fisher test not run (table not 2x2; do you have both groups?)")
}

writeLines(lines, OUT_TXT)
message("[DONE] wrote ", OUT_TXT)

cat(readLines("nmd_heuristic_results.txt"), sep = "\n")

library(dplyr)
library(ggplot2)
library(readr)
library(scales)

# x2 assumed to be the output from 04_nmd_heuristic_compare.R
# containing: gene_symbol, group, nmd_escape (TRUE/FALSE)

gene_summ <- x2 %>%
  group_by(group, gene_symbol) %>%
  summarise(
    prop_escape = mean(nmd_escape),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3)   # optional but recommended

p <- ggplot(gene_summ, aes(x = group, y = prop_escape, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(
    "Discordant"     = "firebrick3",
    "Non-discordant" = "steelblue4"
  )) +
  labs(
    y = "Per-gene proportion of truncations predicted to escape NMD",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

ggsave(
  filename = "Supplementary_Figure_NMD_escape_per_gene.png",
  plot = p,
  width = 6.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  filename = "Supplementary_Figure_NMD_escape_per_gene.pdf",
  plot = p,
  width = 6.5,
  height = 4.5
)

message("[DONE] Wrote Supplementary_Figure_NMD_escape_per_gene.(png/pdf)")


