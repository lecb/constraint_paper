#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
})

API <- "https://gnomad.broadinstitute.org/api"

#KEEP FOR PAPER

# gnomAD v4 is GRCh38
REF     <- "GRCh38"
DATASET <- "gnomad_r4"

`%||%` <- function(a, b) if (is.null(a)) b else a

# -------------------------
# Low-level GraphQL caller
# -------------------------
gql_post <- function(query, variables = list()) {
  res <- httr::POST(
    API,
    encode = "json",
    body = list(query = query, variables = variables),
    httr::accept_json(),
    httr::user_agent("constraint-paper/1.0 (ellie)")
  )
  
  gql_post <- function(query, variables = list(), max_tries = 6, base_sleep = 1) {
    for (i in seq_len(max_tries)) {
      res <- httr::POST(
        API,
        encode = "json",
        body = list(query = query, variables = variables),
        httr::accept_json(),
        httr::user_agent("constraint-paper/1.0 (ellie)")
      )
      
      txt <- httr::content(res, "text", encoding = "UTF-8")
      status <- httr::status_code(res)
      js <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
      
      # If rate-limited, exponential backoff (+ jitter)
      if (status == 429) {
        sleep <- base_sleep * (2^(i - 1)) + runif(1, 0, 0.5)
        message("[429] Rate limited. Sleeping ", round(sleep, 2), "s then retry ", i, "/", max_tries)
        Sys.sleep(sleep)
        next
      }
      
      if (status != 200) {
        return(list(
          ok = FALSE, status = status, text = txt,
          data = NULL,
          errors = if (!is.null(js$errors)) js$errors else list(list(message = paste0("HTTP ", status)))
        ))
      }
      
      if (is.null(js)) {
        return(list(ok = FALSE, status = status, text = txt, data = NULL,
                    errors = list(list(message = "JSON parse failed"))))
      }
      
      if (!is.null(js$errors)) {
        return(list(ok = FALSE, status = status, text = txt, data = js$data %||% NULL, errors = js$errors))
      }
      
      return(list(ok = TRUE, status = status, text = txt, data = js$data, errors = NULL))
    }
    
    # Out of retries
    list(ok = FALSE, status = 429, text = "", data = NULL,
         errors = list(list(message = "HTTP 429 (rate limited): max retries exceeded")))
  }
  
  run_set <- function(syms, label) {
    message("[RUN] ", label, " n=", length(syms))
    res <- map(syms, function(s) {
      Sys.sleep(runif(1, 0.15, 0.35))  # gentle throttle
      fetch_gene_variants(s)
    })
    data  <- bind_rows(map(res, "data"))  %>% mutate(set = label)
    fails <- bind_rows(map(res, "fail"))  %>% mutate(set = label)
    list(data = data, fails = fails)
  }
  
  
  txt <- httr::content(res, "text", encoding = "UTF-8")
  status <- httr::status_code(res)
  
  js <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
  
  if (status != 200) {
    return(list(
      ok = FALSE, status = status, text = txt,
      data = NULL,
      errors = if (!is.null(js$errors)) js$errors else list(list(message = paste0("HTTP ", status)))
    ))
  }
  
  if (is.null(js)) {
    return(list(ok = FALSE, status = status, text = txt, data = NULL,
                errors = list(list(message = "JSON parse failed"))))
  }
  
  if (!is.null(js$errors)) {
    return(list(ok = FALSE, status = status, text = txt, data = js$data %||% NULL, errors = js$errors))
  }
  
  list(ok = TRUE, status = status, text = txt, data = js$data, errors = NULL)
}

err_to_msg <- function(out) {
  if (!is.null(out$errors) && length(out$errors) > 0) {
    m <- vapply(out$errors, function(e) as.character(e$message %||% NA_character_), character(1))
    m <- m[!is.na(m) & nzchar(m)]
    if (length(m) > 0) return(paste(m, collapse = " | "))
  }
  # fallback: include a snippet of the raw response
  snippet <- substr(as.character(out$text %||% ""), 1, 220)
  if (nzchar(snippet)) return(paste0("GraphQL failed (status=", out$status, "): ", snippet))
  paste0("GraphQL failed (status=", out$status, ")")
}

fail_row <- function(gene, stage, msg, extra = NA_character_) {
  tibble(
    gene  = as.character(gene),
    stage = as.character(stage),
    msg   = if (is.na(msg) || !nzchar(as.character(msg))) NA_character_ else as.character(msg),
    extra = if (is.na(extra) || !nzchar(as.character(extra))) NA_character_ else as.character(extra)
  )
}

# -------------------------
# Preferred approach (WORKS on gnomAD browser API):
# query by gene_symbol directly (no gene_search needed)
# -------------------------
Q_GENE_VARIANTS_BY_SYMBOL <- '
query GeneVariantsBySymbol($gene_symbol: String!, $reference_genome: ReferenceGenomeId!, $dataset: DatasetId!) {
  gene(gene_symbol: $gene_symbol, reference_genome: $reference_genome) {
    gene_id
    symbol
    variants(dataset: $dataset) {
      variant_id
      pos
      chrom
      ref
      alt
      hgvsc
      hgvsp
      consequence
      lof
      lof_flags
      flags
      exome { ac an af }
      genome { ac an af }
    }
  }
}
'



fetch_gene_variants <- function(symbol) {
  out <- gql_post(Q_GENE_VARIANTS_BY_SYMBOL,
                  list(gene_symbol = symbol, reference_genome = REF, dataset = DATASET))
  
  if (!out$ok) {
    return(list(data = tibble(), fail = fail_row(symbol, "gene_query", err_to_msg(out))))
  }
  
  g <- out$data$gene
  if (is.null(g)) {
    return(list(data = tibble(), fail = fail_row(symbol, "gene_query", "Gene not found (gene returned NULL)")))
  }
  
  vv <- g$variants
  if (is.null(vv) || length(vv) == 0) {
    # Not necessarily an error; still record so you can QC
    return(list(data = tibble(), fail = fail_row(symbol, "variants", "No variants returned for this gene")))
  }
  
  df <- tibble(
    gene_symbol = g$symbol %||% symbol,
    gene_id     = g$gene_id %||% NA_character_,
    variant_id  = vapply(vv, function(v) v$variant_id %||% NA_character_, character(1)),
    chrom       = vapply(vv, function(v) v$chrom %||% NA_character_, character(1)),
    pos         = suppressWarnings(as.integer(vapply(vv, function(v) v$pos %||% NA, numeric(1)))),
    ref         = vapply(vv, function(v) v$ref %||% NA_character_, character(1)),
    alt         = vapply(vv, function(v) v$alt %||% NA_character_, character(1)),
    consequence = vapply(vv, function(v) paste(v$consequence %||% character(0), collapse = ";"), character(1)),
    hgvsc       = vapply(vv, function(v) v$hgvsc %||% NA_character_, character(1)),
    hgvsp       = vapply(vv, function(v) v$hgvsp %||% NA_character_, character(1)),
    lof         = vapply(vv, function(v) v$lof %||% NA_character_, character(1)),
    lof_flags   = vapply(vv, function(v) paste(v$lof_flags %||% character(0), collapse = ";"), character(1)),
    flags       = vapply(vv, function(v) paste(v$flags %||% character(0), collapse = ";"), character(1)),
    exome_ac    = suppressWarnings(as.integer(vapply(vv, function(v) v$exome$ac %||% NA, numeric(1)))),
    exome_an    = suppressWarnings(as.integer(vapply(vv, function(v) v$exome$an %||% NA, numeric(1)))),
    exome_af    = suppressWarnings(as.numeric(vapply(vv, function(v) v$exome$af %||% NA, numeric(1)))),
    genome_ac   = suppressWarnings(as.integer(vapply(vv, function(v) v$genome$ac %||% NA, numeric(1)))),
    genome_an   = suppressWarnings(as.integer(vapply(vv, function(v) v$genome$an %||% NA, numeric(1)))),
    genome_af   = suppressWarnings(as.numeric(vapply(vv, function(v) v$genome$af %||% NA, numeric(1))))
  )
  
  list(data = df, fail = tibble())
}

# -------------------------
# Driver: read sets from RDS
# -------------------------
if (file.exists("discordant_genes.rds")) discordant_genes <- readRDS("discordant_genes.rds")
if (file.exists("non_discordant_genes.rds")) non_discordant_genes <- readRDS("non_discordant_genes.rds")

stopifnot(exists("discordant_genes"))
stopifnot("Gene" %in% names(discordant_genes))

disc_syms <- unique(as.character(discordant_genes$Gene))
non_syms  <- if (exists("non_discordant_genes") && "Gene" %in% names(non_discordant_genes)) {
  unique(as.character(non_discordant_genes$Gene))
} else {
  character(0)
}

run_set <- function(syms, label) {
  message("[RUN] ", label, " n=", length(syms))
  res <- map(syms, fetch_gene_variants)
  data  <- bind_rows(map(res, "data"))  %>% mutate(set = label)
  fails <- bind_rows(map(res, "fail"))  %>% mutate(set = label)
  list(data = data, fails = fails)
}

disc_res <- run_set(disc_syms, "Discordant")
non_res  <- run_set(non_syms,  "Non-discordant")

all_df  <- bind_rows(disc_res$data, non_res$data)
fail_df <- bind_rows(disc_res$fails, non_res$fails)

write_csv(all_df,  "gnomad_variants_all.csv", na = "")
write_csv(fail_df, "gnomad_variants_failures.csv", na = "")

disc35 <- unique(as.character(discordant_genes$Gene))

v <- readr::read_csv("gnomad_variants_all.csv", show_col_types = FALSE) %>%
  mutate(
    gene_u = toupper(trimws(gene_symbol)),
    consequence_l = tolower(consequence %||% ""),
    is_stop = str_detect(consequence_l, "stop_gained"),
    is_fs   = str_detect(consequence_l, "frameshift_variant"),
    is_trunc = is_stop | is_fs,
    # pLoF from gnomAD LOF annotation if present
    is_plof = tolower(lof %||% "") %in% c("hc", "lc") | is_trunc,
    flags_l = tolower(flags %||% ""),
    in_segdup = str_detect(flags_l, "segdup"),
    in_lcr    = str_detect(flags_l, "\\blcr\\b"),
    in_rep    = in_segdup | in_lcr,
    # totals across exome+genome
    AC = coalesce(exome_ac, 0L) + coalesce(genome_ac, 0L),
    hemi = coalesce(exome_hemizygote_count, genome_hemizygote_count, NA_integer_)
  ) %>%
  filter(gene_u %in% toupper(disc35))

flagaware_summary <- function(df, type = c("plof","trunc")) {
  type <- match.arg(type)
  df1 <- df %>%
    filter(if (type == "plof") is_plof else is_trunc)
  
  df_flagaware <- df1 %>% filter(!in_rep)
  
  df_flagaware %>%
    group_by(Gene = gene_symbol) %>%
    summarise(
      n_rows_type = n(),
      n_rows_removed_by_flags = sum(in_rep, na.rm = TRUE), # will be 0 here; kept for compatibility
      n_variants_flagaware = n_distinct(variant_id),
      AC_sum_flagaware = sum(AC, na.rm = TRUE),
      hom_sum_flagaware = if (all(is.na(hom))) NA_real_ else sum(hom, na.rm = TRUE),
      hemi_sum_flagaware = if (all(is.na(hemi))) NA_real_ else sum(hemi, na.rm = TRUE),
      classes_flagaware = paste(
        unique(c(if (any(is_stop, na.rm=TRUE)) "stop" else NULL,
                 if (any(is_fs,   na.rm=TRUE)) "frameshift" else NULL)),
        collapse = ", "
      ),
      .groups = "drop"
    )
}

plof_sum35  <- flagaware_summary(v, "plof")
trunc_sum35 <- flagaware_summary(v, "trunc")

write_csv(plof_sum35,  "flagaware_plof_summary_discordant35.csv",  na = "")
write_csv(trunc_sum35, "flagaware_trunc_summary_discordant35.csv", na = "")

# Build your combined supp table (same structure as before)
plof2 <- plof_sum35 %>%
  transmute(
    Gene,
    plof_n_rows_type              = n_rows_type,
    plof_n_rows_removed_by_flags  = n_rows_removed_by_flags,
    plof_n_variants_flagaware     = n_variants_flagaware,
    plof_AC_sum_flagaware         = AC_sum_flagaware,
    plof_hom_sum_flagaware        = hom_sum_flagaware,
    plof_hemi_sum_flagaware       = hemi_sum_flagaware,
    plof_classes_flagaware        = classes_flagaware
  )

trunc2 <- trunc_sum35 %>%
  transmute(
    Gene,
    trunc_n_rows_type             = n_rows_type,
    trunc_n_rows_removed_by_flags = n_rows_removed_by_flags,
    trunc_n_variants_flagaware    = n_variants_flagaware,
    trunc_AC_sum_flagaware        = AC_sum_flagaware,
    trunc_hom_sum_flagaware       = hom_sum_flagaware,
    trunc_hemi_sum_flagaware      = hemi_sum_flagaware,
    trunc_classes_flagaware       = classes_flagaware
  )

supp35 <- plof2 %>%
  inner_join(trunc2, by = "Gene") %>%
  mutate(
    plof_n_rows_retained_unflagged  = plof_n_rows_type  - plof_n_rows_removed_by_flags,
    trunc_n_rows_retained_unflagged = trunc_n_rows_type - trunc_n_rows_removed_by_flags,
    trunc_retains_ge1 = trunc_n_rows_retained_unflagged >= 1,
    trunc_retains_ge2 = trunc_n_rows_retained_unflagged >= 2
  ) %>%
  arrange(desc(trunc_n_rows_retained_unflagged), Gene)

dir.create("./tables", showWarnings = FALSE, recursive = TRUE)
write_csv(supp35, "./tables/SuppTable_flagaware_variant_filtering_discordant35.csv", na = "")


message("[DONE] wrote gnomad_variants_all.csv (rows=", nrow(all_df), ")")
message("[DONE] wrote gnomad_variants_failures.csv (rows=", nrow(fail_df), ")")

if (nrow(fail_df) > 0) {
  message("[FAIL SUMMARY]")
  print(fail_df %>% count(stage, msg, sort = TRUE), n = 50)
}

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

plof  <- read_csv("flagaware_plof_summary.csv",  show_col_types = FALSE)
trunc <- read_csv("flagaware_trunc_summary.csv", show_col_types = FALSE)

library(readr)
library(dplyr)
library(tidyr)
library(stringr)

plof  <- read_csv("flagaware_plof_summary.csv",  show_col_types = FALSE)
trunc <- read_csv("flagaware_trunc_summary.csv", show_col_types = FALSE)

# ---- NEW: define the discordant gene set (all 35) from the RDS you already loaded ----
# (Assumes discordant_genes.rds contains the full discordant set with column "Gene")
disc35 <- unique(as.character(discordant_genes$Gene))

# ---- helper to build the supp table for an arbitrary gene set ----
build_flagaware_supp <- function(plof, trunc, genes, label = "discordant") {
  plof_g <- plof %>% filter(Gene %in% genes)
  trunc_g <- trunc %>% filter(Gene %in% genes)
  
  # Sanity: ensure same genes in both inputs for this subset
  if (!setequal(plof_g$Gene, trunc_g$Gene)) {
    missing_in_plof  <- setdiff(trunc_g$Gene, plof_g$Gene)
    missing_in_trunc <- setdiff(plof_g$Gene, trunc_g$Gene)
    stop(
      "Gene mismatch between plof and trunc summaries for ", label, ".\n",
      "Missing in PLoF summary: ", paste(missing_in_plof, collapse = ", "), "\n",
      "Missing in trunc summary: ", paste(missing_in_trunc, collapse = ", ")
    )
  }
  
  plof2 <- plof_g %>%
    transmute(
      Gene,
      plof_n_rows_type              = n_rows_type,
      plof_n_rows_removed_by_flags  = n_rows_removed_by_flags,
      plof_n_variants_flagaware     = n_variants_flagaware,
      plof_AC_sum_flagaware         = AC_sum_flagaware,
      plof_hom_sum_flagaware        = hom_sum_flagaware,
      plof_hemi_sum_flagaware       = hemi_sum_flagaware,
      plof_classes_flagaware        = classes_flagaware
    )
  
  trunc2 <- trunc_g %>%
    transmute(
      Gene,
      trunc_n_rows_type             = n_rows_type,
      trunc_n_rows_removed_by_flags = n_rows_removed_by_flags,
      trunc_n_variants_flagaware    = n_variants_flagaware,
      trunc_AC_sum_flagaware        = AC_sum_flagaware,
      trunc_hom_sum_flagaware       = hom_sum_flagaware,
      trunc_hemi_sum_flagaware      = hemi_sum_flagaware,
      trunc_classes_flagaware       = classes_flagaware
    )
  
  supp <- plof2 %>%
    inner_join(trunc2, by = "Gene") %>%
    mutate(
      plof_n_rows_retained_unflagged  = plof_n_rows_type  - plof_n_rows_removed_by_flags,
      trunc_n_rows_retained_unflagged = trunc_n_rows_type - trunc_n_rows_removed_by_flags,
      trunc_retains_ge1 = trunc_n_rows_retained_unflagged >= 1,
      trunc_retains_ge2 = trunc_n_rows_retained_unflagged >= 2
    ) %>%
    arrange(desc(trunc_n_rows_retained_unflagged), Gene)
  
  # summary stats used in text
  summary_claims <- list(
    label = label,
    n_genes = nrow(supp),
    removed_plof_total  = sum(supp$plof_n_rows_removed_by_flags, na.rm = TRUE),
    removed_trunc_total = sum(supp$trunc_n_rows_removed_by_flags, na.rm = TRUE),
    n_retained_ge1 = sum(supp$trunc_retains_ge1, na.rm = TRUE),
    n_retained_ge2 = sum(supp$trunc_retains_ge2, na.rm = TRUE),
    genes_with_single_retained = paste(supp$Gene[supp$trunc_n_rows_retained_unflagged == 1], collapse = ", ")
  )
  
  list(supp = supp, summary = summary_claims)
}

# ---- Build and write TOP15 (backwards compatible) ----
# If your plof/trunc summaries are *already* top15-only, this will just reproduce the same output.
top15 <- disc35[1:min(15, length(disc35))]
res15 <- build_flagaware_supp(plof, trunc, top15, label = "top15")
dir.create("./tables", showWarnings = FALSE, recursive = TRUE)
write_csv(res15$supp, "./tables/SuppTable_flagaware_variant_filtering_top15.csv", na = "")

# ---- Build and write DISCORDANT35 ----
res35 <- build_flagaware_supp(plof, trunc, disc35, label = "discordant35")
write_csv(res35$supp, "./tables/SuppTable_flagaware_variant_filtering_discordant35.csv", na = "")

# Print summaries to console (useful for manuscript sentences)
res15$summary
res35$summary
