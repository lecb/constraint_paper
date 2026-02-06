#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(httr)
  library(jsonlite)
})



# ============================================================
# 01_generate_relpos.R (CLEAN, SELF-CONTAINED)
#
# Goal:
#   Build trunc_position_domain_nmd_variants.csv for your LOEUF<0.2
#   precision-filtered universe, with BOTH groups:
#     - Background (non-discordant within universe)
#     - Discordant (disc35 / Supplementary_Table1 list)
#
# Key constraints:
#   Your gnomAD LoF bgz lacks protein position / hgvsp, so we MUST use VEP
#   to obtain protein_id + protein_start, then Ensembl lookup for protein length.
#
# Outputs:
#   - trunc_position_domain_nmd_variants.csv
#   - trunc_position_domain_nmd_summary_by_group.csv
#
# Caching:
#   - cache/trunc_pdn_cache.rds
#   - cache/vep_relpos_checkpoint.rds
#   - cache/ensembl_translation_length_cache.rds
# ============================================================

# ----------------------------
# CLI FLAGS
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
HAS_FLAG <- function(x) any(args == x)

FORCE_RECOMPUTE <- HAS_FLAG("--force-recompute")

# Default: if you source in RStudio, behave like a normal run
if (interactive() && length(args) == 0) {
  FORCE_RECOMPUTE <- FALSE
}

message("[MODE] FORCE_RECOMPUTE=", FORCE_RECOMPUTE)

# ----------------------------
# INPUTS
# ----------------------------
UNIVERSE_FILE <- {
  p1 <- file.path("tables","loeuf_lt_0.2_precision_filtered.csv");
  p2 <- file.path("gnomad_lof_discordance_out","loeuf_lt_0.2_precision_filtered.csv");
  if (file.exists(p1)) p1 else p2
}
DISC_FILE     <- file.path("tables", "Supplementary_Table1_discordant_genes.csv")

stopifnot(file.exists(UNIVERSE_FILE))
stopifnot(file.exists(DISC_FILE))

# gnomAD LoF bgz (already on disk in your run logs)

GNOMAD_LOF_BGZ <- {
  candidates <- c(
    file.path(Sys.getenv("PIPELINE_ROOT", ""), "data", "gnomad", "gnomad.v2.1.1.all_lofs.txt.bgz"),
    file.path("data", "gnomad", "gnomad.v2.1.1.all_lofs.txt.bgz"),
    file.path("..", "gnomad.v2.1.1.all_lofs.txt.bgz")
  )

  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("Could not find gnomad.v2.1.1.all_lofs.txt.bgz. Tried: ", paste(candidates, collapse=" | "))
  invisible(hit)}
stopifnot(file.exists(GNOMAD_LOF_BGZ))


message("[DEBUG] GNOMAD_LOF_BGZ resolved to: ", GNOMAD_LOF_BGZ)
message("[DEBUG] exists? ", file.exists(GNOMAD_LOF_BGZ))
message("[DEBUG] readable? ", file.access(GNOMAD_LOF_BGZ, 4) == 0)
# ----------------------------
# OUTPUTS + CACHE DIR
# ----------------------------
OUT_TRUNC_PDN_VARIANTS <- "trunc_position_domain_nmd_variants.csv"
OUT_TRUNC_PDN_SUMMARY  <- "trunc_position_domain_nmd_summary_by_group.csv"

CACHEDIR <- "cache"
dir.create(CACHEDIR, showWarnings = FALSE, recursive = TRUE)

TRUNC_PDN_RDS        <- file.path(CACHEDIR, "trunc_pdn_cache.rds")
VEP_CHECKPOINT_RDS   <- file.path(CACHEDIR, "vep_relpos_checkpoint.rds")
TXLEN_CACHE_RDS      <- file.path(CACHEDIR, "ensembl_translation_length_cache.rds")

# ----------------------------
# LOAD GENE SETS
# ----------------------------
universe_genes <- readr::read_csv(UNIVERSE_FILE, show_col_types = FALSE) %>%
  pull(gene) %>% toupper() %>% unique()

disc35 <- readr::read_csv(DISC_FILE, show_col_types = FALSE) %>%
  pull(Gene) %>% toupper() %>% unique()

message("[relpos] universe genes: ", length(universe_genes))
message("[relpos] discordant genes: ", length(disc35))

# ----------------------------
# QUICK RETURN FROM CACHE
# ----------------------------
if (!FORCE_RECOMPUTE && file.exists(TRUNC_PDN_RDS)) {
  message("[CACHE] Using cached trunc_pdn: ", TRUNC_PDN_RDS)
  trunc_pdn <- readRDS(TRUNC_PDN_RDS)
  
  # write outputs (guarantee they exist)
  readr::write_csv(trunc_pdn, OUT_TRUNC_PDN_VARIANTS)
  message("Wrote: ", OUT_TRUNC_PDN_VARIANTS)
  
  summ <- trunc_pdn %>%
    group_by(group) %>%
    summarise(
      n_variants_total = n(),
      n_with_relpos = sum(!is.na(rel_pos)),
      relpos_median = median(rel_pos, na.rm = TRUE),
      relpos_IQR = IQR(rel_pos, na.rm = TRUE),
      early = sum(pos_bin == "early", na.rm = TRUE),
      middle = sum(pos_bin == "middle", na.rm = TRUE),
      late = sum(pos_bin == "late", na.rm = TRUE),
      prop_late = late / pmax(1, early + middle + late),
      .groups = "drop"
    )
  readr::write_csv(summ, OUT_TRUNC_PDN_SUMMARY)
  message("Wrote: ", OUT_TRUNC_PDN_SUMMARY)
  
  message("[DONE] (cache)")
invisible(NULL)
  invisible(NULL)
}

# ============================================================
# HELPERS
# ============================================================
`%||%` <- function(a, b) if (is.null(a)) b else a

ensembl_null <- function() structure(list(), class = "ensembl_null")

bin_relpos <- function(r) {
  dplyr::case_when(
    is.na(r) ~ NA_character_,
    r < 1/3 ~ "early",
    r < 2/3 ~ "middle",
    TRUE ~ "late"
  )
}

# Ensembl GRCh37 REST base (matches your VEP usage)
ENSEMBL_BASE <- "https://grch37.rest.ensembl.org"
ENSEMBL_UA   <- "constraint-paper/1.0 (ellie)"

.request_json <- function(
    verb = c("GET","POST"),
    url,
    body = NULL,
    ua = ENSEMBL_UA,
    connecttimeout = 15,
    timeout = 90,
    max_tries = 6,
    backoff_base = 0.8,
    quiet = FALSE,
    max_elapsed = 600
) {
  verb <- match.arg(verb)
  t0 <- Sys.time()
  
  for (attempt in seq_len(max_tries)) {
    
    if (as.numeric(difftime(Sys.time(), t0, units = "secs")) > max_elapsed) {
      if (!quiet) message("ABORT (elapsed>", max_elapsed, "s): ", url)
    }
    
    res <- tryCatch({
      cfg <- httr::config(
        connecttimeout = connecttimeout,
        low_speed_time = 30,
        low_speed_limit = 1
      )
      if (verb == "GET") {
        httr::GET(url, httr::accept("application/json"), httr::user_agent(ua),
                  httr::timeout(timeout), cfg)
      } else {
        httr::POST(url, httr::accept("application/json"), httr::content_type_json(),
                   httr::user_agent(ua), httr::timeout(timeout), cfg,
                   body = body, encode = "json")
      }
    }, error = function(e) NULL)
    
    code <- if (!is.null(res)) httr::status_code(res) else NA_integer_
    
    if (!is.null(res) && code == 200L) {
      txt <- httr::content(res, "text", encoding = "UTF-8")
      if (!nzchar(txt)) return(ensembl_null())
      js <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) NULL)
      if (is.null(js)) return(ensembl_null())
      return(js)
    }
    
    transient <- is.na(code) || code %in% c(408, 429, 500, 502, 503, 504)
    if (!transient) {
      if (!quiet) message("HTTP ", code, ": ", url)

    }
    
    if (!quiet) message("Retry ", attempt, "/", max_tries, " (HTTP ", code, "): ", url)
    Sys.sleep((backoff_base * (2^(attempt-1))) + runif(1, 0, 0.6))
  }
  
  if (!quiet) message("Failed after retries: ", url)
}

# VEP region endpoint (GRCh37)
vep_region_post <- function(variants) {
  url <- paste0(
    ENSEMBL_BASE,
    "/vep/homo_sapiens/region",
    "?protein=1",
    "&canonical=1",
    "&hgvs=1"
  )

  .request_json(
    "POST", url,
    body = list(variants = variants),
    timeout = 180, connecttimeout = 30,
    max_tries = 6, max_elapsed = 900,
    quiet = FALSE
  )
}

# Ensembl lookup for translation lengths (protein_id => length)
txlen_cache <- new.env(parent = emptyenv())
if (file.exists(TXLEN_CACHE_RDS)) {
  tmp <- readRDS(TXLEN_CACHE_RDS)
  try(list2env(tmp, envir = txlen_cache), silent = TRUE)
  message("[CACHE] Loaded translation lengths: ", TXLEN_CACHE_RDS)
}

lookup_lengths_bulk <- function(ids, cache_env, batch = 1500L) {
  ids <- unique(sub("\\..*$", "", as.character(ids)))
  ids <- ids[nchar(ids) > 0]
  ids <- ids[grepl("^ENSP[0-9]+$", ids)]
  ids <- ids[!vapply(ids, function(x) exists(x, envir = cache_env, inherits = FALSE), logical(1))]
  if (length(ids) == 0) return(invisible(NULL))

  
  url <- paste0(ENSEMBL_BASE, "/lookup/id")
  for (i in seq(1, length(ids), by = batch)) {
    chunk <- ids[i:min(length(ids), i + batch - 1)]
    js <- .request_json("POST", url, body = list(ids = chunk),
                        timeout = 120, connecttimeout = 20, max_tries = 6, max_elapsed = 600, quiet = FALSE)
    if (is.null(js) || inherits(js, "ensembl_null")) {
      for (id in chunk) assign(id, NA_real_, envir = cache_env)
      next
    }
    for (id in chunk) {
      obj <- js[[id]]
      len <- NA_real_
      if (!is.null(obj) && !is.null(obj$length)) len <- suppressWarnings(as.numeric(obj$length))
      if (!is.finite(len) || len <= 0) len <- NA_real_
      assign(id, len, envir = cache_env)
    }
  }
  invisible(NULL)
}

get_translation_length_cached <- function(ensp_id, cache_env) {
  # Always return length-1 numeric (NA_real_ on failure)
  if (is.null(ensp_id) || is.na(ensp_id) || !nzchar(as.character(ensp_id))) return(NA_real_)

  pid <- sub("\\..*$", "", as.character(ensp_id))

  # cache hit
  if (exists(pid, envir = cache_env, inherits = FALSE)) {
    v <- get(pid, envir = cache_env, inherits = FALSE)
    v <- suppressWarnings(as.numeric(v))
    if (length(v) != 1 || !is.finite(v) || v <= 0) return(NA_real_)
    return(v)
  }

  # Ensembl single-id lookup
  url <- paste0(ENSEMBL_BASE, "/lookup/id/", pid, "?expand=1")
  js <- .request_json("GET", url, quiet = TRUE)

  if (is.null(js) || inherits(js, "ensembl_null") || length(js) == 0) {
    assign(pid, NA_real_, envir = cache_env)
    return(NA_real_)
  }

  len <- NA_real_
  if (!is.null(js$length)) len <- suppressWarnings(as.numeric(js$length))
  if (length(len) != 1 || !is.finite(len) || len <= 0) len <- NA_real_

  assign(pid, len, envir = cache_env)
  return(len)
}

canonical_flag <- function(x) {
  # Always return length-1 logical; default FALSE
  if (is.null(x) || length(x) == 0 || is.na(x)[1]) return(FALSE)
  if (is.logical(x)) return(isTRUE(x[1]))

  x <- toupper(as.character(x[1]))
  x %in% c("1", "TRUE", "YES", "Y")
}

extract_canonical_hit <- function(rec) {
  tcs <- rec$transcript_consequences
  if (is.null(tcs) || length(tcs) == 0) return(NULL)

  ord <- order(vapply(tcs, function(tc) canonical_flag(tc$canonical), logical(1)), decreasing = TRUE)
  tcs <- tcs[ord]

  for (tc in tcs) {
    pid <- as.character(tc$protein_id %||% "")
    aa  <- suppressWarnings(as.numeric(tc$protein_start %||% NA_real_))

    if (nzchar(pid) && is.finite(aa)) {
      pid <- sub("\\..*$", "", pid)
      return(list(protein_id = pid, aa_pos = aa))
    }
  }
  NULL
}

# checkpoint helpers
recover_vep_checkpoint <- function(path) {
  if (!file.exists(path)) return(ensembl_null())

  ck <- readRDS(path)
  if (!is.list(ck) || is.null(ck$done) || is.null(ck$rel_tbl)) {

  }
  ck$rel_tbl <- tibble::as_tibble(ck$rel_tbl)
  if (!all(c("protein_id","aa_pos") %in% names(ck$rel_tbl))) {
    ck$rel_tbl <- tibble::tibble(protein_id=character(0), aa_pos=numeric(0))
  }
  ck
}

save_vep_checkpoint <- function(done, rel_tbl, path) {
  saveRDS(list(done = done, rel_tbl = rel_tbl), path)
}

# gnomAD row => VEP "region" variant string (chr pos id ref alt)
to_vep_variant <- function(chrom, pos, ref, alt, id = ".") {
  chrom <- sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  pos   <- suppressWarnings(as.integer(pos))
  ref   <- as.character(ref)
  alt   <- as.character(alt)
  
  if (!is.finite(pos) || pos <= 0) return(ensembl_null())
  if (!nzchar(chrom) || !nzchar(ref) || !nzchar(alt)) return(ensembl_null())
  
  paste(chrom, pos, id, ref, alt, ".", ".", ".", sep = " ")
}

# ============================================================
# 1) LOAD + FILTER GNOMAD LOF (to your universe genes + truncating)
# ============================================================
message("[bgz] Reading header + streaming filtered rows...")

# read header
hdr <- readLines(pipe(sprintf("zcat -f %s | head -n 1", shQuote(GNOMAD_LOF_BGZ))), n = 1)
stopifnot(length(hdr) == 1)
cn <- strsplit(hdr, "\t", fixed = TRUE)[[1]] %>% trimws()

# v2.1.1 all_lofs header includes:
# chrom pos ref alt most_severe_consequence gene_ids gene_symbols transcript_ids
stopifnot(all(c("chrom","pos","ref","alt","most_severe_consequence","gene_symbols") %in% cn))

# awk-filter to only universe genes (by gene_symbols)
tmp_genes <- tempfile(fileext = ".txt")
writeLines(sort(unique(universe_genes)), tmp_genes)

gene_col_idx <- match("gene_symbols", cn)

awk_cmd <- sprintf(
  "awk -F'\t' 'NR==FNR{a[$1]=1;next} FNR==1{print;next} ($%d in a){print}' %s <(zcat -f %s)",
  gene_col_idx,
  shQuote(tmp_genes),
  shQuote(GNOMAD_LOF_BGZ)
)
cmd <- sprintf("bash -lc %s", shQuote(awk_cmd))
gn <- readr::read_tsv(pipe(cmd), show_col_types = FALSE, progress = TRUE)

message("[bgz] Filtered rows (universe genes): ", nrow(gn))
if (nrow(gn) == 0) stop("No rows found for universe genes in gnomAD LoF table.")

gn <- gn %>%
  mutate(
    geneU = toupper(gene_symbols),
    cons  = tolower(most_severe_consequence)
  ) %>%
  filter(geneU %in% universe_genes) %>%
  filter(str_detect(cons, "stop_gained|frameshift_variant"))

message("[bgz] Rows after truncation filter: ", nrow(gn))
if (nrow(gn) == 0) stop("No truncating rows for universe genes after filter.")

# ============================================================
# 2) VEP to obtain protein_id + protein_start (aa_pos)
# ============================================================
message("[VEP] Building VEP variant strings...")
vep_vars <- vapply(seq_len(nrow(gn)), function(i) {
  to_vep_variant(gn$chrom[i], gn$pos[i], gn$ref[i], gn$alt[i])
}, character(1))

keep_ok <- !is.na(vep_vars) & nzchar(vep_vars)
gn <- gn[keep_ok, , drop = FALSE]
vep_vars <- vep_vars[keep_ok]
vep_vars <- unique(vep_vars)

message("[VEP] Unique variants to annotate: ", length(vep_vars))
if (length(vep_vars) == 0) stop("No usable variants to send to VEP.")

# checkpoint resume
ck <- recover_vep_checkpoint(VEP_CHECKPOINT_RDS)
done <- ck$done
rel_tbl <- ck$rel_tbl

if (length(done) != length(vep_vars)) {
  done <- rep(FALSE, length(vep_vars))
  rel_tbl <- tibble::tibble(protein_id=character(0), aa_pos=numeric(0))
}

message("[VEP] Resume: done ", sum(done), "/", length(done), " | hits n=", nrow(rel_tbl))

VEP_BATCH <- if (interactive()) 80L else 200L
VEP_THROTTLE_S <- if (interactive()) 0.12 else 0.25

for (i in seq(1, length(vep_vars), by = VEP_BATCH)) {
  idx <- i:min(length(vep_vars), i + VEP_BATCH - 1L)
  idx <- idx[!done[idx]]
  if (length(idx) == 0) next
  
  js <- vep_region_post(vep_vars[idx])
  if (is.null(js) || length(js) == 0) {
    message("[VEP] Empty response for batch starting ", i, " (marking done)")
    done[idx] <- TRUE
    save_vep_checkpoint(done, rel_tbl, VEP_CHECKPOINT_RDS)
    Sys.sleep(VEP_THROTTLE_S)
    next
  }
  
  picks <- list()
  for (rec in js) {
    hit <- extract_canonical_hit(rec)
    if (!is.null(hit)) picks[[length(picks) + 1L]] <- hit
  }
  
  if (length(picks) > 0) {
    rel_tbl <- bind_rows(
      rel_tbl,
      tibble::tibble(
        protein_id = vapply(picks, `[[`, character(1), "protein_id"),
        aa_pos      = vapply(picks, `[[`, numeric(1), "aa_pos")
      ))
  }
  
  done[idx] <- TRUE
  
  # checkpoint every ~500 variants
  if (sum(done) %% 500L == 0L) {
    save_vep_checkpoint(done, rel_tbl, VEP_CHECKPOINT_RDS)
    message("  [checkpoint] done=", sum(done), "/", length(done), " | hits n=", nrow(rel_tbl))
  }
  
  Sys.sleep(VEP_THROTTLE_S)
}

save_vep_checkpoint(done, rel_tbl, VEP_CHECKPOINT_RDS)
message("[VEP] Finished. Total protein_id+aa_pos hits: ", nrow(rel_tbl))

if (nrow(rel_tbl) == 0) stop("VEP returned no usable canonical protein hits; cannot compute rel_pos.")

# ============================================================
# 3) Translation lengths => rel_pos
# ============================================================
ids <- unique(rel_tbl$protein_id)
invisible(lookup_lengths_bulk(ids, cache_env = txlen_cache, batch = 1500L))
try(saveRDS(as.list(txlen_cache), TXLEN_CACHE_RDS), silent = TRUE)

prot_len <- vapply(rel_tbl$protein_id, function(pid) {
  get_translation_length_cached(pid, cache_env = txlen_cache)
}, numeric(1))

rel_pos <- rel_tbl$aa_pos / prot_len
ok <- is.finite(rel_pos) & rel_pos >= 0 & rel_pos <= 1 &
  is.finite(rel_tbl$aa_pos) & rel_tbl$aa_pos > 0 &
  is.finite(prot_len) & prot_len > 0

rel_tbl2 <- rel_tbl[ok, , drop = FALSE] %>%
  mutate(prot_len_aa = prot_len[ok],
         rel_pos = rel_pos[ok])

message("[relpos] Rows with finite rel_pos after VEP: ", nrow(rel_tbl2))
if (nrow(rel_tbl2) == 0) stop("No finite rel_pos could be computed after protein length lookup.")

# ============================================================
# 4) Build final trunc_pdn table (two groups)
#    NOTE: gnomAD LoF file itself does not retain gene per VEP hit easily.
#          So we build a positional distribution table for the universe using VEP hits.
#          Group assignment is done by mapping gene membership (disc35 vs rest).
#
#    Practically: we only need rel_pos distributions by group; gene label per variant is not essential.
#    But if you want Gene per row, you must annotate each VEP hit back to gene via VEP outputs.
#    Here we do the simpler and robust version: sample is "universe trunc variants",
#    and we stratify by gene set membership at the gene level is not possible without gene mapping.
#
#    HOWEVER: You *already* got 2292 rows and group counts printing earlier — that implies your
#    local run is treating the universe trunc set as "Background" plus "Discordant" by Gene.
#    To preserve that behavior, we do gene mapping by using the original filtered gnomAD rows
#    and joining on variant identity (chrom/pos/ref/alt) back from VEP inputs.
# ============================================================

# Reconstruct variant key from vep string: "chrom pos id ref alt ..."
parse_vep_key <- function(v) {
  parts <- strsplit(v, " ", fixed = TRUE)[[1]]
  if (length(parts) < 5) return(ensembl_null())
  tibble(
    chrom = parts[1],
    pos   = suppressWarnings(as.integer(parts[2])),
    ref   = parts[4],
    alt   = parts[5]
  )
}

# Build a table of (chrom,pos,ref,alt) in same order as VEP hits.
# We don’t have a direct mapping from each VEP hit back to a single row,
# so we join using (chrom,pos,ref,alt) from the original gn.
# We will rebuild keys from gn itself:
gn_key <- gn %>%
  transmute(
    chrom = as.character(chrom),
    pos   = suppressWarnings(as.integer(pos)),
    ref   = as.character(ref),
    alt   = as.character(alt),
    geneU = geneU
  ) %>%
  distinct()

# Create a pool of rel_pos rows by assigning each rel_pos to any matching key.
# We don’t have per-hit key stored, so we approximate by repeating rel_pos
# across all matching keys in gn_key is NOT OK.
#
# Therefore: we need to re-run VEP in a way that retains the input variant key per record.
# Ensembl VEP region returns 'input' field per record: e.g. "1 123 . A G . . ."
# We'll re-derive from checkpointed VEP by re-querying only the variants once here,
# but cheaply: we can do a single pass over vep_vars (already built) and request in batches
# and store mapping (input -> protein_id, aa_pos). That is still VEP, but we can cache it.
#
# To keep this script reliable and avoid more complexity:
# We'll build trunc_pdn with a *variant_id* key as "chrom-pos-ref-alt" and map gene via gn_key.
# ============================================================

MAP_RDS <- file.path(CACHEDIR, "vep_input_to_proteinpos_map.rds")

if (!FORCE_RECOMPUTE && file.exists(MAP_RDS)) {
  message("[CACHE] Using cached VEP input map: ", MAP_RDS)
  vep_map <- readRDS(MAP_RDS)
} else {
  message("[VEP MAP] Building input->(protein_id,aa_pos) map ...")
  vep_map <- tibble::tibble(
    chrom = character(0),
    pos   = integer(0),
    ref   = character(0),
    alt   = character(0),
    protein_id = character(0),
    aa_pos = numeric(0)
  )
  done2 <- rep(FALSE, length(vep_vars))
  VEP_BATCH2 <- if (interactive()) 80L else 200L
  
  for (i in seq(1, length(vep_vars), by = VEP_BATCH2)) {
    idx <- i:min(length(vep_vars), i + VEP_BATCH2 - 1L)
    js <- vep_region_post(vep_vars[idx])
    if (is.null(js) || length(js) == 0) {
      Sys.sleep(VEP_THROTTLE_S)
      next
    }
    
    rows <- list()
    for (rec in js) {
      inp <- rec$input %||% ""
      hit <- extract_canonical_hit(rec)
      if (!nzchar(inp) || is.null(hit)) next
      
      # parse the input back to key
      parts <- strsplit(inp, " ", fixed = TRUE)[[1]]
      if (length(parts) < 5) next
      rows[[length(rows) + 1L]] <- tibble(
        chrom = as.character(parts[1]),
        pos   = suppressWarnings(as.integer(parts[2])),
        ref   = as.character(parts[4]),
        alt   = as.character(parts[5]),
        protein_id = as.character(hit$protein_id),
        aa_pos = as.numeric(hit$aa_pos)
    )
    }
    
    if (length(rows) > 0) vep_map <- bind_rows(vep_map, bind_rows(rows))
    Sys.sleep(VEP_THROTTLE_S)
  }
  
  vep_map <- vep_map %>% filter(!is.na(pos), chrom != "", ref != "", alt != "", protein_id != "", is.finite(aa_pos))
  saveRDS(vep_map, MAP_RDS)
  message("[CACHE] Saved: ", MAP_RDS, " (rows=", nrow(vep_map), ")")
}

# protein lengths for vep_map
invisible(lookup_lengths_bulk(unique(vep_map$protein_id), cache_env = txlen_cache, batch = 1500L))
try(saveRDS(as.list(txlen_cache), TXLEN_CACHE_RDS), silent = TRUE)

vep_map <- vep_map %>%
  mutate(
    prot_len_aa = vapply(protein_id, function(pid) get_translation_length_cached(pid, txlen_cache), numeric(1)),
    rel_pos = aa_pos / prot_len_aa
  ) %>%
  filter(is.finite(rel_pos), rel_pos >= 0, rel_pos <= 1)

message("[VEP MAP] Rows with finite rel_pos: ", nrow(vep_map))
if (nrow(vep_map) == 0) stop("VEP mapping produced no finite rel_pos rows.")

# Join geneU from gnomAD filtered table
trunc_pdn <- vep_map %>%
  left_join(gn_key, by = c("chrom","pos","ref","alt")) %>%
  filter(!is.na(geneU), geneU %in% universe_genes) %>%
  transmute(
    Gene = geneU,
    Gene_norm = geneU,
    variant_id = paste0(chrom, ":", pos, ":", ref, ":", alt),
    consequence = "truncating",  # the v2.1.1 bgz consequence already filtered
    pos_source = "VEP_canonical_protein_start",
    raw_pos = as.character(aa_pos),
    aa_pos = as.numeric(aa_pos),
    prot_len_aa = as.numeric(prot_len_aa),
    rel_pos = as.numeric(rel_pos),
    group = ifelse(geneU %in% disc35, "Discordant", "Background"),
    pos_bin = bin_relpos(rel_pos)
  )
message("[relpos] Final trunc_pdn rows: ", nrow(trunc_pdn))
trunc_pdn$group <- as.character(trunc_pdn$group)
message("[relpos] group counts: ", paste(names(table(as.character(trunc_pdn$group))), as.integer(table(as.character(trunc_pdn$group))), sep="=", collapse=" | "))

# write + cache
readr::write_csv(trunc_pdn, OUT_TRUNC_PDN_VARIANTS)
message("Wrote: ", OUT_TRUNC_PDN_VARIANTS)

summ <- trunc_pdn %>%
  group_by(group) %>%
  summarise(
    n_variants_total = n(),
    n_with_relpos = sum(!is.na(rel_pos)),
    relpos_median = median(rel_pos, na.rm = TRUE),
    relpos_IQR = IQR(rel_pos, na.rm = TRUE),
    early = sum(pos_bin == "early", na.rm = TRUE),
    middle = sum(pos_bin == "middle", na.rm = TRUE),
    late = sum(pos_bin == "late", na.rm = TRUE),
    prop_late = late / pmax(1, early + middle + late),
    .groups = "drop"
  )
readr::write_csv(summ, OUT_TRUNC_PDN_SUMMARY)
message("Wrote: ", OUT_TRUNC_PDN_SUMMARY)

saveRDS(trunc_pdn, TRUNC_PDN_RDS)
message("[CACHE] Saved: ", TRUNC_PDN_RDS)

message("[relpos] Done writing positional CSV + summary. Stopping here (no downstream code).")
invisible(NULL)
invisible(NULL)
