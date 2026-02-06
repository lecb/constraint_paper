#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ----------------------------
# INPUTS
# ----------------------------
EXPORT_DIR <- "discordant_genes"       # folder containing 35 LoF-only gene CSVs (gnomAD browser exports)
DISCORDANT_RDS <- "discordant_genes.rds"
OUTDIR <- "tables"

# IMPORTANT: “rare homozygous” analysis should NOT count common truncating artefacts
# (e.g. WNK1 / PPFIA1 common frameshifts in alt/low-confidence contexts).
# Use AF filter from browser export.
RARE_AF_CUTOFF <- 1e-3   # 0.001; set to 1e-4 if you want ultra-rare

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 1) Load discordant gene list (35)
# ----------------------------
disc <- readRDS(DISCORDANT_RDS)
stopifnot("Gene" %in% names(disc))
disc35 <- unique(as.character(disc$Gene))
disc35_u <- toupper(trimws(disc35))

message("[INFO] discordant genes in RDS = ", length(disc35))

# ----------------------------
# 2) Load all LoF-only exports
# ----------------------------
files <- list.files(EXPORT_DIR, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
message("[INFO] CSV files found in ", EXPORT_DIR, " = ", length(files))

read_gene_export <- function(path) {
  x <- suppressMessages(read_csv(path, show_col_types = FALSE))
  
  # Identify gene column (varies depending on browser export)
  gene_col <- intersect(names(x), c("Gene", "gene", "gene_symbol", "Symbol", "symbol"))
  if (length(gene_col) >= 1) {
    x <- x %>% mutate(Gene = as.character(.data[[gene_col[[1]]]]))
  } else {
    # fallback: infer from filename
    x <- x %>% mutate(Gene = tools::file_path_sans_ext(basename(path)))
  }
  
  x %>% mutate(.src_file = basename(path))
}

raw <- map_dfr(files, read_gene_export) %>%
  mutate(Gene_u = toupper(trimws(Gene)))

# Keep only the discordant 35 genes (in case extra files in folder)
raw <- raw %>% filter(Gene_u %in% disc35_u)

message("[INFO] rows loaded (discordant-only) = ", nrow(raw))
message("[INFO] genes present in exports       = ", n_distinct(raw$Gene_u))

# Warn if file set doesn’t match the discordant set
present_genes <- sort(unique(raw$Gene_u))
missing_genes <- setdiff(sort(unique(disc35_u)), present_genes)
extra_genes   <- setdiff(present_genes, sort(unique(disc35_u)))
if (length(missing_genes) > 0) warning("[WARN] Missing exports for: ", paste(missing_genes, collapse = ", "))
if (length(extra_genes) > 0) warning("[WARN] Extra exports present (not in RDS): ", paste(extra_genes, collapse = ", "))

# ----------------------------
# 3) Robust column name mapping
# ----------------------------
pick_col <- function(df, candidates) {
  hit <- intersect(names(df), candidates)
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

colmap <- list(
  variant_id   = c("gnomAD ID", "Variant ID", "variant_id"),
  vep          = c("VEP Annotation", "VEP consequence", "Consequence", "consequence"),
  tx_conseq    = c("Transcript Consequence", "Transcript consequence", "transcript_consequence"),
  flags        = c("Flags", "flags"),
  AC           = c("Allele Count", "AC", "ac"),
  AF           = c("Allele Frequency", "AF", "af"),
  NHOM         = c("Homozygote Count", "Homozygotes", "nhomalt", "NHOMALT"),
  NHEMI        = c("Hemizygote Count", "Hemizygotes", "hemi_count", "NHEMI")
)

vcol    <- pick_col(raw, colmap$variant_id)
vepcol  <- pick_col(raw, colmap$vep)
tccol   <- pick_col(raw, colmap$tx_conseq)  # optional
fcol    <- pick_col(raw, colmap$flags)
accol   <- pick_col(raw, colmap$AC)
afcol   <- pick_col(raw, colmap$AF)
hcol    <- pick_col(raw, colmap$NHOM)       # optional (usually present)
hemicol <- pick_col(raw, colmap$NHEMI)      # optional (usually present)

# Require these minimum columns
need_present <- c(!is.null(vcol), !is.null(vepcol), !is.null(fcol), !is.null(accol), !is.null(afcol))
if (!all(need_present)) {
  stop(
    "Missing one or more required columns in exports.\n",
    "Found: variant_id=", vcol %||% "NONE",
    " | VEP Annotation=", vepcol %||% "NONE",
    " | flags=", fcol %||% "NONE",
    " | AC=", accol %||% "NONE",
    " | AF=", afcol %||% "NONE",
    "\nTip: run names(read_csv('discordant_genes/<onefile>.csv')) to inspect headers."
  )
}

message("[INFO] Using columns: ",
        "variant_id=", vcol, " | vep=", vepcol, " | flags=", fcol,
        " | AC=", accol, " | AF=", afcol,
        if (!is.null(hcol)) paste0(" | NHOM=", hcol) else "",
        if (!is.null(hemicol)) paste0(" | NHEMI=", hemicol) else "")

# ----------------------------
# 4) Build unified working table and de-duplicate transcript rows -> variant-level
#    IMPORTANT: apply rarity filter for homozygote claims
# ----------------------------
df_var <- raw %>%
  transmute(
    Gene   = Gene,
    Gene_u = Gene_u,
    variant_id = as.character(.data[[vcol]]),
    
    # use VEP annotation as the primary consequence source (most reliable)
    vep = as.character(.data[[vepcol]] %||% ""),
    
    # optional: include transcript consequence for extra context (not used for trunc call)
    tx_consequence = if (!is.null(tccol)) as.character(.data[[tccol]] %||% "") else "",
    
    flags = as.character(.data[[fcol]] %||% ""),
    AC    = suppressWarnings(as.numeric(.data[[accol]] %||% 0)),
    AF    = suppressWarnings(as.numeric(.data[[afcol]] %||% NA_real_)),
    NHOM  = if (!is.null(hcol)) suppressWarnings(as.numeric(.data[[hcol]] %||% 0)) else NA_real_,
    NHEMI = if (!is.null(hemicol)) suppressWarnings(as.numeric(.data[[hemicol]] %||% 0)) else NA_real_
  ) %>%
  mutate(
    vep_l = tolower(vep %||% ""),
    is_trunc = str_detect(vep_l, "stop_gained|frameshift_variant"),
    
    flags_l = tolower(flags %||% ""),
    in_bad  = str_detect(flags_l, "segdup") | str_detect(flags_l, "\\blcr\\b"),
    
    # rarity filter for homozygote/hemizygote “rare genotype” claims
    is_rare = is.na(AF) | AF <= RARE_AF_CUTOFF
  ) %>%
  # Apply rarity filter globally for this table (keeps your downstream summaries consistent
  # with the “rare homozygous truncating genotypes” statement)
  filter(is_rare) %>%
  group_by(Gene_u, Gene, variant_id) %>%
  summarise(
    # take max across transcript rows (prevents double-counting)
    AC    = max(AC, na.rm = TRUE),
    AF    = suppressWarnings(min(AF, na.rm = TRUE)),
    NHOM  = if (all(is.na(NHOM))) NA_real_ else max(NHOM, na.rm = TRUE),
    NHEMI = if (all(is.na(NHEMI))) NA_real_ else max(NHEMI, na.rm = TRUE),
    
    # “bad region” if any row for the variant has flags
    in_bad   = any(in_bad, na.rm = TRUE),
    
    # trunc if any row for the variant is stop/frameshift by VEP annotation
    is_trunc = any(is_trunc, na.rm = TRUE),
    
    .groups = "drop"
  )

stopifnot(all(c("Gene","variant_id","AC","AF","NHOM","NHEMI","in_bad","is_trunc") %in% names(df_var)))
message("[INFO] df_var rows = ", nrow(df_var),
        " | genes = ", n_distinct(df_var$Gene_u),
        " | rare AF cutoff = ", RARE_AF_CUTOFF)

# ----------------------------
# 5) Summarise per gene (pLoF and truncating), flag-aware
# ----------------------------
summarise_set <- function(df, which = c("plof","trunc")) {
  which <- match.arg(which)
  
  # “plof” set: LoF-only export (already LoF-filtered upstream)
  df0 <- if (which == "plof") df else df %>% filter(is_trunc)
  
  pre <- df0 %>%
    group_by(Gene) %>%
    summarise(
      n_rows_type = n(),
      n_rows_removed_by_flags = sum(in_bad, na.rm = TRUE),
      .groups = "drop"
    )
  
  post <- df0 %>%
    filter(!in_bad) %>%
    group_by(Gene) %>%
    summarise(
      n_variants_flagaware = n_distinct(variant_id),
      AC_sum_flagaware = sum(AC, na.rm = TRUE),
      hom_sum_flagaware = if (all(is.na(NHOM))) NA_real_ else sum(NHOM, na.rm = TRUE),
      hemi_sum_flagaware = if (all(is.na(NHEMI))) NA_real_ else sum(NHEMI, na.rm = TRUE),
      classes_flagaware = if (which == "plof") "LoF-only export" else "stop/frameshift",
      .groups = "drop"
    )
  
  full_join(pre, post, by = "Gene") %>%
    mutate(
      n_rows_type = coalesce(n_rows_type, 0),
      n_rows_removed_by_flags = coalesce(n_rows_removed_by_flags, 0),
      n_variants_flagaware = coalesce(n_variants_flagaware, 0),
      AC_sum_flagaware = coalesce(AC_sum_flagaware, 0),
      classes_flagaware = coalesce(classes_flagaware, "")
    )
}

plof_sum35  <- summarise_set(df_var, "plof")
trunc_sum35 <- summarise_set(df_var, "trunc")

# Write the summaries
write_csv(plof_sum35,  file.path(OUTDIR, "flagaware_plof_summary_discordant35.csv"),  na = "")
write_csv(trunc_sum35, file.path(OUTDIR, "flagaware_trunc_summary_discordant35.csv"), na = "")

# ----------------------------
# 6) Build combined supplementary table (same structure as your top15 table)
# ----------------------------
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

write_csv(supp35, file.path(OUTDIR, "SuppTable_flagaware_variant_filtering_discordant35.csv"), na = "")

message("[DONE] wrote:")
message(" - ", file.path(OUTDIR, "flagaware_plof_summary_discordant35.csv"))
message(" - ", file.path(OUTDIR, "flagaware_trunc_summary_discordant35.csv"))
message(" - ", file.path(OUTDIR, "SuppTable_flagaware_variant_filtering_discordant35.csv"))
message("[DONE] n genes in supp35 = ", nrow(supp35))

# ----------------------------
# 7) Print homozygote hits for Results sentence (rare + flag-aware)
# ----------------------------
hom_hits <- supp35 %>%
  select(Gene, trunc_hom_sum_flagaware) %>%
  filter(!is.na(trunc_hom_sum_flagaware), trunc_hom_sum_flagaware > 0) %>%
  arrange(desc(trunc_hom_sum_flagaware))

message("\n[HOMOZYGOTE HITS] truncating (stop/frameshift), rare (AF <= ",
        RARE_AF_CUTOFF, "), and flag-aware (!segdup/!LCR):")
print(hom_hits, n = 50)


#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(purrr)
})

EXPORT_DIR <- "discordant_genes"
DISCORDANT_RDS <- "discordant_genes.rds"
OUTDIR <- "tables"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Load discordant genes (35)
# ----------------------------
disc <- readRDS(DISCORDANT_RDS)
disc35 <- unique(as.character(disc$Gene))

# ----------------------------
# Read all exports
# ----------------------------
files <- list.files(EXPORT_DIR, pattern="\\.csv$", full.names=TRUE)

raw <- map_dfr(files, ~{
  x <- read_csv(.x, show_col_types = FALSE)
  gene <- tools::file_path_sans_ext(basename(.x))
  x$Gene <- gene
  x
})

raw <- raw %>% filter(Gene %in% disc35)

# ----------------------------
# Identify ancestry AC columns
# ----------------------------
ac_cols <- grep("^Allele Count ", names(raw), value = TRUE)

stopifnot(length(ac_cols) > 0)

# ----------------------------
# Restrict to truncating variants only
# ----------------------------
df <- raw %>%
  mutate(
    vep_l = tolower(`VEP Annotation`),
    is_trunc = str_detect(vep_l, "stop_gained|frameshift_variant")
  ) %>%
  filter(is_trunc)

# ----------------------------
# Count ancestries per variant
# ----------------------------
ancestry_counts <- df %>%
  rowwise() %>%
  mutate(
    n_ancestries_with_ac = sum(c_across(all_of(ac_cols)) > 0, na.rm = TRUE)
  ) %>%
  ungroup()

# Variant-level summary table
supp <- ancestry_counts %>%
  transmute(
    Gene,
    variant_id = `gnomAD ID`,
    n_ancestries_with_ac
  )

# Write full table
write_csv(
  supp,
  file.path(OUTDIR, "SuppTable_ancestry_distribution_discordant35.csv"),
  na = ""
)

# ----------------------------
# Print key manuscript statistic
# ----------------------------
median_anc <- median(supp$n_ancestries_with_ac, na.rm = TRUE)
mean_anc   <- mean(supp$n_ancestries_with_ac, na.rm = TRUE)

message("\n[DONE]")
message("Variants analysed = ", nrow(supp))
message("Median ancestries per trunc variant = ", round(median_anc, 2))
message("Mean ancestries per trunc variant   = ", round(mean_anc, 2))
