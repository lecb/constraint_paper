#!/usr/bin/env Rscript
# make_supp_sheet.R
#
# Build an Excel supplementary workbook from specific CSV/TSV files found anywhere
# under a root directory, writing each file into its own worksheet named
# "Supplementary Table S1" ... "Supplementary Table S11".
#
# Usage:
#   Rscript make_supp_sheet.R \
#     --root "/path/to/project" \
#     --out Supplementary_Tables.xlsx \
#     --prefer newest
#
# Notes:
# - Searches subdirectories (recursive=TRUE)
# - If duplicate filenames are found, chooses one based on --prefer (newest|largest|first)
# - Fixes duplicate/blank column names so openxlsx doesn't crash.

suppressPackageStartupMessages({
  library(readr)
  library(openxlsx)
})

# ----------------------------
# Config: filename -> sheet name
# ----------------------------
file_to_sheet <- c(
  "Supplementary_Table1_discordant_genes.csv" = "Supplementary Table S1",
  "sensitivity_counts.csv" = "Supplementary Table S2",
  "threshold_sensitivity_summary.csv" = "Supplementary Table S3",
  "SuppTable_splice_free_stop_frameshift_counts.csv" = "Supplementary Table S4",
  "SuppTable_flagaware_variant_filtering_discordant35.csv" = "Supplementary Table S5",
  "SuppTable_ancestry_distribution_discordant35.csv" = "Supplementary Table S6",
  "Supp_combined_universe_flags_disease.csv" = "Supplementary Table S7",
  "Supplementary_Table_ClinVar_PLP_truncating_discordant_genes_GE2stars.tsv" = "Supplementary Table S8",
  "SuppTable_domain_structured_truncation_features.csv" = "Supplementary Table S9",
  "SuppTable_NMD_escape_summary.csv" = "Supplementary Table S10",
  "SuppTable_truncation_position_tests.tsv" = "Supplementary Table S11"
)

# ----------------------------
# Arg parsing (no extra deps)
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop("Missing value after ", flag)
  args[i + 1]
}

root     <- get_arg("--root", default = NULL)
out_xlsx <- get_arg("--out", default = "Supplementary_Tables.xlsx")
prefer   <- get_arg("--prefer", default = "newest") # newest | largest | first

if (is.null(root)) {
  stop("Usage: Rscript make_supp_sheet.R --root /path --out file.xlsx [--prefer newest|largest|first]")
}

if (!prefer %in% c("newest", "largest", "first")) {
  stop("--prefer must be one of: newest, largest, first")
}

root <- normalizePath(root, mustWork = TRUE)

# ----------------------------
# Helpers
# ----------------------------
is_csv <- function(p) grepl("\\.csv(\\.gz)?$", tolower(p))
is_tsv <- function(p) grepl("\\.tsv(\\.gz)?$", tolower(p))

read_table_auto <- function(path) {
  # Keep everything as character so IDs don’t become scientific notation etc.
  if (is_tsv(path)) {
    readr::read_tsv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE,
      na = character()
    )
  } else if (is_csv(path)) {
    readr::read_csv(
      path,
      col_types = readr::cols(.default = readr::col_character()),
      progress = FALSE,
      na = character()
    )
  } else {
    stop("Unsupported file type: ", path)
  }
}

choose_one <- function(paths, prefer = "newest") {
  if (length(paths) == 1) return(paths[[1]])
  if (prefer == "first") return(sort(paths)[1])
  
  info <- file.info(paths)
  if (prefer == "largest") return(paths[[which.max(info$size)]])
  # newest
  paths[[which.max(info$mtime)]]
}

make_unique_colnames <- function(nms) {
  # Trim whitespace and replace blanks
  nms <- trimws(nms)
  nms[nms == ""] <- "X"
  
  # Make exact duplicates unique
  nms <- make.unique(nms, sep = "_dup")
  
  # Protect against case-insensitive collisions (e.g., "Gene" vs "gene")
  lower <- tolower(nms)
  if (any(duplicated(lower))) {
    idx <- ave(seq_along(nms), lower, FUN = seq_along)
    nms <- ifelse(idx > 1, paste0(nms, "_case", idx), nms)
  }
  nms
}

# ----------------------------
# Find targets under root
# ----------------------------
targets <- names(file_to_sheet)

all_files <- list.files(
  root,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE,
  include.dirs = FALSE
)
all_files <- all_files[is_csv(all_files) | is_tsv(all_files)]

basenames <- basename(all_files)

found <- setNames(vector("list", length(targets)), targets)
for (fn in targets) {
  found[[fn]] <- all_files[basenames == fn]
}

missing <- targets[lengths(found) == 0]
dupes   <- targets[lengths(found) > 1]

if (length(missing) > 0) {
  message("WARNING: missing files (not found under root):")
  for (fn in missing) message("  - ", fn)
}

if (length(dupes) > 0) {
  message("\nWARNING: multiple matches found for these filenames:")
  for (fn in dupes) {
    message("  - ", fn)
    for (p in found[[fn]]) message("      ", p)
  }
  message("\n(Using --prefer ", prefer, " to choose a single match per filename.)")
}

# ----------------------------
# Write workbook
# ----------------------------
wb <- openxlsx::createWorkbook()

# Order sheets S1..S11 by numeric suffix
sheet_nums <- as.integer(sub(".*S", "", unname(file_to_sheet)))
ordered_fns <- names(file_to_sheet)[order(sheet_nums)]

for (fn in ordered_fns) {
  paths <- found[[fn]]
  if (length(paths) == 0) next
  
  chosen <- choose_one(paths, prefer = prefer)
  sheet  <- unname(file_to_sheet[[fn]])
  
  message("Adding ", chosen, "  -->  '", sheet, "'")
  
  # Wrap each sheet in a tryCatch so errors tell you exactly which file failed
  tryCatch({
    df <- read_table_auto(chosen)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    names(df) <- make_unique_colnames(names(df))
    
    openxlsx::addWorksheet(wb, sheetName = sheet)
    
    # If there are zero columns (rare edge case), write an empty placeholder
    if (ncol(df) == 0) {
      openxlsx::writeData(wb, sheet = sheet, x = data.frame(Note = "No columns detected"), withFilter = FALSE)
    } else {
      openxlsx::writeDataTable(
        wb,
        sheet = sheet,
        x = df,
        withFilter = TRUE,
        tableStyle = "TableStyleMedium2"
      )
      # Auto column widths (can be slow for very wide tables)
      openxlsx::setColWidths(wb, sheet = sheet, cols = 1:ncol(df), widths = "auto")
    }
  }, error = function(e) {
    stop("Failed on sheet: ", sheet, " | file: ", chosen, "\nReason: ", conditionMessage(e))
  })
}

openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("\nDone. Wrote: ", normalizePath(out_xlsx, mustWork = FALSE))
