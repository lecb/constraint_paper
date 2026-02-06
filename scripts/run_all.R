## ---- portable root bootstrap ----
args <- commandArgs(trailingOnly = TRUE)
if ("--root" %in% args) {
  ROOT <- args[match("--root", args) + 1]
  Sys.setenv(PIPELINE_ROOT = ROOT)
  setwd(ROOT)
} else {
  ROOT <- normalizePath(getwd(), mustWork = TRUE)
  Sys.setenv(PIPELINE_ROOT = ROOT)
  setwd(ROOT)
}
## -----------------------------------

#!/usr/bin/env Rscript

# Cached, reproducible pipeline runner for constraint paper.
# Usage:
#   Rscript repo/scripts/run_all.R
#   FORCE_RECOMPUTE=TRUE Rscript repo/scripts/run_all.R   # rebuild all steps

suppressPackageStartupMessages({
  library(tools)
})

PIPE <- new.env(parent = globalenv())

force_recompute <- identical(Sys.getenv("FORCE_RECOMPUTE"), "TRUE")

message("[run_all] Starting pipeline...")
message("[run_all] FORCE_RECOMPUTE=", force_recompute)

# Ensure expected output dirs exist (many scripts assume these)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("cache", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("gnomad_lof_discordance_out", showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("gnomad_lof_discordance_out", "cache"), showWarnings = FALSE, recursive = TRUE)
dir.create("positional_tests_out", showWarnings = FALSE, recursive = TRUE)
dir.create("pext", showWarnings = FALSE, recursive = TRUE)

should_run <- function(outputs) {
  # If no declared outputs, always run (but try to declare outputs for caching!)
  if (length(outputs) == 0) return(TRUE)
  if (force_recompute) return(TRUE)
  !all(file.exists(outputs))
}

run_step <- function(name, rel_script_path, outputs = character(), optional = FALSE) {
  script_path <- file.path(ROOT, "scripts", rel_script_path)
  if (!file.exists(script_path)) {
    msg <- paste0("[run_all] Missing script: ", rel_script_path)
    if (optional) { message(msg, " (skip)"); return(invisible(FALSE)) }
    stop(msg)
  }

  if (!should_run(outputs)) {
    message("\n[run_all] Skip: ", name, " (outputs exist)")
    return(invisible(TRUE))
  }

  message("\n[run_all] Running: ", name, " -> ", rel_script_path)
  ok <- TRUE
  tryCatch(
    source(script_path, local = PIPE),
    error = function(e) {
      ok <<- FALSE
      msg <- paste0("[run_all] ERROR in ", name, ": ", conditionMessage(e))
      if (optional) message(msg) else stop(msg)
    }
  )

  if (ok && length(outputs)) {
    missing <- outputs[!file.exists(outputs)]
    if (length(missing)) {
      msg <- paste0("[run_all] WARNING: step ran but some outputs missing for ", name,
                    ":\n  - ", paste(missing, collapse = "\n  - "))
      if (optional) message(msg) else stop(msg)
    }
  }
  invisible(ok)
}

# ----------------------------
# Core steps (paper-critical)
# ----------------------------
run_step(
  name = "constraint_universe (S1 + QC + Table1 top15)",
  rel_script_path = "core/01_constraint_universe.R",
  outputs = c(
    "tables/Supplementary_Table1_discordant_genes.csv",
    "tables/loeuf_lt_0.2_precision_filtered.csv"
  )
)

run_step(
  name = "relpos / trunc_position caches",
  rel_script_path = "core/02_generate_relpos.R",
  outputs = c(
    "cache/trunc_pdn_cache.rds",
    "trunc_position_domain_nmd_variants.csv",
    "trunc_position_domain_nmd_summary_by_group.csv"
  )
)

run_step(
  name = "pext analysis (Fig3 + pext tables)",
  rel_script_path = "core/05_pext_analysis.R",
  outputs = c(
    "figures/SuppFig_pext_expression_attenuation_discordant35_2panel.png",
    "tables/SuppTable_pext_variant_level_discordant_vs_bg.tsv.gz",
    "tables/SuppTable_pext_gene_level_discordant_vs_bg.tsv.gz"
  )
)

run_step(
  name = "threshold + positional sensitivity (S2/S3 + S11)",
  rel_script_path = "core/06_threshold_and_position_sensitivity.R",
  outputs = c(
    "tables/threshold_sensitivity_summary.csv",
    "tables/sensitivity_counts.csv",
    "tables/SuppFig_trunc_position_ECDF_and_DeltaCDF_NATIVE.png",
    "tables/SuppTable_truncation_position_tests.tsv"
  )
)

run_step(
  name = "SuppFig S1 robustness superfigure",
  rel_script_path = "core/07_make_SuppFig_S1_robustness.R",
  outputs = c(
    "gnomad_lof_discordance_out/SuppFig_super_discordance_robustness.png"
  )
)

run_step(
  name = "Main figs + key supp superfigs (Fig1/Fig2/Supp S2)",
  rel_script_path = "core/09_make_main_and_key_supp_figs.R",
  outputs = c(
    "manuscript_figures/main/Figure2_discordant_vs_rest.png",
    "manuscript_figures/main/SuperFig_discordant_summary.png"
  )
)

message("\n[run_all] Core pipeline complete.")


message("\n[run_all] Optional pipeline removed.")

message("\n[run_all] Done.")

# ----------------------------
# Final: collect outputs into single dirs
# ----------------------------
message("\n[run_all] Collecting outputs into final/figures and final/tables ...")
system("bash scripts/core/99_collect_outputs.sh", intern = FALSE)
