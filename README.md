Constraint paper pipeline

Reproducible analysis pipeline for the manuscript
“Constraint does not equal penetrance: population loss-of-function variation in ultra-constrained human genes.”

This repository contains the full R pipeline used to generate all main figures and tables for the manuscript.

TL;DR — run everything
git clone git@github.com:lecb/constraint_paper.git
cd constraint_paper
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages("renv"); renv::restore()'
bash scripts/bootstrap_inputs.sh
Rscript scripts/run_all.R


Final outputs appear in:

final/figures
final/tables

Contents

Quick start

Required inputs

Running the pipeline

Outputs

Reproducibility

Repository layout

Troubleshooting

Quick start
1) Clone the repository
git clone git@github.com:lecb/constraint_paper.git
cd constraint_paper

2) Restore the R environment (renv)

This project uses renv to lock package versions.

Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages("renv"); renv::restore()'


This installs the exact R package versions used in the paper.

Required inputs

Large data files are not included in the repository due to size/licensing.

The pipeline expects the following files in the repo root:

gnomad.v4.1.constraint_metrics.tsv
gnomad_variants_all.csv
trunc_position_domain_nmd_variants.csv


To help users, a bootstrap script is provided:

bash scripts/bootstrap_inputs.sh


You may instead symlink your own local copies.

Running the pipeline

Run the full cached pipeline:

Rscript scripts/run_all.R


Force recompute all steps:

FORCE_RECOMPUTE=TRUE Rscript scripts/run_all.R

Outputs

At the end of the run, all figures and tables are collected into:

final/
├── figures/
└── tables/


Intermediate outputs may also appear in:

figures/
tables/
manuscript_figures/
cache/
logs/

Reproducibility

This project is designed for full computational reproducibility.

R package versions are pinned in renv.lock

The pipeline is deterministic given fixed inputs and seeds

All analyses can be regenerated with a single command

Repository layout
scripts/
  run_all.R                     Main pipeline entrypoint
  core/                         Paper-critical scripts
  optional/                     Additional analyses (not run by default)
  _legacy/                      Archived historical scripts

renv.lock                       Pinned R package environment
renv/                           renv infrastructure
final/                          Collected outputs (after pipeline run)

Troubleshooting
“Package was built under R version …”

These warnings are normal and safe to ignore.

Missing input files

Confirm required files exist (see Required inputs section).

Pipeline finishes but outputs are missing

Check the final collection script:

scripts/core/99_collect_outputs.sh

Citation

If you use this code, please cite the associated manuscript.
Full citation details will be added upon publication.