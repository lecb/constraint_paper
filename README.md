Constraint does not equal penetrance

Reproducible analysis pipeline

This repository contains the full reproducible analysis pipeline used for the manuscript:

“Constraint does not equal penetrance: population loss-of-function variation in ultra-constrained human genes.”

The pipeline processes gnomAD constraint and variant data, performs downstream analyses, generates all manuscript figures and tables, and collects them into a single output directory.

Overview

Running the pipeline will automatically:

Build the ultra-constrained gene universe

Generate truncation position and NMD analyses

Perform pext expression analyses

Run robustness and sensitivity analyses

Generate all main and supplementary figures

Collect all outputs into a single directory

Final outputs appear in:

final/figures/
final/tables/

Requirements

Tested on:

macOS and Linux

R ≥ 4.3 (recommended 4.4+)

All R packages are installed automatically via renv.

Quick start (from a clean machine)
1) Clone the repository
git clone https://github.com/YOUR_USERNAME/constraint-paper.git
cd constraint-paper

2) Install R dependencies
Rscript -e 'install.packages("renv"); renv::restore()'


This installs the exact package versions used to generate the manuscript.

3) Prepare required input data

Due to size and licensing, large gnomAD-derived inputs are not included in the repository.

Create the input directory:

mkdir -p data/gnomad


Place the following files in the repository root (or symlink them):

File	Description
gnomad.v4.1.constraint_metrics.tsv	gnomAD v4.1 constraint metrics
gnomad_variants_all.csv	Pre-extracted gnomAD pLoF variants
trunc_position_domain_nmd_variants.csv	Truncating variant positional/NMD annotations

These are the only required external inputs to run the pipeline.

4) Run the full pipeline
Rscript scripts/run_all.R


The pipeline uses caching, so reruns are fast.

Output

After completion, all outputs are collected into:

final/
├── figures/
└── tables/


This directory contains all figures and tables required for the manuscript.

Re-running the pipeline

Force a complete rebuild:

FORCE_RECOMPUTE=TRUE Rscript scripts/run_all.R

Repository structure
scripts/
  core/        → core analysis scripts
  run_all.R    → main pipeline entrypoint

figures/       → intermediate figures
tables/        → intermediate tables
cache/         → cached objects

final/         → collected outputs (figures + tables)


The pipeline automatically:

Runs all core analyses

Generates manuscript figures

Collects outputs into a single directory

Reproducibility

This repository ensures reproducibility via:

deterministic caching

pinned R package versions (renv)

scripted figure generation

single-command execution

Running the pipeline from scratch reproduces all outputs.

Running the pipeline in practice

Typical workflow:

git pull
Rscript scripts/run_all.R

Contact

For questions about the pipeline or manuscript, please contact the corresponding author.