# Constraint paper pipeline
 
**Constraint does not equal penetrance: population loss-of-function variation in ultra-constrained human genes.**

This repository contains the full R pipeline used to generate all main figures and tables for the manuscript.

---

## Quick start

```bash
git clone git@github.com:lecb/constraint_paper.git
cd constraint_paper
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages("renv"); renv::restore()'
bash scripts/bootstrap_inputs.sh
Rscript scripts/run_all.R
```

Final outputs appear in:

```
final/figures
final/tables
```

---

## Contents
- Required inputs
- Running the pipeline
- Outputs
- Reproducibility
- Repository layout
- Troubleshooting
- Citation

---

## Required inputs

Large datasets are not committed to the repository (size/licensing).  
The pipeline expects the following files in the repo root:

- gnomad.v4.1.constraint_metrics.tsv
- gnomad_variants_all.csv
- trunc_position_domain_nmd_variants.csv

Helper script:

```bash
bash scripts/bootstrap_inputs.sh
```

You may also symlink local copies.

---

## Running the pipeline

Run the cached pipeline:

```bash
Rscript scripts/run_all.R
```

Force recomputation:

```bash
FORCE_RECOMPUTE=TRUE Rscript scripts/run_all.R
```

---

## Outputs

All figures and tables are collected into:

```text
final/
├── figures/
└── tables/
```

Intermediate artefacts may also appear in:

- cache/
- logs/
- figures/
- tables/
- manuscript_figures/

---

## Reproducibility

This project uses renv to pin package versions via renv.lock.

Restore the environment:

```r
renv::restore()
```

---

## Repository layout

```text
scripts/
  run_all.R                      Main pipeline entrypoint
  core/                          Paper-critical analyses
  optional/                      Additional analyses (not run by default)
  _legacy/                       Archived development scripts

renv.lock                         Pinned R package environment
renv/                             renv infrastructure
```

---

## Troubleshooting

"package … was built under R version …" warnings are expected and usually safe to ignore.

If GitHub says "fetch first", run:

```bash
git pull --rebase origin main
git push
```

---

## Citation

If you use this code, please cite the associated manuscript.
