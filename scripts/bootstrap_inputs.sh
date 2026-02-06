#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-.}"
cd "$ROOT"

mkdir -p data/{gnomad,clinvar,pext,domains}

echo "[bootstrap] Place or download required inputs into repo/data/ ..."
echo "TODO: add curl/wget commands + checksums for:"
echo "  - gnomAD v4.1 constraint table -> data/gnomad/..."
echo "  - pext GRCh38 -> data/pext/..."
echo "  - ClinVar VCF+TBI -> data/clinvar/clinvar.vcf.gz(.tbi)"
echo "  - InterPro/BioMart domains -> data/domains/domains_interpro_biomart.csv"
echo ""
echo "For now, if you already have them, copy them into these locations."

