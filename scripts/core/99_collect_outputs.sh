#!/usr/bin/env bash
set -euo pipefail

ROOT="${PIPELINE_ROOT:-$(pwd)}"
cd "$ROOT"

OUT_ROOT="final"
FIG_OUT="${OUT_ROOT}/figures"
TAB_OUT="${OUT_ROOT}/tables"

mkdir -p "$FIG_OUT" "$TAB_OUT"

SEARCH_DIRS=(
  "figures"
  "tables"
  "manuscript_figures"
  "gnomad_lof_discordance_out"
  "positional_tests_out"
  "results"
  "output"
)

copy_with_path() {
  local src="$1"
  local dest_dir="$2"
  local rel="${src#./}"
  local flat="${rel//\//__}"
  local dest="${dest_dir}/${flat}"

  if [[ -e "$dest" ]]; then
    local base="${dest%.*}"
    local ext="${dest##*.}"
    local i=2
    while [[ -e "${base}__${i}.${ext}" ]]; do i=$((i+1)); done
    dest="${base}__${i}.${ext}"
  fi

  cp -p "$src" "$dest"
}

for d in "${SEARCH_DIRS[@]}"; do
  [[ -d "$d" ]] || continue
  while IFS= read -r -d '' f; do
    copy_with_path "$f" "$FIG_OUT"
  done < <(find "$d" -type f \( -iname "*.png" -o -iname "*.pdf" -o -iname "*.svg" \) -print0)
done

for d in "${SEARCH_DIRS[@]}"; do
  [[ -d "$d" ]] || continue
  while IFS= read -r -d '' f; do
    copy_with_path "$f" "$TAB_OUT"
  done < <(find "$d" -type f \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" -o -iname "*.xlsx" -o -iname "*.docx" \) -print0)
done

echo "[collect] final outputs:"
echo "  ${FIG_OUT}  ($(find "$FIG_OUT" -type f | wc -l | tr -d " ") files)"
echo "  ${TAB_OUT}  ($(find "$TAB_OUT" -type f | wc -l | tr -d " ") files)"
