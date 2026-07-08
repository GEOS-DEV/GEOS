#!/usr/bin/env bash
set -euo pipefail

DOC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${DOC_DIR}"

REFRESH_DATA=0
if [[ "${1:-}" == "--refresh-data" ]]; then
  REFRESH_DATA=1
  shift
fi
if [[ "$#" -ne 0 ]]; then
  echo "Usage: $0 [--refresh-data]" >&2
  exit 2
fi

determine_repo_root() {
  if [[ -n "${GEOS_REPO_ROOT:-}" ]]; then
    printf '%s\n' "${GEOS_REPO_ROOT}"
  elif git -C "${DOC_DIR}" rev-parse --show-toplevel >/dev/null 2>&1; then
    git -C "${DOC_DIR}" rev-parse --show-toplevel
  elif [[ -d "${DOC_DIR}/../../../../src" && -d "${DOC_DIR}/../../../../scripts" ]]; then
    cd "${DOC_DIR}/../../../.." && pwd
  else
    echo "Could not determine GEOS repository root." >&2
    echo "Set GEOS_REPO_ROOT=/path/to/geos or rerun without --refresh-data if generated/data is already present." >&2
    exit 2
  fi
}

mkdir -p generated/data generated/tex
if [[ "${REFRESH_DATA}" -eq 1 || ! -f generated/data/geos_mpm_extracted.json ]]; then
  REPO_ROOT="$(determine_repo_root)"
  python3 tools/extract_mpm_metadata.py \
    --repo "${REPO_ROOT}" \
    --out generated/data/geos_mpm_extracted.json
fi
python3 tools/render_generated_tex.py \
  --in generated/data/geos_mpm_extracted.json \
  --out generated/tex

rm -f geos_mpm_manual.aux geos_mpm_manual.toc geos_mpm_manual.out \
      geos_mpm_manual.idx geos_mpm_manual.ind geos_mpm_manual.ilg \
      geos_mpm_manual.log geos_mpm_manual.pdf geos_mpm_manual_normalized.pdf

pdflatex -interaction=nonstopmode -halt-on-error geos_mpm_manual.tex
makeindex geos_mpm_manual.idx
pdflatex -interaction=nonstopmode -halt-on-error geos_mpm_manual.tex
pdflatex -interaction=nonstopmode -halt-on-error geos_mpm_manual.tex

# The direct pdfTeX output is the deliverable. If a downstream workflow requires
# additional PDF normalization, run it as a separate post-processing step and verify
# the normalized PDF with pdfinfo/rendering before replacing this file.
