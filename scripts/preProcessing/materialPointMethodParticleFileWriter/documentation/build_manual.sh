#!/usr/bin/env bash
set -euo pipefail

python3 make_manual.py

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
