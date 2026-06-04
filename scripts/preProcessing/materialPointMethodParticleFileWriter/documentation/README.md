
# GEOS-MPM Starter Manual source bundle

This bundle contains a LaTeX starter manual generated from the uploaded GEOS-MPM
source archive.

## Regenerate narrative + generated appendices

```bash
python3 make_manual.py
pdflatex -interaction=nonstopmode geos_mpm_manual.tex
makeindex geos_mpm_manual
pdflatex -interaction=nonstopmode geos_mpm_manual.tex
pdflatex -interaction=nonstopmode geos_mpm_manual.tex
```

## Refresh code-derived data

`generated/geos_mpm_extracted.json` was produced from the uploaded source tree.
`tools/update_generated_tables.py` is the extraction script used for this run and
should be pointed at or adapted for the current checkout before regenerating.

The `linked_reports/` directory contains the verification, validation, and
examples report PDFs found in the archive, plus a placeholder for future
LLNL-specific material-model documentation.
