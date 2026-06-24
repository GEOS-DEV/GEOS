#!/usr/bin/env bash
set -euo pipefail

stamp="$(date +%y%m%d%H%M)"

exec "$HOME/geos/scripts/make_mpm_review_bundle.sh" \
  --lean \
  -o "$HOME/geos/${stamp}_geos_mpm_bundle.zip"
