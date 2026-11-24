#!/bin/bash
# Run from the repository root in the macOS CI job.

set -euo pipefail
log="$RUNNER_TEMP/load-tpl.log"
set +e
{
  curl -fsSL -o "$RUNNER_TEMP/crane.tgz" \
    https://github.com/google/go-containerregistry/releases/download/v0.20.6/go-containerregistry_Darwin_arm64.tar.gz
  tar -xzf "$RUNNER_TEMP/crane.tgz" -C "$RUNNER_TEMP" crane
  chmod +x "$RUNNER_TEMP/crane"
  xattr -d com.apple.quarantine "$RUNNER_TEMP/crane" 2>/dev/null || true
  image="geosx/macos-26-apple-clang:${DOCKER_IMAGE_TAG}"
  "$RUNNER_TEMP/crane" export "$image" "$RUNNER_TEMP/tpl.tar"
  sudo mkdir -p /opt/GEOS
  sudo /usr/bin/python3 scripts/macos_extract_tpl.py "$RUNNER_TEMP/tpl.tar" /opt/GEOS
  test -f /opt/GEOS/spack-generated.cmake
  sudo chown -R "$(id -un)" /opt/GEOS
} >"$log" 2>&1
status=$?
set -e
if [[ "$status" -ne 0 ]]; then
  tail -n 40 "$log" | while IFS= read -r line; do
    echo "::error title=Load TPL image::$line"
  done
  exit "$status"
fi
