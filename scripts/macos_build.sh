#!/bin/bash
# Run from the repository root in the macOS CI job.

set -euo pipefail
brew install cmake open-mpi openblas python@3.14 binutils
curl -fsSL -o "$RUNNER_TEMP/sccache.tgz" \
  https://github.com/mozilla/sccache/releases/download/v0.17.0/sccache-v0.17.0-aarch64-apple-darwin.tar.gz
tar -xzf "$RUNNER_TEMP/sccache.tgz" -C "$RUNNER_TEMP" --strip-components=1 \
  sccache-v0.17.0-aarch64-apple-darwin/sccache
chmod +x "$RUNNER_TEMP/sccache"
xattr -d com.apple.quarantine "$RUNNER_TEMP/sccache" 2>/dev/null || true
# sccache reads SCCACHE_CONF first. On macOS the default is
# ~/Library/Application Support, so a Linux-style ~/.config file is ignored
# and the cache stays on the runner's local disk.
export SCCACHE_CONF="${HOME}/.config/sccache/config"
mkdir -p "$(dirname "$SCCACHE_CONF")"
cat > "$SCCACHE_CONF" <<EOF
basedirs = ["${GITHUB_WORKSPACE}"]
[cache.gcs]
rw_mode = "READ_WRITE"
cred_path = "${GOOGLE_GHA_CREDS_PATH}"
bucket = "geos-dev"
key_prefix = "sccache"
EOF
"$RUNNER_TEMP/sccache" --show-stats || true
"$(brew --prefix python@3.14)/bin/python3.14" -m pip install --upgrade --break-system-packages virtualenv
# Caliper records the Open MPI Cellar path from the TPL build. Homebrew on this
# runner may install a different version, so point that path at the one just installed.
if [[ ! -d /opt/homebrew/Cellar/open-mpi/5.0.11 ]]; then
  sudo mkdir -p /opt/homebrew/Cellar/open-mpi
  sudo ln -sfn "$(readlink -f /opt/homebrew/opt/open-mpi)" /opt/homebrew/Cellar/open-mpi/5.0.11
fi
# Apple Clang rejects -mtune=native, which the host-config adds for every Clang.
sed -i '' 's/ -march=native -mtune=native//' /opt/GEOS/spack-generated.cmake
log="$RUNNER_TEMP/geos-build.log"
set +e
set -o pipefail
{
  python3 scripts/config-build.py \
    -hc /opt/GEOS/spack-generated.cmake \
    -bt Release \
    --buildpath "$RUNNER_TEMP/geos-build" \
    -DCMAKE_C_COMPILER_LAUNCHER="$RUNNER_TEMP/sccache" \
    -DCMAKE_CXX_COMPILER_LAUNCHER="$RUNNER_TEMP/sccache"
  cmake --build "$RUNNER_TEMP/geos-build" -j 2 --target geosx
} 2>&1 | tee "$log"
status=$?
set +o pipefail
mkdir -p "${GITHUB_WORKSPACE}/.sccache-runtime"
"$RUNNER_TEMP/sccache" --show-adv-stats > "${GITHUB_WORKSPACE}/.sccache-runtime/stats.txt" || true
if [[ "$status" -ne 0 ]]; then
  grep -E 'error:|CMake Error' "$log" | tail -n 8 | while IFS= read -r line; do
    echo "::error title=Build GEOS::$line"
  done || true
  exit "$status"
fi
