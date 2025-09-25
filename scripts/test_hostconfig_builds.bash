#!/usr/bin/env bash
# test-builds.sh
# - Usage: ./test-builds.sh -hc <hostconfig_directory> [-f filter]
# - Finds host-configs in the given dir (non-recursive) whose filename includes the filter
# - For each: configure (via config-build.py), build, then install
# - Logs each phase; prints success/failure summary
set -uo pipefail

# --- defaults ---
HOSTCONFIG_DIR=""
FILTER="$(hostname -s)"
LOG_DIR="_build_logs"
BUILD_TYPES=("Debug" "Release" "RelWithDebInfo")
JOBS="$(/usr/bin/getconf _NPROCESSORS_ONLN 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)"
CONFIG_WRAPPER="python3 scripts/config-build.py"
WRAPPER_EXTRA_ARGS="--ninja"

usage() {
  cat <<EOF
Usage: $0 -hc <hostconfig_directory> [-f filter]
  -hc   Directory containing host-config .cmake files (non-recursive) [required]
  -f    Filter string to match in filenames (default: \$(hostname -s) -> $FILTER)

Examples:
  $0 -hc host-configs/LLNL
  $0 -hc host-configs/LLNL -f toss_4
EOF
  exit 1
}

# --- parse flags ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -hc) HOSTCONFIG_DIR="$2"; shift 2 ;;
    -f)  FILTER="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "[ERROR] Unknown option: $1"; usage ;;
  esac
done

[[ -z "$HOSTCONFIG_DIR" ]] && { echo "[ERROR] -hc is required"; usage; }
[[ -d "$HOSTCONFIG_DIR" ]] || { echo "[ERROR] Directory '$HOSTCONFIG_DIR' not found"; exit 1; }

# ------------------------------------------------
ts() { date "+%Y-%m-%d %H:%M:%S"; }
lower() { echo "$1" | tr '[:upper:]' '[:lower:]'; }
mkdir -p "$LOG_DIR" || { echo "[ERROR] Cannot create $LOG_DIR"; exit 1; }

# 1) discover host-configs (non-recursive, match filter in filename)
mapfile -t HOSTCONFIGS < <(find "$HOSTCONFIG_DIR" -maxdepth 1 -type f -name "*.cmake" \
  | awk -v hn="$FILTER" -F/ '{ if ($NF ~ hn) print }' | sort)

if ((${#HOSTCONFIGS[@]} == 0)); then
  echo "[ERROR] No *.cmake host-configs in '$HOSTCONFIG_DIR' matching filter '$FILTER'"
  exit 1
fi

echo "[INFO] $(ts) Found ${#HOSTCONFIGS[@]} host-config(s) matching '$FILTER':"
printf '  - %s\n' "${HOSTCONFIGS[@]}"

declare -A RESULT
SUCCESS_COUNT=0
FAIL_COUNT=0

run_one() {
  local hc="$1" bt="$2"
  local hc_file="$(basename "$hc")"
  local hc_base="${hc_file%.cmake}"
  local bt_lc; bt_lc="$(lower "$bt")"

  # wrapper default: ./build-<hc_base>-<bt_lc>
  local build_dir="./build-${hc_base}-${bt_lc}"

  local log_cfg="${LOG_DIR}/${hc_base}__${bt}__configure.log"
  local log_bld="${LOG_DIR}/${hc_base}__${bt}__build.log"
  local log_ins="${LOG_DIR}/${hc_base}__${bt}__install.log"

  echo
  echo "[INFO] $(ts) === ${hc_base} :: ${bt} ==="
  echo "[INFO] Expected build dir: $build_dir"
  echo "[INFO] Logs: $log_cfg, $log_bld, $log_ins"

  # 2) configure (wrapper creates/cleans build dir)
  local cfg_cmd="$CONFIG_WRAPPER -hc \"$hc\" -bt \"$bt\" $WRAPPER_EXTRA_ARGS"
  echo "[INFO] Running config: $cfg_cmd"
  eval $cfg_cmd >"$log_cfg" 2>&1
  local cfg_rc=$?
  if (( cfg_rc != 0 )); then
    echo "[FAIL] $(ts) Configure failed (rc=$cfg_rc). See $log_cfg"
    return $cfg_rc
  fi

  # 3) build
  local bld_cmd="cmake --build \"$build_dir\" --parallel \"$JOBS\""
  echo "[INFO] Running build: $bld_cmd"
  eval $bld_cmd >"$log_bld" 2>&1
  local bld_rc=$?
  if (( bld_rc != 0 )); then
    echo "[FAIL] $(ts) Build failed (rc=$bld_rc). See $log_bld"
    return $bld_rc
  fi

  # 4) install
  local ins_cmd="cmake --build \"$build_dir\" --target install"
  echo "[INFO] Running install: $ins_cmd"
  eval $ins_cmd >"$log_ins" 2>&1
  local ins_rc=$?
  if (( ins_rc != 0 )); then
    echo "[FAIL] $(ts) Install failed (rc=$ins_rc). See $log_ins"
    return $ins_rc
  fi

  echo "[OK]   $(ts) Build + Install succeeded"
  return 0
}

# main loop
for hc in "${HOSTCONFIGS[@]}"; do
  for bt in "${BUILD_TYPES[@]}"; do
    key="$(basename "${hc%.cmake}")__${bt}"
    if run_one "$hc" "$bt"; then
      RESULT["$key"]="OK";   ((SUCCESS_COUNT++))
    else
      RESULT["$key"]="FAIL"; ((FAIL_COUNT++))
    fi
  done
done

# summary
echo
echo "================ Build Summary ================"
printf "%-60s  %s\n" "Config (host-config __ buildtype)" "Result"
printf -- "------------------------------------------------------------  ------\n"
for k in "${!RESULT[@]}"; do
  printf "%-60s  %s\n" "$k" "${RESULT[$k]}"
done
echo "------------------------------------------------------------  ------"
echo "Successes: $SUCCESS_COUNT"
echo "Failures : $FAIL_COUNT"
echo "Logs in  : $LOG_DIR"
echo "Builds in: (current directory)"
echo "==============================================="
(( FAIL_COUNT == 0 )) || exit 1
