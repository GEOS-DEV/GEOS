#!/usr/bin/env bash
# Flux batch smoke test for the GEOS-MPM colliding-disks Tuolumne case.
#
# Submit from the GEOS repo root with:
#   flux batch tuolemneSmokeTest/runSmoke.sh
#
# Or override Flux scheduler options on submission:
#   flux batch -B mahem -q pdebug -t 30m tuolemneSmokeTest/runSmoke.sh
#
# Useful runtime overrides:
#   GEOSX=/path/to/geosx flux batch tuolemneSmokeTest/runSmoke.sh
#   REPO_ROOT=/usr/WS1/$USER/GEOS flux batch tuolemneSmokeTest/runSmoke.sh
#   CASE_DIR=/usr/WS1/$USER/GEOS/tuolemneSmokeTest flux batch tuolemneSmokeTest/runSmoke.sh
#   RUN_DIR=/p/lustre5/$USER/geosxTests/my_run flux batch tuolemneSmokeTest/runSmoke.sh
#   TRACE_DATA_MIGRATION=1 flux batch tuolemneSmokeTest/runSmoke.sh
#   RANKS=1 XPAR=1 YPAR=1 ZPAR=1 flux batch tuolemneSmokeTest/runSmoke.sh

#flux: -N 1
#flux: -q pdebug
#flux: -t 30m
#flux: -B mahem

set -Eeuo pipefail

msg() {
  printf '[runSmoke_tuolemne] %s\n' "$*"
}

quote_cmd() {
  printf '%q ' "$@"
  printf '\n'
}

run_cmd() {
  msg "running: $(quote_cmd "$@")"
  "$@"
}

abs_dir_if_exists() {
  local path="$1"
  if [ -d "$path" ]; then
    (cd "$path" && pwd -P)
  else
    return 1
  fi
}

has_case_files() {
  local dir="$1"
  [ -f "$dir/$INPUT_XML_NAME" ] && [ -f "$dir/$PARTICLE_FILE_NAME" ]
}

find_case_dir() {
  local candidates=()
  local cwd script_dir lc_user
  cwd=$(pwd -P)
  script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
  lc_user=${USER:-$(id -un)}

  if [ -n "${CASE_DIR:-}" ]; then
    candidates+=("$CASE_DIR")
  fi

  if [ -n "${REPO_ROOT:-}" ]; then
    candidates+=("$REPO_ROOT/tuolemneSmokeTest")
    candidates+=("$REPO_ROOT/tuolumneSmokeTest")
    candidates+=("$REPO_ROOT")
  fi

  # When run under flux batch, BASH_SOURCE often points at a temporary copy of
  # the script, but PWD is usually the submission directory. Check PWD first.
  candidates+=("$cwd/tuolemneSmokeTest")
  candidates+=("$cwd/tuolumneSmokeTest")
  candidates+=("$cwd")

  # This works when the script is run directly rather than copied by flux batch.
  candidates+=("$script_dir")
  candidates+=("$script_dir/tuolemneSmokeTest")
  candidates+=("$script_dir/tuolumneSmokeTest")

  # Common LLNL checkout locations used in this MPM workflow.
  candidates+=("$HOME/geosx/tuolemneSmokeTest")
  candidates+=("$HOME/geosx/tuolumneSmokeTest")
  candidates+=("/usr/WS1/${lc_user}/GEOS/tuolemneSmokeTest")
  candidates+=("/usr/WS1/${lc_user}/GEOS/tuolumneSmokeTest")
  candidates+=("/g/g19/${lc_user}/geosx/tuolemneSmokeTest")
  candidates+=("/g/g19/${lc_user}/geosx/tuolumneSmokeTest")

  local dir real seen=""
  for dir in "${candidates[@]}"; do
    [ -n "$dir" ] || continue
    real=$(abs_dir_if_exists "$dir" 2>/dev/null || true)
    [ -n "$real" ] || continue
    case " $seen " in
      *" $real "*) continue ;;
    esac
    seen="$seen $real"
    if has_case_files "$real"; then
      printf '%s\n' "$real"
      return 0
    fi
  done

  printf 'Searched these candidate directories for %s and %s:\n' "$INPUT_XML_NAME" "$PARTICLE_FILE_NAME" >&2
  for dir in "${candidates[@]}"; do
    printf '  %s\n' "$dir" >&2
  done
  return 1
}

CASE_NAME=${CASE_NAME:-collidingDisks}
INPUT_XML_NAME=${INPUT_XML_NAME:-mpm_collidingDisks.xml}
PARTICLE_FILE_NAME=${PARTICLE_FILE_NAME:-mpmParticleFile_collidingDisks}

LC_USER=${USER:-$(id -un)}
TIMESTAMP=$(date +%Y%m%d-%H%M%S)
JOB_ID=${FLUX_JOB_ID:-manual}
RUN_ROOT=${RUN_ROOT:-/p/lustre5/${LC_USER}/geosxTests}
RUN_DIR=${RUN_DIR:-${RUN_ROOT}/tuolemneSmokeTest_${CASE_NAME}_${TIMESTAMP}_${JOB_ID}}
LOG_FILE=${LOG_FILE:-${RUN_DIR}/runSmoke_${CASE_NAME}_${TIMESTAMP}.log}
GEOS_LOG=${GEOS_LOG:-${RUN_DIR}/geosx_${CASE_NAME}_${TIMESTAMP}.log}
ERRORS_FILE=${ERRORS_FILE:-${RUN_DIR}/geosx_errors_${CASE_NAME}_${TIMESTAMP}.yaml}
OUTPUT_DIR=${OUTPUT_DIR:-${RUN_DIR}/geosx_output}

GEOSX=${GEOSX:-/g/g19/${LC_USER}/geosx/build-tuolumne-toss_4_x86_64_ib_cray-llvm-amdgpu@6.4.2-rocm@6.4.2-mpm-minimal-tpl-release/bin/geosx}
NODES=${NODES:-1}
RANKS=${RANKS:-4}
XPAR=${XPAR:-2}
YPAR=${YPAR:-2}
ZPAR=${ZPAR:-1}
OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
TRACE_DATA_MIGRATION=${TRACE_DATA_MIGRATION:-0}
MPIBIND_VERBOSE=${MPIBIND_VERBOSE:-1}

mkdir -p "$RUN_DIR" "$OUTPUT_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

finish() {
  local rc=$?
  msg "finished with exit code: ${rc}"
  msg "main log: ${LOG_FILE}"
  msg "GEOS log: ${GEOS_LOG}"
  msg "errors file: ${ERRORS_FILE}"
  msg "run directory: ${RUN_DIR}"
  if [ -f "$GEOS_LOG" ]; then
    msg "last 120 lines of GEOS log:"
    tail -120 "$GEOS_LOG" || true
  fi
  exit "$rc"
}
trap finish EXIT

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
SUBMIT_OR_RUN_PWD=$(pwd -P)
CASE_SOURCE_DIR=$(find_case_dir)
REPO_ROOT_DETECTED=$(cd "$CASE_SOURCE_DIR/.." && pwd -P)

msg "starting Tuolumne colliding-disks smoke test"
msg "date: $(date)"
msg "host: $(hostname)"
msg "script dir: ${SCRIPT_DIR}"
msg "submit/run PWD: ${SUBMIT_OR_RUN_PWD}"
msg "case source dir: ${CASE_SOURCE_DIR}"
msg "repo root detected: ${REPO_ROOT_DETECTED}"
msg "run dir: ${RUN_DIR}"
msg "GEOSX: ${GEOSX}"

if [ ! -x "$GEOSX" ]; then
  msg "ERROR: GEOSX executable not found or not executable: ${GEOSX}"
  msg "Set GEOSX=/path/to/geosx and resubmit."
  exit 2
fi

# Always stage the XML and particle file into the run directory. GEOS then runs
# from the staged directory, so the ParticleMesh particleFile path is local and
# independent of where Flux copied this batch script.
run_cmd cp -p "$CASE_SOURCE_DIR/$INPUT_XML_NAME" "$RUN_DIR/$INPUT_XML_NAME"
run_cmd cp -p "$CASE_SOURCE_DIR/$PARTICLE_FILE_NAME" "$RUN_DIR/$PARTICLE_FILE_NAME"
run_cmd ls -lh "$RUN_DIR/$INPUT_XML_NAME" "$RUN_DIR/$PARTICLE_FILE_NAME"

cd "$RUN_DIR"

# Runtime environment for Tuolumne GPU-aware MPI/HIP runs.
export MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-1}
export HSA_XNACK=${HSA_XNACK:-1}
export OMP_NUM_THREADS
export GEOS_MPM_MACHINE=${GEOS_MPM_MACHINE:-tuolumne}
ulimit -c unlimited || true

msg "resource/launch settings: NODES=${NODES}, RANKS=${RANKS}, XPAR=${XPAR}, YPAR=${YPAR}, ZPAR=${ZPAR}, OMP_NUM_THREADS=${OMP_NUM_THREADS}"
if [ $(( XPAR * YPAR * ZPAR )) -ne "$RANKS" ]; then
  msg "ERROR: x/y/z partition product must equal RANKS. Got $(( XPAR * YPAR * ZPAR )) vs ${RANKS}."
  exit 2
fi

{
  echo '===== basic system info ====='
  date
  hostname
  uname -a
  echo
  echo '===== paths ====='
  echo "script dir: ${SCRIPT_DIR}"
  echo "submit/run PWD: ${SUBMIT_OR_RUN_PWD}"
  echo "case source dir: ${CASE_SOURCE_DIR}"
  echo "repo root detected: ${REPO_ROOT_DETECTED}"
  echo "run dir: ${RUN_DIR}"
  echo
  echo '===== selected environment ====='
  env | sort | grep -E '^(AMDGPU|CHAI|FI_|FLUX|GEOS|HIP|HSA|LD_LIBRARY_PATH|MPICH|OMP|PATH|ROCM|ROCR|UMPIRE|USER)=' || true
  echo
  echo '===== modules ====='
  type module >/dev/null 2>&1 && module list 2>&1 || true
  echo
  echo '===== git state ====='
  git -C "$REPO_ROOT_DETECTED" rev-parse --short HEAD 2>/dev/null || true
  git -C "$REPO_ROOT_DETECTED" status --short 2>/dev/null || true
  git -C "$REPO_ROOT_DETECTED" diff --stat 2>/dev/null || true
  echo
  echo '===== GEOS executable ====='
  ls -lh "$GEOSX" || true
  ldd "$GEOSX" 2>/dev/null | head -200 || true
  echo
  echo '===== GEOS help probe ====='
  "$GEOSX" --help 2>&1 | head -160 || true
  echo
  echo '===== Flux info ====='
  command -v flux || true
  flux version 2>/dev/null || true
  flux resource list 2>/dev/null || true
  flux jobs -a 2>/dev/null | head -80 || true
  echo
  echo '===== ROCm/HIP info ====='
  command -v hipconfig >/dev/null 2>&1 && hipconfig --full 2>&1 | head -240 || true
  command -v rocm-smi >/dev/null 2>&1 && rocm-smi 2>&1 | head -160 || true
  command -v rocminfo >/dev/null 2>&1 && rocminfo 2>&1 | grep -E 'Name:|Marketing Name:|gfx|HSA Agents' | head -160 || true
  echo
  echo '===== staged input summary ====='
  ls -lh
  echo
  grep -n 'SolidMechanics_MPM\|damageFieldPartitioning\|needsNeighborList\|initialDt\|maxTime\|ParticleMesh\|particleFile' "$INPUT_XML_NAME" || true
  echo
  echo '===== particle file header/sample ====='
  head -40 "$PARTICLE_FILE_NAME" || true
} | tee "$RUN_DIR/debug_context_${TIMESTAMP}.txt"

FLUX_CMD=( flux run --env=* --exclusive -N "$NODES" -n "$RANKS" )
if [ "$MPIBIND_VERBOSE" = "1" ]; then
  FLUX_CMD+=( --setopt=mpibind=verbose:1 )
fi

GEOS_ARGS=(
  -i "$INPUT_XML_NAME"
  -x "$XPAR"
  -y "$YPAR"
  -z "$ZPAR"
  -n "${CASE_NAME}_tuolumne_smoke"
  -o "$OUTPUT_DIR"
  -e "$ERRORS_FILE"
)

if [ "$TRACE_DATA_MIGRATION" = "1" ]; then
  GEOS_ARGS+=( --trace-data-migration )
fi

if [ -n "${GEOS_EXTRA_ARGS:-}" ]; then
  # shellcheck disable=SC2206
  EXTRA_ARGS=( ${GEOS_EXTRA_ARGS} )
  GEOS_ARGS+=( "${EXTRA_ARGS[@]}" )
fi

msg "GEOS command: $(quote_cmd "${FLUX_CMD[@]}" "$GEOSX" "${GEOS_ARGS[@]}")"
printf '%q ' "${FLUX_CMD[@]}" "$GEOSX" "${GEOS_ARGS[@]}" > "$RUN_DIR/geosx_command_${TIMESTAMP}.sh"
printf '\n' >> "$RUN_DIR/geosx_command_${TIMESTAMP}.sh"
chmod +x "$RUN_DIR/geosx_command_${TIMESTAMP}.sh"

set +e
if command -v stdbuf >/dev/null 2>&1; then
  stdbuf -oL -eL "${FLUX_CMD[@]}" "$GEOSX" "${GEOS_ARGS[@]}" 2>&1 | tee "$GEOS_LOG"
  GEOS_RC=${PIPESTATUS[0]}
else
  "${FLUX_CMD[@]}" "$GEOSX" "${GEOS_ARGS[@]}" 2>&1 | tee "$GEOS_LOG"
  GEOS_RC=${PIPESTATUS[0]}
fi
set -e

msg "GEOS exit code: ${GEOS_RC}"

msg "post-run file summary:"
find "$RUN_DIR" -maxdepth 2 -type f -printf '%s %p\n' | sort -n | tail -100 || true

if [ -f "$ERRORS_FILE" ]; then
  msg "GEOS errors file contents:"
  sed -n '1,240p' "$ERRORS_FILE" || true
fi

if grep -qiE 'ERROR|Segmentation fault|Signal no\.|StackTrace|device-side assert|hipError|HSA_STATUS|roc|MPI_ABORT|aborting' "$GEOS_LOG"; then
  msg "diagnostic matches found in GEOS log:"
  grep -niE 'ERROR|Segmentation fault|Signal no\.|StackTrace|device-side assert|hipError|HSA_STATUS|roc|MPI_ABORT|aborting' "$GEOS_LOG" | tail -120 || true
fi

exit "$GEOS_RC"
