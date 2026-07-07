#!/usr/bin/env bash
#SBATCH -J sipLayerCompare
#SBATCH -A mahem
#SBATCH -p pdebug
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 00:45:00

set -euo pipefail

if [[ -n "${SURFACE_POLYMER_SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR=$(cd "$SURFACE_POLYMER_SCRIPT_DIR" && pwd)
elif [[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR=$(cd "$SLURM_SUBMIT_DIR" && pwd)
else
  SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
fi

if [[ ! -f "$SCRIPT_DIR/pfw_input_surfaceInformedPolymerLayer.py" ]]; then
  echo "Could not locate surfaceInformedPolymerLayer sources in: $SCRIPT_DIR" >&2
  echo "Run from the surfaceInformedPolymerLayer directory or set SURFACE_POLYMER_SCRIPT_DIR." >&2
  exit 2
fi

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  cd "$SCRIPT_DIR"
  exec sbatch --export=ALL,SURFACE_POLYMER_SCRIPT_DIR="$SCRIPT_DIR" "$SCRIPT_DIR/runTestAndBuildReport.sh" "$@"
fi

PYTHON=${PFW_PYTHON:-${PYTHON:-/usr/tce/bin/python3}}
export PFW_PYTHON="$PYTHON"

SRC=$SCRIPT_DIR
PFW=$(cd "$SRC/../.." && pwd)
RUNROOT=${SURFACE_POLYMER_RUNROOT:-/p/lustre1/homel1/geosxTests/verification}
OUT=$SRC/output/surfaceInformedPolymerLayer
VARIANT_WALLTIME=${SURFACE_POLYMER_VARIANT_WALLTIME:-00:45:00}
STAMP=$(date +%Y%m%d_%H%M%S)
LOGROOT=$SRC/output/surfaceInformedPolymerLayer_clean_logs/$STAMP

rm -rf "$OUT" \
  "$SRC/output/surfaceInformedPolymerLayer__continuum_tension" \
  "$SRC/output/surfaceInformedPolymerLayer__continuum_compression" \
  "$SRC/output/surfaceInformedPolymerLayer__cohesive_tension" \
  "$SRC/output/surfaceInformedPolymerLayer__cohesive_compression" \
  "$RUNROOT/surfaceInformedPolymerLayer__continuum_tension" \
  "$RUNROOT/surfaceInformedPolymerLayer__continuum_compression" \
  "$RUNROOT/surfaceInformedPolymerLayer__cohesive_tension" \
  "$RUNROOT/surfaceInformedPolymerLayer__cohesive_compression"

mkdir -p "$OUT" "$LOGROOT"
exec > >(tee "$LOGROOT/all.log") 2>&1

echo "logs: $LOGROOT"
echo "started: $(date)"
echo "slurm job: $SLURM_JOB_ID"

declare -a VARIANTS RUN_DIRS RUN_SCRIPTS

stage_variant() {
  local variant=$1 label=$2 case_name=$3 model=$4 loading=$5 strain=$6
  local cid="surfaceInformedPolymerLayer__${variant}"
  local rundir="$RUNROOT/$cid"

  echo "staging $variant"
  (
    cd "$SRC"
    SURFACE_POLYMER_CASE_NAME="$case_name" \
    SURFACE_POLYMER_VARIANT_LABEL="$label" \
    SURFACE_POLYMER_MODEL_KIND="$model" \
    SURFACE_POLYMER_LOADING="$loading" \
    SURFACE_POLYMER_FINAL_GLOBAL_STRAIN="$strain" \
    "$PYTHON" "$PFW/mpm_vv_case_runner.py" \
      --suite verification \
      --case-id "$cid" \
      --input pfw_input_surfaceInformedPolymerLayer.py \
      --source-dir "$SRC" \
      --output-prefix "$case_name" \
      --python "$PYTHON" \
      --force \
      --no-submit \
      --no-post \
      --walltime "$VARIANT_WALLTIME"
  ) 2>&1 | tee "$LOGROOT/${variant}_stage.log"

  mkdir -p "$rundir/siloFiles/data"

  shopt -s nullglob
  local scripts=("$rundir"/*_runGEOS.sh)
  shopt -u nullglob
  if [[ ${#scripts[@]} -ne 1 ]]; then
    echo "expected one runGEOS script for $variant, found ${#scripts[@]}" >&2
    exit 2
  fi

  sed -i 's/^srun /srun --exclusive /' "${scripts[0]}"

  VARIANTS+=("$variant")
  RUN_DIRS+=("$rundir")
  RUN_SCRIPTS+=("${scripts[0]}")
}

stage_variant continuum_tension "Continuum tension" surfaceInformedPolymerLayer_continuum_tension continuum tension 0.10
stage_variant continuum_compression "Continuum compression" surfaceInformedPolymerLayer_continuum_compression continuum compression 0.08
stage_variant cohesive_tension "Cohesive-zone tension" surfaceInformedPolymerLayer_cohesive_tension cohesive tension 0.10
stage_variant cohesive_compression "Cohesive-zone compression" surfaceInformedPolymerLayer_cohesive_compression cohesive compression 0.08

pids=()
for i in "${!VARIANTS[@]}"; do
  variant=${VARIANTS[$i]}
  rundir=${RUN_DIRS[$i]}
  run_script=${RUN_SCRIPTS[$i]}
  echo "launching $variant"
  (
    cd "$rundir"
    bash "$run_script"
  ) > "$LOGROOT/${variant}_geos.log" 2>&1 &
  pids+=("$!")
done

failed=0
for i in "${!pids[@]}"; do
  if wait "${pids[$i]}"; then
    echo "${VARIANTS[$i]} complete"
  else
    rc=$?
    echo "${VARIANTS[$i]} failed with rc=$rc; see $LOGROOT/${VARIANTS[$i]}_geos.log"
    failed=1
  fi
done
[[ $failed -eq 0 ]]

"$PYTHON" - "$SRC" "$OUT" <<'PY'
import json
import os
import sys
from pathlib import Path

src = Path(sys.argv[1])
out = Path(sys.argv[2])
variants = [
    ("continuum_tension", "Continuum tension", "surfaceInformedPolymerLayer_continuum_tension"),
    ("continuum_compression", "Continuum compression", "surfaceInformedPolymerLayer_continuum_compression"),
    ("cohesive_tension", "Cohesive-zone tension", "surfaceInformedPolymerLayer_cohesive_tension"),
    ("cohesive_compression", "Cohesive-zone compression", "surfaceInformedPolymerLayer_cohesive_compression"),
]
subcases = []
for name, label, case_name in variants:
    cid = f"surfaceInformedPolymerLayer__{name}"
    jobs_path = src / "output" / cid / f"{case_name}_jobs.json"
    jobs = json.loads(jobs_path.read_text()) if jobs_path.is_file() else {}
    jobs["geos_job_id"] = os.environ.get("SLURM_JOB_ID", "")
    jobs["status"] = "completed-single-allocation"
    subcases.append(
        {
            "name": name,
            "label": label,
            "case_id": cid,
            "case_name": case_name,
            "returncode": 0,
            "jobs": jobs,
        }
    )
manifest = {
    "case_id": "surfaceInformedPolymerLayer",
    "suite": "verification",
    "source_dir": str(src),
    "input": "pfw_input_surfaceInformedPolymerLayer.py",
    "single_allocation_job_id": os.environ.get("SLURM_JOB_ID", ""),
    "subcases": subcases,
}
out.mkdir(parents=True, exist_ok=True)
(out / "surfaceInformedPolymerLayer_jobs.json").write_text(json.dumps(manifest, indent=2))
PY

mkdir -p "$OUT/mplconfig" "$OUT/tmp"

echo "postprocessing"
MPLCONFIGDIR="$OUT/mplconfig" TMPDIR="$OUT/tmp" \
"$PYTHON" "$SRC/postProcess_surfaceInformedPolymerLayer.py" \
  --suite verification \
  --source-dir "$SRC" \
  --output-dir "$OUT" \
  --case-id surfaceInformedPolymerLayer \
  --python "$PYTHON" \
  --no-visit \
  2>&1 | tee "$LOGROOT/post.log"

echo "building report"
cd "$PFW"
verification/runSuite report \
  --case surfaceInformedPolymerLayer \
  --python "$PYTHON" \
  2>&1 | tee "$LOGROOT/report.log"

echo "finished: $(date)"
echo "output: $OUT"
echo "report: $PFW/verification/verification_suite_report/verification_suite_report.pdf"
echo "logs: $LOGROOT"
