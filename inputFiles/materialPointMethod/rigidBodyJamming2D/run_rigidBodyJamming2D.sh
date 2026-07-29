#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
cd "$SCRIPT_DIR"

GEOSX_EXECUTABLE=${GEOSX_EXECUTABLE:-geosx}

python3 generateRigidBodyJamming2DParticles.py

if [ -n "${MPI_LAUNCHER:-}" ]; then
  # shellcheck disable=SC2206
  launcher=( ${MPI_LAUNCHER} )
elif command -v srun >/dev/null 2>&1; then
  launcher=(srun -n 6)
elif command -v mpirun >/dev/null 2>&1; then
  launcher=(mpirun -n 6)
else
  echo "rigidBodyJamming2D requires an MPI launcher; set MPI_LAUNCHER explicitly." >&2
  exit 2
fi

"${launcher[@]}" "$GEOSX_EXECUTABLE" \
  -i mpm_rigidBodyJamming2D.xml \
  -x 2 -y 3 -z 1

python3 checkRigidBodyJamming2D.py --require-run-output
