#!/bin/bash
# Run the eigensolver comparison of the documentation example. Each run writes a log.
# usage: run_comparison.sh GEOSX GEOS_REPOSITORY OUTDIR
GEOSX=${1:?path to geosx}
REPO=${2:?path to the GEOS repository}
OUT=${3:?output directory}
DECKS=$REPO/inputFiles/solidMechanics
mkdir -p "$OUT"
cp "$DECKS"/modalFreeBlock_*.xml "$DECKS"/modalFreeFreeBeam_*.xml "$OUT"/
cd "$OUT"

run() {  # name deck sed-expression
  sed -e "$3" "$2" > "$1.xml"
  start=$(date +%s.%N)
  "$GEOSX" -i "$1.xml" > "$1.log" 2>&1
  end=$(date +%s.%N)
  echo "$1 $(echo "$end - $start" | bc)" >> wall.txt
}
: > wall.txt

# Free block, 16 modes
run block_arnoldi              modalFreeBlock_arnoldi.xml 's/XX//'
run block_arnoldi_deflated     modalFreeBlock_arnoldi.xml 's/modalDeflateRigidBodyModes="0"/modalDeflateRigidBodyModes="1"/'
run block_lobpcg               modalFreeBlock_lobpcg.xml  's/modalDeflateRigidBodyModes="1"/modalDeflateRigidBodyModes="0"/'
run block_lobpcg_deflated      modalFreeBlock_lobpcg.xml  's/XX//'

# Free-free beam, 10 and 20 modes
beam() {  # name solver deflate modes
  run "$1" modalFreeFreeBeam_arnoldi_smoke.xml \
      "s/modalSolverType=\"arnoldi\"/modalSolverType=\"$2\" modalDeflateRigidBodyModes=\"$3\" modalMaxIterations=\"400\" logLevel=\"1\"/;s/modalNumModes=\"10\"/modalNumModes=\"$4\"/"
}
beam beam10_arnoldi          arnoldi 0 10
beam beam10_arnoldi_deflated arnoldi 1 10
beam beam10_lobpcg           lobpcg  0 10
beam beam10_lobpcg_deflated  lobpcg  1 10
beam beam20_arnoldi          arnoldi 0 20
beam beam20_lobpcg_deflated  lobpcg  1 20
