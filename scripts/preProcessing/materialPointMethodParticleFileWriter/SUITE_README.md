# GEOS-MPM pfw suites

This directory contains a stdlib-only suite runner for the Material Point Method particle file writer (pfw).  It is designed for the current `feature/mpm` branch and the Dane MPM minimal-TPL build.

## Build GEOS on Dane with the MPM minimal-TPL host config

From the GEOS repository root on Dane:

```bash
python3 scripts/config-build.py \
  -hc host-configs/LLNL/dane-toss_4_x86_64_ib-gcc@12.1.1-mpm-minimal-tpl.cmake \
  -bt Release \
  -br /usr/workspace/$USER/geos-builds \
  -ir /usr/workspace/$USER/geos-installs

cmake --build /usr/workspace/$USER/geos-builds/build-dane-toss_4_x86_64_ib-gcc@12.1.1-mpm-minimal-tpl-release -j
```

Set `GEOS_EXECUTABLE` to the resulting executable before running suites, for example:

```bash
export GEOS_EXECUTABLE=/usr/workspace/$USER/geos-builds/build-dane-toss_4_x86_64_ib-gcc@12.1.1-mpm-minimal-tpl-release/bin/geosx
export GEOS_BANK=<your_lc_bank>
```

Adjust the build and install roots to your workspace policy.  The suite runner writes `userDefs_$USER.py` into each prepared run directory, so no machine-local `userDefs` file needs to be committed.

## Parent-suite commands

The parent wrappers are:

```bash
examples/runExamplesSuite
verification/runVerificationSuite
validation/runValidationSuite
```

Common commands:

```bash
# Static compile/deprecated-keyword check only.
./verification/runVerificationSuite preflight

# List cases that will be discovered.  Template inputs with XXXX/YYYY/ZZZZ are skipped by default.
./verification/runVerificationSuite list

# Prepare run directories and invoke pfw, but do not submit GEOS jobs.
./verification/runVerificationSuite run \
  --run-root /p/lustre1/$USER/geos-mpm-suite-runs \
  --geos-path "$GEOS_EXECUTABLE" \
  --bank "$GEOS_BANK" \
  --partition pdebug \
  --output-type silo \
  --keep-going

# Prepare and submit generated GEOS batch scripts.
./verification/runVerificationSuite run \
  --run-root /p/lustre1/$USER/geos-mpm-suite-runs \
  --geos-path "$GEOS_EXECUTABLE" \
  --bank "$GEOS_BANK" \
  --partition pdebug \
  --output-type silo \
  --submit \
  --auto-restart \
  --keep-going

# Refresh status/report after jobs have run.
./verification/runVerificationSuite status \
  --run-root /p/lustre1/$USER/geos-mpm-suite-runs \
  --geos-path "$GEOS_EXECUTABLE"
```

Reports are written by default to:

```text
<run-root>/<suite-name>/suite_report.md
<run-root>/<suite-name>/suite_report.json
```

The report includes case IDs, source paths, status, generated XML/particle files, scheduler job IDs when pfw reports them, presence of `reactionHistory.csv` and `boxAverageHistory.csv`, and the leading comment-block description from each input file.

## Per-problem commands

Each problem directory containing `pfw_input_*.py` now has a `runProblem` wrapper.  For example:

```bash
cd verification/Ftable
./runProblem preflight
./runProblem run --run-root /p/lustre1/$USER/geos-mpm-suite-runs --geos-path "$GEOS_EXECUTABLE" --bank "$GEOS_BANK" --submit
./runProblem status --run-root /p/lustre1/$USER/geos-mpm-suite-runs
```

Per-problem wrappers discover only the direct `pfw_input_*.py` files in that directory.  Parent-suite wrappers discover recursively.

## Case selection

Use fnmatch-style case filters:

```bash
./verification/runVerificationSuite run --only 'Ftable/*' --only 'stressControl/*'
./verification/runVerificationSuite run --exclude 'expandingRing/*' --exclude 'expandingBar/*'
```

Template inputs containing `XXXX`, `YYYY`, or `ZZZZ` are skipped by default.  Add `--include-templates` only when a template has first been concretized or the script intentionally handles parameter substitution.

## VisIt rendering

`pfw_visit_render.py` can render initial/final frames after a case has Silo/VTK output:

```bash
visit -nowin -cli -s pfw_visit_render.py \
  --run-dir /p/lustre1/$USER/geos-mpm-suite-runs/verification/Ftable/elasticBlockUni \
  --variable Damage \
  --view xy
```

Frames are written to `<run-dir>/visit_frames`.  If the requested variable is not present, the script tries several common MPM variables (`Damage`, `Dmg`, `MaterialType`, `ParticleType`, `Density`, `Velocity`).

## Files added for the suite workflow

- `pfw_suite.py`: suite discovery, run preparation, pfw invocation, status scan, Markdown/JSON report generation, and run script emission.
- `pfw_visit_render.py`: VisIt CLI renderer for initial/final frames.
- `examples/runExamplesSuite`, `verification/runVerificationSuite`, `validation/runValidationSuite`: parent-suite wrappers.
- `runProblem` in each problem directory containing direct `pfw_input_*.py` files.

## Modernization checks applied in this branch copy

- Replaced old geometry keyword `v=` with the current `vel=` keyword in pfw inputs.
- Fixed three syntax-level input issues that blocked compile preflight.
- Fixed a `particleFileWriter.py` typo where NumPy array normalization called `_normalize_ndarray` instead of `normalize_ndarray`.
