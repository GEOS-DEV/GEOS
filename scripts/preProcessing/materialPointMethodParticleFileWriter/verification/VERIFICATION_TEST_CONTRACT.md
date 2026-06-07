# GEOS-MPM verification test-folder contract

Modern verification tests are folder scoped.  A folder owns the physics setup,
the run variants, post-processing logic, report text, and figure assets needed to
build its section of the verification-suite report.

Required files for a modern folder:

- `pfw_input_<test>.py` - one clean, single-purpose PFW input file.  It should
  contain only non-default settings and should express mesh size through
  `refine`, `cpp`, and `ppc` variables when practical.  Series are defined by a
  small number of environment variables rather than by duplicating input files.
- `runTest` - executable Python script that calls `mpm_vv_folder_runner.run_folder`.
  The script defines the variants, case names, environment overrides, and the
  folder default case id.
- `postProcess_<test>.py` - executable Python post-processor.  It should be safe
  to run even if one variant fails, write logs/summary files, and produce the
  quantitative metrics used in the report.
- `test.tex` - report section fragment.  It should describe the setup, model and
  feature under test, expected/analytical solution, error metric, and generated
  results include.
- `setup_*.png` and any other figure assets needed by `test.tex`.

Recommended generated post-processing outputs:

- `<test>_metrics.csv` with one row per diagnostic sample.
- `<test>_summary.json` with variant-level pass/fail-relevant metrics.
- `<test>_results.tex` for inclusion by `test.tex`.
- Plot PNGs for histories/errors.
- Optional VisIt PNGs under `visit_frames/`, using a consistent variable-specific
  color-map convention across the suite.

Source inputs should not include the temporary fast-debug override blocks that
were formerly appended by the suite harness.  Fast staging changes belong in the
suite runner or in the staged working copy only, not in the canonical input.

The first migrated legacy folders following this contract are:

- `Ftable` - elastic F-table boundary switch.
- `stressControl` - elastic stress-control tracking variants.
- `vonMisesJ` - uniaxial Von Mises plasticity.
- `periodicBoundaries` - two-dimensional periodic free advection.
- `diskThruPartitions` - 3 by 3 partition-crossing free flight.
- `spinningDisk` - PIC/FLIP/FMPM2 angular-momentum and energy drift.

Legacy inputs retained in migrated folders are renamed `legacy_pfw_input_*.py` so
automatic discovery does not treat them as active verification cases.  They can
be used as references when future split tests are developed.
