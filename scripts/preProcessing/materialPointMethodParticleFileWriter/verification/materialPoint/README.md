# Optional material-point driver verification

This directory contains verification cases for the isolated material-point driver
`geos_mpm_material_point_driver`.  The driver is intentionally independent of the
production MPM solver and is built only when the optional CMake/setupMPM switch is
enabled.

Build the driver, for example:

```bash
./setupMPM --dane --enable-material-point-driver
```

Run the local Von Mises uniaxial-stress verification:

```bash
cd scripts/preProcessing/materialPointMethodParticleFileWriter/verification/materialPoint
./runTest
```

By default the test uses `GEOS_MPM_MATERIAL_POINT_DRIVER` when it is set.
Otherwise it reads the current PFW `userDefs_<user>.py` file and uses an
executable named `geos_mpm_material_point_driver` next to `userDefs.geosPath`
(the configured `geosx` executable).  You can still override the executable
explicitly:

```bash
python3 pfw_input_vonMisesJ_uniaxial_stress_material_point.py \
  --driver /path/to/geos_mpm_material_point_driver
```

The Von Mises test prescribes axial strain increments in `xx` and solves for
strain increments in `yy`, `zz`, `yz`, `xz`, and `xy` to keep those stress
components at zero.  The output CSV from the compiled driver is compared with the
small-strain uniaxial elastic-perfectly-plastic Von Mises solution.

Run the Triaxial-style Graphite pressure/orientation sweep:

```bash
./runGraphiteTriaxialTest
```

The default Graphite sweep runs five material orientations,
`theta = 0, 30, 45, 60, 90` degrees, each at the seven confining pressures
`P0 = 0, 1, 3, 5, 10, 20, 30` GPa.  It initializes the stress hydrostatically to
`[-P0, -P0, -P0, 0, 0, 0]`, holds the out-of-plane `xx` strain at zero, controls
the expanding periodic `yy` direction to `-P0`, applies axial true-strain-rate
control in `zz`, and keeps the shear strain components fixed at zero.  The plots
use the compression-positive quantity `-sigma_zz - P0` versus accumulated axial
true strain `eps_zz`.

The Graphite sweep is a plotting/regression-data generator rather than a strict
pass/fail unit test.  By default it asks the compiled driver to continue after
stress-control nonconvergence, writing a full curve with non-converged rows
flagged in the CSV and summary.  The driver also has a scalar bracket/bisection
fallback for one-component stress control, which is the default plane-strain
Graphite protocol here.  Use `--driver-stress-failure-policy stop` to recover the
older partial-curve behavior, or add `--strict --fail-fast` if you want the
wrapper to return a nonzero exit code on the first failed case.  The Python
wrapper captures the compiled driver's stdout/stderr to each case's `.driver.log`
file by default, which keeps LvArray/CHAI allocation and free messages out of the
terminal while preserving them for debugging.

The Graphite sweep writes one processed CSV, one JSON summary, one combined
five-panel plot, and one per-orientation plot:

```text
graphite_triaxial_sweep_material_point/graphite_triaxial_sweep_processed.csv
graphite_triaxial_sweep_material_point/graphite_triaxial_sweep_summary.json
graphite_triaxial_sweep_material_point/graphite_triaxial_response_5panel.png
graphite_triaxial_sweep_material_point/graphite_triaxial_response_theta_000.png
...
```

For a quick input-generation check without running 35 material-point solves:

```bash
./runGraphiteTriaxialTest --prepare-only
```
