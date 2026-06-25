# Isolated GEOS-MPM material-point driver

This package provides a non-invasive material-point control layer for single-particle or single-element constitutive tests. It lives entirely under the PFW scripts tree and does not add a branch to `SolidMechanicsMPM` or change any MPM solver code.

The Python layer handles:

- material selection through `pfw_materials.py` or inline material dictionaries;
- material direction input as a normal/fiber vector or a GEOS/PFW-style row-wise `3 x 3` frame;
- material direction updates with `fixed`, `rotation`, `fiber`, `normal`, `graphite`, or `mpmCofactor` kinematics;
- mixed strain/stress control with finite-difference Newton iterations over stress-controlled components;
- state restoration for trial evaluations by cloning the point state before every residual evaluation;
- internal-energy and accumulated stress-power integration using `sigma : D / rho`; and
- CSV output suitable for fitting and optimization workflows.

The built-in `elastic` backend is for smoke tests and verification of the driver mechanics. Production GEOS-MPM constitutive testing should use the optional compiled executable `geos_mpm_material_point_driver`, generated from `src/coreComponents/constitutive/materialPointDriver` with `GEOS_ENABLE_MPM_MATERIAL_POINT_DRIVER=ON`. With `setupMPM`, this remains default-off and can be requested with `./setupMPM --dane --enable-material-point-driver` or `./setupMPM --tuolumne --enable-material-point-driver`.

Typical Python smoke test:

```bash
cd scripts/preProcessing/materialPointMethodParticleFileWriter
python3 pfw_material_point.py --write-example /tmp/material_point_case.json
python3 pfw_material_point.py --case /tmp/material_point_case.json --output /tmp/material_point_history.csv
```

Write sidecar files for the compiled GEOS driver:

```bash
cd scripts/preProcessing/materialPointMethodParticleFileWriter
python3 pfw_material_point.py \
  --case /tmp/material_point_case.json \
  --export-compiled-driver-prefix /tmp/material_point_case

# After building the optional C++ target:
/tmp/material_point_case.run.sh
```

Kroonblawd/Figure-1-style Graphite input generation:

```bash
cd scripts/preProcessing/materialPointMethodParticleFileWriter
python3 verification/materialPoint/pfw_input_graphite_kroonblawd_material_point.py \
  --work-dir /tmp/graphite_theta75_p30 \
  --theta-deg 75 \
  --pressure-gpa 30
```

The Kroonblawd-style example sets lab `z` as the compression direction, uses the first material-frame row as the graphite basal-plane normal, keeps the run isothermal at 300 K, and accumulates stress work while retaining zero thermal energy by default.
