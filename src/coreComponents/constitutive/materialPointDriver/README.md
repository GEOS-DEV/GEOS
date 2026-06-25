# Optional GEOS-MPM material-point driver

This directory contains a standalone executable target for single-point MPM
constitutive tests.  It is intentionally isolated from the production MPM solver:

- it is not built unless `GEOS_ENABLE_MPM_MATERIAL_POINT_DRIVER=ON`;
- it links to the existing `constitutive` target and dispatches through
  `ConstitutivePassThruMPM< ContinuumBase >`;
- it does not modify `SolidMechanicsMPM`, particle/grid kernels, MPM events, or
  the default MPM CMake targets.

The driver input is deliberately simple.  A Python front end in
`scripts/preProcessing/materialPointMethodParticleFileWriter/pfw_material_point.py`
imports material dictionaries from `pfw_materials.py`, writes a GEOS
`<Constitutive>` XML fragment, writes a load-path CSV, and launches this
executable.

Minimal CMake usage:

```bash
cmake -DGEOS_ENABLE_MPM_MATERIAL_POINT_DRIVER=ON ...
cmake --build <build-dir> --target geos_mpm_material_point_driver
```

With the bundled `setupMPM` helper, the driver is still default-off.  Enable and build it explicitly with:

```bash
./setupMPM --dane --enable-material-point-driver
./setupMPM --tuolumne --enable-material-point-driver
```

`setupMPM` forwards `-DGEOS_ENABLE_MPM_MATERIAL_POINT_DRIVER=ON` only when this option is requested, and otherwise forwards `OFF` so a default MPM setup does not build the driver.

The load-path CSV has columns

```text
dt,xx,yy,zz,yz,xz,xy
```

and the six component values are interpreted according to the six modes passed
through `--control`.  Supported modes are `strain`, `strainRate`,
`trueStrainRate`, and `stress`.
