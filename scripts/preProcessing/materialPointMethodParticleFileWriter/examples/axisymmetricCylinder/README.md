# Axisymmetric cylinder interpolation comparison

This smoke case applies the same small radial contraction pulse to a solid
elastic cylinder: radius 1, length 1, density 1000, initial ring volume pi.
The default uses a CPDI core and CPTI skin. The optional B-spline case uses
quad glyphs and the one-sided cubic basis at the axis. No CZ is present.

Run from an empty run directory with the usual PFW `userDefs` configuration:

```bash
AXIS_CYLINDER_INTERPOLATION=bspline \
  /path/to/materialPointMethodParticleFileWriter/examples/axisymmetricCylinder/runProblem --generate-only
```

Set `GEOS_EXECUTABLE` and omit `--generate-only` to run the generated XML.
Set `AXIS_CYLINDER_UPDATE=FMPM` for FMPM(2); FLIP is the default. Use separate
run directories for each variant. Without `AXIS_CYLINDER_INTERPOLATION`, the
original mixed CPDI/CPTI example is generated.

This is not an analytical dynamics benchmark. Check finite stress/velocity,
particle retention near the axis, hoop strain, and timestep/mesh convergence
on the target build. B-spline smoothing and point quadrature can change the
response relative to CPDI/CPTI. Source-level tests check mapping consistency,
near-axis regularity and hydrostatic force balance; they do not establish
full-timestep, actual MPI, restart or GPU correctness.

The axis closure does not modify Cartesian F-table boundary handling. Its
weights apply to both radial and axial components before mass/force assembly;
it is not the existing post-mapping velocity reflection applied in a new place.
