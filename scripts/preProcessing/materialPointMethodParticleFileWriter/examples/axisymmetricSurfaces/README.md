# Axisymmetric surface calculations

Run from a fresh output directory with your usual PFW userDefs:

```sh
export GEOS_EXECUTABLE=/path/to/geosx
/path/to/materialPointMethodParticleFileWriter/examples/axisymmetricSurfaces/runProblem
```

Use `--generate-only` to stop after PFW. The runner records `geos.log` and the
postprocessor reports particle volume and first/last globally summed field-0
area. `runProblem` belongs beside the other PFW examples and does not submit a
scheduler job. Use your normal launcher if GEOS requires srun/mpirun.

Environment options:

| Variable | Values/default |
|---|---|
| AXIS_SURFACE_CASE | cylinder (default), sphere |
| AXIS_SURFACE_INTERPOLATION | mixed (default), cpdi, cpti, bspline |
| AXIS_SURFACE_TENSION | 1 (default), 0 for geometry only |
| AXIS_SURFACE_RESOLUTION | 24 (default), cells per unit length; increase to 48 |

Cylinder: R=H=1, volume pi, lateral area 2*pi, gamma=.01, initial pressure .01.
Axial symmetry planes exclude end caps from the free surface. The mixed case
uses CPDI core and CPTI skin, with zero input normals and offsets so the solver
must reconstruct them. Only the outer cylindrical face is flagged.

Sphere: R=1, volume 4*pi/3, area 4*pi, initial pressure .02. Mixed mode uses
CPDI bulk with cut CPTI surface triangles. The cutter uses analytical normals
in PFW; GEOS overwrites them in its first reconstruction pass.

Both use a liquid with density 1, bulk modulus 10 and zero viscosity. Initial
Laplace pressure is a continuum reference, not an exact discrete equilibrium.
Expect residual waves; compare resolutions and inspect gridSurfaceTensionForce,
gridSurfaceArea, reconstructed particle normals/positions and velocity. The
surface-tension capillary timestep restriction remains active. The example is
a smoke/refinement diagnostic, not a validated equilibrium benchmark.

The sphere approximation depends on particle refinement. `postProcess.py`
checks initial volume, prints the initial area error and requires finite
positive area when a GEOS log is present; it does not certify dynamics.
