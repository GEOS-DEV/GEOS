# Axisymmetric pressure regions

The default case is a hollow elastic cylinder with inner radius a=10 mm,
outer radius b=20 mm and length H=40 mm. Axial symmetry planes impose plane
strain; radial grid faces are outflow, well outside the material. Inner pressure
ramps to 100 kPa and outer pressure to 20 kPa over 5 ms, then holds to 7 ms.
E=100 MPa, nu=0.25, density=2000 kg/m3. CPDI quads resolve the meridional plane.
Pressures are positive in compression; stresses are negative in compression.

## Generate and run

From a new run directory, with normal PFW userDefs and GEOS_EXECUTABLE set:

```bash
/path/to/materialPointMethodParticleFileWriter/examples/axisymmetricPressureRegions/runProblem --generate-only
# Or generate, run directly, and report history values:
/path/to/materialPointMethodParticleFileWriter/examples/axisymmetricPressureRegions/runProblem
```

The runner does not submit a scheduler job. Use the generated XML and particle
files with the usual site launcher if needed. Set AXIS_PRESSURE_CASE=sphere
before generation and postprocessing to select the spherical example.

## Which nodes receive pressure?

- BoreholePressure: r < 10.5 mm, independent of axial y. The cutoff lies just
  outside the inner material radius so inner-surface grid support is selected.
- ConfiningPressure box: outside [-1,19.5 mm] x [-1,1] m in (r,y). The lower
  radial and both axial bounds are beyond all grid nodes, including buffers,
  so only r > 19.5 mm receives the outer pressure. The three-entry box XML
  retains z=0 placeholders, ignored in axisymmetry.
- AXIS_PRESSURE_CASE=sphere: material is the positive-r half of a radius-20-mm
  meridional disk. Pressure ramps to 20 kPa where r*r+y*y > (19.5 mm)^2. This
  describes external pressure around an origin-centered sphere after revolution.

Nodes exactly on a cutoff are not selected. Cutoffs are fixed spatial values
tested against current grid coordinates. For moving/deforming surfaces choose
cutoffs that continue to include the intended boundary support; they do not
track the surface. Avoid placing a cutoff across a surface that should remain
unloaded. For the hollow cylinder choose a < r0 < r1 < b, with sufficient
separation relative to interpolation support and expected deformation.

Confining pressure overrides borehole pressure in an overlap; it does not add.
The solver stores one region of each kind. Later events of the same kind update
that shared region and its held pressure. Events require useEvents=1; completed
axisymmetric ramps hold their exact endPressure, including after restart.

## Pressure, traction and interpretation

These events use virtual-fluid background stress -p I, including hoop stress.
Particle volumes already include 2*pi*r. No surface flags, normals, particle
traction records, extra circumference factor or thickness are needed.
This is a spatial hydrostatic-pressure mask, not a named-surface load: any
free/interface boundary represented in the selected nodal support can be loaded.
The radial bands also intersect the cylinder end planes, where axial symmetry
supplies the plane-strain reaction. For finite open-ended cylinders, check
whether those end-band loads match the intended boundary condition.

Generic PFW per-area surface-traction integration remains unsupported in
axisymmetry. Direct GEOS particleSurfaceTraction values represent already
integrated full-ring forces, not pressure or force per radian.

## Static reference and validation

For the hollow cylinder, the small-strain Lame plane-strain solution is

```
A = (pi*a*a - po*b*b)/(b*b-a*a)
B = a*a*b*b*(pi-po)/(b*b-a*a)
sigma_rr(r)   = A - B/(r*r)
sigma_hoop(r) = A + B/(r*r)
sigma_yy     = 2*nu*A
u_r(r) = (1+nu)/E * ((1-2*nu)*A*r + B/r)
```

Here pi and po denote inner and outer pressure, not the circle constant.
postProcess.py compares box-history material-average stresses with these
ring-volume-weighted references. For the sphere the reference is sigma=-po I.
Neither comparison is an automatic pass/fail: transient oscillations, cutoff
placement, curved-boundary particle resolution and finite strains affect it.

Local validation compiles actual production pressure selectors, force assembly,
input validators and event handlers using CPU adapters with strict warnings and
address/undefined-behavior sanitizers. It checks mixed CPDI/CPTI inner/outer
cylindrical resultants, hoop balance, pressure work, region precedence, endpoint
holds and Cartesian regressions. Complete writer runs use 1/3 emulated ranks
and compare particle bytes and XML. A full GEOS build, dynamics benchmark,
actual MPI and restart roundtrip remain target-machine checks.
