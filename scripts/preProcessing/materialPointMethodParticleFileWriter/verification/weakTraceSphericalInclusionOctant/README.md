# Spherical inclusion: one octant and three symmetry planes

This is the positive octant (one eighth of the full 3D body), not a 2D quadrant.
The planes x=0, y=0, and z=0 use the existing GEOS Symmetry option. The inclusion
and outer sphere are centered at the origin, so this is a genuine mirror-symmetric
geometry. The C0 triangular interface need not align with these planes.

## Pressure problem (default)

- Inclusion radius a=0.22; specimen radius b=0.72.
- Matrix/inclusion Young moduli 1 and 4; common nu=0.25; density 1.
- Virtual-fluid `ConfiningPressure` ramps to p=0.004 over 8 time units and holds;
  stop time 12. This uses the origin-centered
  `confiningPressureSphereRadius` mode, with the default radius halfway between
  the inclusion and specimen radii. The pressure background is therefore kept
  away from the weak material interface while still loading the outer curved
  surface region; no Moving box face is used.
- The three positive computational faces are Outflow and are separated from the
  specimen by three empty physical cells by default. They are not pressure walls.
- Total grid cells per direction: 24, including the pre-existing one-cell padding
  at both ends. L=b/(1-gap/(N-2)), so L-b=gap*h exactly. The negative physical
  limits are zero, while the negative padded grid nodes lie at -h.
- Two bulk particles per cell direction; two-pass simplex refinement level 2
  subdivides only particles whose regenerated surface query returns a surface
  flag. Native tetrahedron output is used for cut particles, boxes elsewhere.
  Cartesian surface refinement remains selectable for diagnostics, but the
  default avoids bulk Cartesian particle refinement.
- The benchmark narrows surface tagging to
  `min(PFW skinDepth, WEAK_TRACE_SURFACE_SKIN_CELLS*h)` with default
  `WEAK_TRACE_SURFACE_SKIN_CELLS=0.5`. This keeps the weak/interface markers
  close to the actual spherical surfaces instead of flagging several ppc layers.
- PIC is the reference-test default. XPIC2 can be selected after establishing the
  baseline; this patch does not add a separate projection within XPIC's filters.

Generate separate run directories so one variant cannot overwrite the other:

```bash
# Use your normal PFW userDefs_<username>.py/configuration in each run directory.
# PFW and GEOS paths below are examples to be replaced by your repository paths.
PFW=/path/to/geos/scripts/preProcessing/materialPointMethodParticleFileWriter
cp "$PFW/verification/weakTraceSphericalInclusionOctant/pfw_input_weakTraceSphericalInclusionOctant.py" pfw_input_case.py

WEAK_TRACE_VARIANT=traceContactGroups python "$PFW/particleFileWriter.py" pfw_input_case.py
# Then use your rebuilt CPU executable/generated launch script:
/path/to/geos_executable -i mpm_case.xml -x 1 -y 1 -z 1
```

Repeat in a different directory with `WEAK_TRACE_VARIANT=singleField`.
Both variants retain the same geometry, flags, positions, reference volumes,
and refinement; only contact-group labels and solver selection differ.
The generated input's `weakInterfaceTraceBoundaryCoupling=1` enables the patch.

To reproduce the earlier coarse bulk-refined run, set
`WEAK_TRACE_CPP=14 WEAK_TRACE_PPC=2 WEAK_TRACE_SURFACE_SKIN_CELLS=1e9`.
To use Cartesian surface refinement, set
`WEAK_TRACE_SURFACE_REFINEMENT_METHOD=cartesian`.
For XPIC2 use `WEAK_TRACE_UPDATE_METHOD=XPIC WEAK_TRACE_UPDATE_ORDER=2`.
For a smaller initial diagnostic, use `WEAK_TRACE_CPP=10 WEAK_TRACE_PPC=2`.
To move the virtual-fluid transition surface, set
`WEAK_TRACE_CONFINING_SPHERE_RADIUS`; it must stay between the inclusion and
outer specimen radii.
For a batch submission instead of the debug queue, set
`WEAK_TRACE_RUN_DEBUG=0 WEAK_TRACE_PARTITION=pbatch`.
Do not start a large convergence study before validating the small CPU run:
this implementation is a bounded dense QR reference solve, not a scalable solve.

For solver MPI decomposition, set WEAK_TRACE_XPAR, WEAK_TRACE_YPAR, and
WEAK_TRACE_ZPAR when generating inputs. CPP must be divisible by each count.
Use the same -x/-y/-z values on the GEOS command line. PFW's MPI process count
is separate: its original rank-sliced generation remains unchanged.
For example CPP=14 permits 2x2x2 GEOS partitions without changing grid spacing
relative to the serial CPP=14 case. Each local grid must still meet the existing
B-spline support requirements, including required ghost nodes.

## Three-dimensional static reference

```bash
python "$PFW/verification/weakTraceSphericalInclusionOctant/spherical_reference.py" \
  --output spherical_static_reference.csv
```

The script derives the concentric-sphere solution, not the old plane-strain
circular-inclusion formula. With K=E/[3(1-2nu)] and G=E/[2(1+nu)]:

```
u_in(r) = Ai*r
u_matrix(r) = Am*r + B/r^2
sigma_rr,in = sigma_tt,in = 3*Ki*Ai
sigma_rr,matrix = 3*Km*Am - 4*Gm*B/r^3
sigma_tt,matrix = 3*Km*Am + 2*Gm*B/r^3
```

Coefficients satisfy displacement and radial-traction continuity at r=a and
sigma_rr(b)=-p. The CSV has separate inclusion and matrix rows at r=a: hoop
stress is allowed to jump. This is a SMALL-STRAIN STATIC reference. A finite-rate
ramp or an incompletely equilibrated solution need not equal it. Check kinetic
energy, a slower-ramp repeat, and hold-time convergence before reporting errors.
No simulated stress results accompany this patch.

For a particle at radial direction n, evaluate `sigma_rr=n^T sigma n` and the
mean tangential stress `(trace(sigma)-sigma_rr)/2`. Use both material-side traces
and shared plotting limits. The supplied tracers cover directions close to the
symmetry planes/axes and the interior diagonal; a tracer follows its selected
particle, so inspect its recorded material ID and actual position.

## Separate FTable intersection smoke test

```bash
WEAK_TRACE_LOADING=ftableCube WEAK_TRACE_VARIANT=traceContactGroups \
  python "$PFW/particleFileWriter.py" pfw_input_case.py
```

This changes the physical problem to a cube of side 0.72 cut by a spherical
interface of radius 1.2*L. The interface now intersects all three positive Moving
faces as well as all three negative symmetry planes. A diagonal cosine FTable
compresses to 0.998 and reverses to 1, then holds. No transverse velocities are
prescribed and no Contact faces or virtual-fluid pressure are enabled.

This problem intentionally has NO concentric-sphere stress reference. It tests
BC endpoint targets, buffer dependencies, intersections and reversal. Do not
compare its stresses to `spherical_reference.py`.

## Diagnostics and acceptance

Enable the solver's MPMSubroutines log category (or the appropriate GEOS logLevel)
to retain lines beginning `C0 boundary-compatible QR`. They contain component,
time, cycle, rank, row counts, residuals and fallback counts. For example:

```bash
python "$PFW/verification/weakTraceSphericalInclusionOctant/check_log.py" run.log \
  --minimum-time 12 --maximum-jump 1e-10 --maximum-bc-residual 1e-11
```

The solver already fails if its absolute-plus-relative acceptance tests fail.
This script is a post-run completeness check, not proof of stress accuracy.
Use `--maximum-fallback-dofs 0` for a strict no-fallback diagnostic, not as an
assumption that every coarse geometry will produce zero fallback.

Inspect `gridWeakInterfaceTraceVelocityJumpPost`, `gridWeakInterfaceTraceSkipReason`,
`gridWeakInterfaceTraceForce`, and `gridWeakInterfaceTraceBoundaryForce` near all
three intersection curves. The new BoundaryForce is only the ADDITIONAL reaction
associated with the mortar correction (and internal fallback ties); the old
reactionHistory output still describes the preceding legacy BC operation. It is
not, by itself, the complete coupled reaction history. Sum owned nodes only.
The source deliberately keeps the original boundary acceleration/dVelocity
conventions and applies no second BC overwrite after solving.

Acceptance still required on the full GEOS build: serial run, serial vs 2x2x2,
restart through loading, stress/velocity profiles near and away from symmetry
planes, time/grid/refinement sensitivity, and a full mirrored-domain comparison.
The standalone helper and PFW tests do not substitute for these checks.

## Automated initialization tests

```bash
python test_initialization.py --output /tmp/octant-pfw --communicator mpi
# Only when intentionally testing without an MPI implementation:
python test_initialization.py --output /tmp/octant-pfw --communicator emulated
```

The complete writer is exercised at 1, 3 and 7 process counts. Tests check
identical files across process counts, paired single-field/mortar clouds,
positive native/legacy volumes, region flags and three-plane containment.
Emulated collectives are test-only file exchanges, not an MPI transport test.
