# weakTraceCouetteShear

No-slip platen-fluid Couette shear benchmark for the prescribed-surface C0
vertex-mortar treatment.

The folder dispatches three variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact groups tied by the C0 mortar projection.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.

The domain is plane strain and periodic only in the flow direction (`x`). The
geometry is a three-layer stack: a lower platen, a viscous liquid layer, and an
upper platen. The top and bottom boundaries receive opposite tangential
velocities. The `traceContactGroups` variant places both platen-fluid surfaces
in the C0 mortar pair `{0, 1}` so the test exercises no-slip transfer across
the weak discontinuity instead of tying two elastic layers.

The diagnostic goal is to compare the final fluid shear rate and stress against
the single-field baseline. In the baseline, platen and fluid momentum should be
smeared across the mixed cells near both interfaces. With the C0 interface, the
fluid velocity gradient should be closer to the Couette value
`2 * plate_velocity / fluid_gap`, and the Liquid-model shear stress should be
approximately `fluid_viscosity * plate_velocity / fluid_gap` with the current
GEOS symmetric-gradient convention.

The input requests particle material type, surface flag, owning particle rank,
velocity, stress, and deformation gradient in the Silo output so the interface
location, MPI partitioning, and fluid shear strain can be inspected in VisIt.
The single-field variant still writes the weak-surface markers so surface
refinement and diagnostics use the same particles as the C0 variant; weak
projection and contact remain disabled for that baseline.

## Parallel periodic C0 run

The default particle type is cubic `SinglePointBSpline` and the default periodic
layout has three MPI partitions in `x` and three in `y`. A representative XPIC2
C0 run is:

```bash
WEAK_TRACE_PARTICLE_TYPE=1 \
WEAK_TRACE_UPDATE_METHOD=XPIC \
WEAK_TRACE_UPDATE_ORDER=2 \
WEAK_TRACE_XPAR=3 \
WEAK_TRACE_YPAR=3 \
WEAK_TRACE_NI=120 \
WEAK_TRACE_NJ=120 \
WEAK_TRACE_PPC=4 \
WEAK_TRACE_SURFACE_REFINEMENT_LEVEL=4 \
WEAK_TRACE_PROJECTION_ITERATIONS=30 \
WEAK_TRACE_PROJECTION_SCALE=0.6 \
WEAK_TRACE_GAP_STABILIZATION=0.0 \
WEAK_TRACE_MAX_GAP_CORRECTION_VELOCITY=0.0 \
./runTest --variant singleField --variant traceContactGroups --force --jobs 1
```

`WEAK_TRACE_XPAR` and `WEAK_TRACE_YPAR` set the MPI decomposition, while
`WEAK_TRACE_NI` and `WEAK_TRACE_NJ` can hold the global mesh fixed across
layouts. The global cell counts must be divisible by the corresponding
partition counts. Periodic C0-mortar runs require `WEAK_TRACE_XPAR >= 2`; cubic B-spline runs
have the same restriction in the ordinary MPM mapping. Three or more partitions
are recommended so the periodic halo does not rely on a same-rank wrap.

For a decomposition check, keep `WEAK_TRACE_NI=120` and `WEAK_TRACE_NJ=120`
fixed and run, for example, `3x1`, `4x1`, and `4x2`. Compare the velocity
profile, integrated reaction, global reconstructed segment/constraint counts,
and the logged post-projection velocity jump.

## Current MPI implementation scope

The plane-strain C0 mortar implementation keeps particles and bulk grid fields
distributed. It all-gathers the particle-carried interface markers, reconstructs
each background cell on one deterministic owning rank, all-gathers the resulting
codimension-one segments, and replicates only the mortar support state. Every
rank then executes the same deterministic forward/reverse projection sweeps and
applies the resulting correction to its local owner and ghost copies.

Periodically equivalent grid nodes and interface crossings are identified by
canonical global structured-grid keys, while distance fits and segment lengths
use minimum-image displacements. This joins the horizontal interface across the
periodic seam instead of treating the two lateral ends as free endpoints.

This first correctness-oriented path is CPU-only, plane strain, and requires one
MPI partition in `z`. Its communication and storage scale with the number of
surface markers and mortar-support nodes because those codimension-one data are
replicated. A later 3-D implementation should replace the replication with
neighbor marker exchange and a distributed interface solve.
