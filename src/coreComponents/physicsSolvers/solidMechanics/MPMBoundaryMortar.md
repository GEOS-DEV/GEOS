# Boundary-compatible C0 mortar: opt-in CPU implementation

Performance update (2026-09-25): see [MPMBoundaryMortarPerformance.md](MPMBoundaryMortarPerformance.md)
for native/LAPACK QR backends and concurrent component workers. The mathematics
and acceptance criteria below are unchanged. Historical test limitations below
refer to the original boundary-compatible patch, not a new full GEOS validation.

## Baseline and status

This is an incremental patch against the exact contents of
`2609231646_geos_mpm_bundle.zip`, whose metadata identifies commit
`46b1dc01edcdda855091cc97bce795555b427a2b` with the supplied worktree overlay.
The two-pass PFW and native reference-geometry additions are ALREADY in that
baseline. Do not reapply their earlier patches.

The production algebra helper and particle writer have been tested separately.
The GEOS adapter has been inspected but has NOT been compiled in a complete GEOS
build. No GEOS mechanics run, actual MPI transport run, GPU run, or restart test
has been performed for this patch. Treat the included octant input as a test to
run, not a verified mechanics result.

## Enable the new path

The spherical test input supplies these options. For another PFW input using
prescribed groups and explicit weak surfaces:

```python
pfw["enableWeakInterfaceTraceProjection"] = 1
pfw["weakInterfaceTraceBoundaryCoupling"] = 1
pfw["weakInterfaceTraceProjectionScale"] = 1.0
pfw["weakInterfaceTraceSuppressNodalContact"] = 1
pfw["weakInterfaceTraceLinearSolverTolerance"] = 1.0e-10
pfw["weakInterfaceTraceAbsoluteTolerance"] = 1.0e-12
pfw["weakInterfaceTraceQRRankTolerance"] = 1.0e-12
pfw["weakInterfaceTraceQRMaxEntries"] = 16000000
```

The existing trace configuration still applies: single-point cubic B-spline,
explicit surface geometry, prescribed contact-group fields, `damageFieldPartitioning=0`,
selected `weakInterfaceTracePairs`, and the existing interpolation restrictions.

`weakInterfaceTraceBoundaryCoupling=0` is the default and leaves the original
GS/fallback and BC ordering enabled. With coupling enabled, the GS iteration cap
and denominator regularization are NOT used: the projection is not a relaxed
sweep or a ridge-regularized tie. A rank threshold is not an interface compliance.

The new options are passed through by PFW without any new geometry queries or
changes to rank-sliced particle generation. Old inputs omit them. Rebuild GEOS
before using generated XML containing the new attributes. Regenerate an XSD
from that executable if your workflow validates input against a saved schema.

## Supported and deliberately unsupported combinations

Supported scope is the existing grid-aligned six-face geometry, with arbitrary
orientation of the intersecting material interface. Symmetry and bilateral
Moving/FTable faces are handled. A Moving face needs the existing FTableBC or
stress-control kinematics. Optional prescribed transverse components retain the
existing component precedence at edges/corners. Free tangential components are
not made dependent on reflected interior values.

The code explicitly rejects unilateral Contact faces and an already-applied
legacy material-contact/CZ impulse on participating nodes/fields. It does not
silently claim to preserve a contact law applied before the mortar correction.
A coupled active-set/friction treatment is future work. Nodal contact for the
selected bonded pair is already suppressed by the existing trace-pair option.

Other boundaries of this patch:

- CPU reference solve only; not a distributed or GPU linear-solver replacement.
- Prescribed explicit-surface groups only; no automatic residual/LR or DFG path.
- No new CPDI/CPTI interpolation; existing R-vectors, native particle geometry,
  mass, volume, density, material state, and cohesive update formulas are unchanged.
- No general immersed or arbitrarily rotated EXTERNAL boundary surface. The
  mortar surface itself may intersect the six grid-aligned faces at any angle.
- `prescribedBoundaryFTable` is supported; superimposed `prescribedFTable`
  particle deformation is rejected in this mode.
- Existing higher-order FMPM trace restriction remains. PIC is the octant test
  default. XPIC2 is selectable, but its filtered particle motion and evolving
  interface geometry still need complete-step verification.
- Restricting integration to the physical grid box does not solve the separate
  problem of identifying an arbitrary free-surface termination inside that box.

## Geometry: physical faces, legitimate edges, symmetry images

The new branch uses the physical boundary planes at `m_xGlobalMin/Max`, one cell
inside each nonperiodic padded grid face. The background-cell tetrahedra already
align with those planes, so skipping exterior cells is an exact restriction;
there is no new polygon mesher or integration refinement. The same interior
mortar vertices/functions and quadrature rules are retained.

A reconstructed open edge on a physical face is accepted as a legitimate
termination. An incident exterior facet cannot invalidate an interior row,
because exterior facets are not assembled in this branch. Physical buffer nodes
are nevertheless retained in the B-spline trace.

Near a symmetry face, marker points and normals are reflected for the nodal
surface fit. Combinations handle intersections of two or three symmetry planes.
These are fit samples only: no mirrored particle mass, volume, or exterior
integration facets are created. Moving faces do not imply reflected geometry.
Symmetry still requires physically mirrored material geometry and properties;
it is not a way to assert symmetry of a non-symmetric specimen.

## Algebra: preserve boundary dependencies in the mortar correction

For a normal Moving or symmetry condition, the existing GEOS relation is

```
v_boundary = U
v_buffer   = -v_interior + 2 U
```

Consequently, an allowed increment obeys

```
delta_v_boundary = 0
delta_v_buffer   = -delta_v_interior
```

The helper represents raw velocities as `v = v_D + Z q`, using a signed-affine
dependency map. It includes both the donor and buffer mapped masses in
`M_reduced = Z^T M Z`. Dependencies are formed separately for each material
field and Cartesian component. No nodal common velocity is imposed merely
because two materials touch the same symmetry plane.

The implementation preserves the original GEOS essential-BC update, including
its acceleration/dVelocity semantics, then finds the minimum-mass correction
compatible with those BCs and the accepted mortar rows. This is a deliberately
bounded adapter change: it is not a claim that the original BC operation itself
was replaced by a full minimum-mass projection from the unconstrained trial.

Mass/support checks occur after reducing dependencies. A zero-mass buffer can
contribute through its mass-carrying donor instead of automatically losing its
basis coefficient. Common support filtering remains part of the method. A row
that requests necessary nodal fallback is removed, its ties are included in the
dependency map, and the remaining rows are reduced again before the solve.
There is no post-solve nodal fallback overwrite.

For a reduced row operator D, solve the row-scaled rectangular system

```
B y = b,  minimizing ||y||,  B = S D M_reduced^(-1/2).
```

All backends use column-pivoted Householder QR of `B^T`, rather than explicitly
forming the poorly conditioned Gram matrix. It computes a minimum-norm
correction for the retained numerical rank, then checks ALL accepted original
rows, including dependent ones. Incompatible or insufficiently resolved
constraints cause an error, not a silent least-squares compromise.

Acceptance requires absolute-plus-relative RMS AND maximum average-jump tests,
plus a boundary-dependency residual check. The absolute tolerance has velocity
units. A small initial residual may satisfy the absolute test without satisfying
an arbitrarily tiny relative ratio; inspect both. No fallback to GS occurs if QR
fails, exceeds its memory cap, or does not pass acceptance.

The factor is dense and bounded by `weakInterfaceTraceQRMaxEntries`; 16 million
entries corresponds to 128 million bytes for the factor alone, not total memory.
Factorization time also grows rapidly. Start with the small included problem;
a sparse/distributed solver is needed for a large interface network.

## Timestep staging and FTable

In the opt-in branch:

1. Retain the existing explicit grid dynamics/contact stage.
2. Evaluate table/controller kinematics and apply the existing essential BCs once.
3. Reconstruct/restrict the mortar; close boundary dependencies; gather state.
4. Solve compatible corrections; update velocity, momentum, and active-field
   acceleration/dVelocity consistently.
5. Skip the later duplicate table/BC application, but retain profile output.
6. After owner-to-ghost synchronization, verify that the solved grid velocities
   were not changed before G2P.

This audit is a grid-stage check, not a proof about later XPIC/FMPM filtering or
particle surface-marker transport. The original diagonal FTable endpoint rate
and exact finite grid-stretch ratio are retained. Table/controller history is
not advanced inside QR or separately for each Cartesian component.

## MPI ownership and limits

The existing gathered canonical-key interface representation is retained.
Required nodes are closed under both donor-to-buffer and buffer-to-donor
relations before the nodal-state gather, including edge/corner combinations.
Owner-preferred deduplication counts a mass once, rather than summing partition
copies. Physical buffer nodes are distinct from duplicate MPI ghosts.

Each scalar component is a complete dense solve. Up to three ranks now solve
x/y/z concurrently; `weakInterfaceTraceComponentRanks=1` restores rank-zero
scheduling. Every component status is exchanged before any result is committed.
Every rank receives the same solution and reconstructs local entries;
ranks with no local interface still participate in collectives. A post-sync
check is collective even on an empty rank. Missing donor records and the MPI
count limit fail explicitly. Actual MPI transport and GEOS field migration are
not established by the standalone/emulated tests supplied with this patch.

## Diagnostics and reaction accounting

The logs identify geometry/fallback counts, component/time/cycle, accepted rows,
rank/free columns, QR pivot scale, final relative/absolute jumps, boundary
residual, and fallback degrees of freedom. Use the MPMSubroutines log category
or the corresponding solver log level. Do not infer convergence from a run
merely reaching its final time.

The existing `gridWeakInterfaceTraceForce` represents the TOTAL compatible
velocity-correction force in this mode. The new
`gridWeakInterfaceTraceBoundaryForce` is the ADDITIONAL BC/internal-fallback
reaction needed to balance the raw mortar impulse. Fixed components can have a
reaction even when their velocity correction is zero. With rank redundancy,
individual dual contributions need not be unique.

The old `reactionHistory` still describes the preceding legacy BC update; it is
not automatically augmented with this extra reaction. Inspect/sum the new force
on owned nodes as appropriate. Do not label the old history as the complete
coupled wall reaction. No momentum loss should be diagnosed without including
external boundary impulses, and moving-wall work must be separated from internal
projection energy removal.

Fallback diagnostics distinguish a deliberate tie from an accepted mortar row;
a zero jump at a fallback node is not evidence that its original row survived.
The complete-step audit checks solution preservation after synchronization, not
only the QR library residual.

## Files and tests

- `kernels/MPMBoundaryMortar.hpp`: dependency-free production algebra and grid-face helpers.
- `SolidMechanicsMPM.cpp/.hpp`: options, geometry restriction, MPI adapter, staging and diagnostics.
- `tests/testMPMBoundaryMortar.cpp`: standalone production-header unit tests.
- `tests/testBoundaryMortarReference.py`: independent SVD and ordering comparisons.
- PFW `verification/weakTraceSphericalInclusionOctant`: paired pressure inputs,
  separate FTable smoke, static 3D sphere reference and log/initialization checks.

From the GEOS root (C++14 compiler, NumPy and SciPy for the independent reference):

```bash
python src/coreComponents/physicsSolvers/solidMechanics/tests/testBoundaryMortarReference.py \
  --output /tmp/boundary-mortar-algebra
```

For the spherical input and run procedure, see its dedicated README. Its positive
octant pressure problem has three true symmetry planes and three positive-side
void gaps. The FTable cube variant is a different physical problem; it is not
compared to the spherical pressure reference.

Required integration acceptance remains: a complete GEOS CPU build, serial
pressure run, serial/2x2x2 comparison, restart, full/mirrored-domain comparison,
final particle kinematics, and stress/traction convergence near all three
symmetry intersections. None is replaced by the standalone matrix tests.
