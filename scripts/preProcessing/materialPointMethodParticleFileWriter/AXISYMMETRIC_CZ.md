# Axisymmetric cohesive zones

This patch follows the ordinary axisymmetric contact patch and its PFW enum
fix. It adds prescribed interfaces between exactly two body contact groups,
0 and 1, using the shared explicit timestep and existing CZ reference grid.
The particle record is unchanged. The reference-event fix adds an internal
restart field recording whether that event has already applied its reference.

## Supported configuration

- `formulation='Axisymmetric'`, no swirl, native quad/triangle geometry.
- Q1 point, CPDI quads and CPTI triangles, including mixed CPDI/CPTI bodies.
- `UncoupledCohesiveZone`, `CoupledCohesiveZone`, `BicrystalCohesiveZone`.
  These are traction-separation laws. Use material symmetry compatible with
  axisymmetry; arbitrary anisotropic crystal orientations can break it.
- Explicit `SurfaceNormal` and `SurfacePosition`, `SurfaceFlag=3` on the
  interface, and matching `CZTag` / `CohesiveZoneRegion tag`.
- Set `useEvents=1` on the solver so `ReferenceCohesiveZones` activates the
  interface. With the default zero, the event never executes; surface flags
  and tags alone do not establish a bond. Its displacement options
  `Nodal`, `TypeA` and `TypeB` retain their meanings; `computeNormalsAndPositions`
  must be zero. `TypeB` intersects a normal ray with the actual unscaled
  material shape, so it also works for a triangle.
- `preventCZInterpenetration=1` lets ordinary contact supply compression while
  the CZ supplies opening and shear forces. Once all cached CZ nodes of a
  particle fail, its cohesive flag is cleared and ordinary contact applies.

The solver checks the global number of groups at activation. PFW checks each
cohesive particle's group, tag, meridional vectors and required columns. A
single particle has one CZ tag and one cached reference mapping; arbitrary
multiple field pairs or intersecting interfaces are not introduced here.

`PolymerCohesiveZone` and `SurfaceInformedPolymerCohesiveZone` are explicitly
rejected in axisymmetry. Their finite-layer reconstruction needs in-surface
hoop stretch as a constitutive input; multiplying force by ring area does not
supply that strain or its constitutive coupling. That is the next CZ stage.
DFG, active-CZ subdivision, swirl, Weibull/localized hoop fracture and
solver-reconstructed surfaces remain outside this patch.

## One-time reference activation

Each `ReferenceCohesiveZones` event requests initialization once, on its first
active timestep. An internal `referenceApplied` flag records that request and
is saved in restart files. It is not an XML option. A later, separate event
can deliberately request a new reference configuration for the same region.
The existing finite-end behavior is retained: the region is disabled when
the event ends. General event scheduling and pressure/temperature ramps are
unchanged. This correction applies to Cartesian and axisymmetric CZ events.

Previously, the handler identified the first timestep solely by a window
around `startTime` with width `dt`. A growing timestep could admit two calls
and rebuild the reference after separation had begun. A shrinking timestep
could skip this window. The dispatcher already prevents completed events
from running again, but this event remains active until its end time.
Its `hasStarted` flag is also insufficient: dependent events set that flag
when their start time is scheduled, before any configured delay has elapsed.

After installing this source fix, rebuild and rerun a damaged case from time
zero. An old restart does not contain the new activation flag, and a reference
already overwritten by the bug cannot be recovered by this change.

## Reference and current area

At activation the solver fits the existing locally planar meridional interface
from particle surface normals/positions. For each reference node I it computes

    A_I,CZ = integral_Gamma 2*pi*r*N_I ds.

The line is clipped to the Q1 support rectangle and r>=0, then split at the
radial/axial nodal-hat kinks. Two-point Gauss integration on each segment is
exact: Q1 restricted to the line times radius is at most cubic. Both Mesh and
BruteForce selections use this same native-2D integral. No artificial
thickness, radius floor or division by radius is used.

A radial-axis node has finite area for an axial interface: its support extends
to positive radius. A cylindrical interface at r=0 has zero area. Thus using
`2*pi*nodeRadius*planarArea` is incorrect. For the optional region setting
`czVolumeNormalization=1`, the denominator is also ring-weighted:

    V_I = integral_(r>=0) 2*pi*r*N_I dA,
    V_I(r_I=0) = 2*pi*h_r^2*h_y/6.

Away from the axis it is `2*pi*r_I*h_r*h_y`. The numerator remains the existing
mapped particle-volume estimate. Assign the region's tag throughout the
adjacent bodies if using this normalization, not only to the surface layer.
Specify `czVolumeNormalization` on the `CohesiveZoneRegion`; this patch does
not change the existing event-versus-region parameter plumbing.

This is exact integration of a reconstructed local line, not exact clipping
of every physical particle face or a finite curved interface. Finite interface
ends and corners retain the existing occupancy-normalization approximation;
check convergence of their total area and reactions under refinement.

Store the existing particle activation deformation gradient F_CZ. At each
step use the full 3x3 relative deformation and cofactor:

    H = F * inverse(F_CZ),
    a_I = average(cofactor(H)) * (A_I,CZ * n_I,CZ).

The cofactor's meridional block includes hoop stretch. Under diagonal
`H=diag(lambda_r,lambda_y,lambda_theta)`, an axial interface scales by
`lambda_r*lambda_theta`, and a cylindrical one by `lambda_y*lambda_theta`.
The existing law resolves opening and meridional slip, and multiplies traction
by the current area once. Mass and force already represent full rings.

The patch corrects the old absolute-F area update for all formulations. That
update was inconsistent when CZ activation occurred at nonidentity F and when
F was reset. Reference normals are now unit vectors before multiplying area,
and all nine mapped cofactor components are initialized and normalized.

## Resets, damage and outputs

The existing reset relation

    F_CZ,new = F_CZ,old * inverse(F_old) * F_new

preserves H and therefore displacement reconstruction and current CZ area.
It applies with either `geometryReset='FOnly'` or `'FAndGeometry'`. The latter
still resets material glyph geometry to F_new*G0; it does not replace the
cached cohesive anchor. Integration-domain scaling and oversized-domain resets
no longer erase cohesive surface normals/positions. Active or damaged CZ
surface particles are still excluded from subdivision because cached mapping
and constitutive-history transfer are not implemented.

A particle contributing to one failed CZ node and one intact node now retains
its intact force and cohesive flag. Previously, the failed node could disable
the whole particle. The existing contact handoff is particle-based: a partially
failed interface remains cohesive within the support of particles that still
have intact CZ nodes. This patch does not introduce per-node/per-pair contact
state or an independently tracked crack front.

`gridCohesiveArea` now reports current area vectors (full-ring units), rather
than the reference unit normals previously written under that name.
`gridCohesiveForce` contains force, despite the internal temporary array's
traction name. Both are diagnostics at cached reference-grid locations.
Silo/VTK geometry remains a meridional section; no revolved mesh is added.
Initialization reports the sum of the two sides' reference areas.

## Validation and target acceptance

Run from the PFW directory:

```bash
python3 tests/run_axisymmetric_cz_cpp_tests.py
python3 tests/run_cz_event_cpp_tests.py
python3 tests/run_axisymmetric_cz_writer_tests.py --output /tmp/axis-cz-writer
```

The standalone C++ test compiles production area routines and extracted
production CZ state/law, reference-F, reset, TypeB-anchor and force-gather
code with strict warnings, AddressSanitizer and UBSan. It checks exact axial
and cylindrical measures, 60 inclined-line numerical comparisons, axis
limits, 120 noncommuting reset/delayed-activation cases, quad/triangle anchors,
324 mixed-mode law cases, equal/opposite forces and partial/complete failure.
The event test compiles the production event class and extracts the scheduler
and CZ handler. It checks growing/shrinking timesteps, dependency delays,
restart-wrapper state persistence, finite end times and a separate event
requesting a new reference. Its repository/serialization adapters do not
constitute a full GEOS restart test.
The writer test runs elastic/debond variants with FLIP/FMPM at 1 and 3 emulated
ranks, checking generated XML, groups, CZ tags, interface vectors and ring
volumes, with identical particle bytes across rank counts.
The writer checks require `useEvents=1` in the emitted XML and reject omitted
or disabled activation switches before generation. Input surface flags mark
candidates; inspect runtime cohesive flags, nonzero initialized reference
area and subsequent CZ forces to confirm actual activation. In the example,
`particleStressYY` is axial tension and `particleStressZZ` is hoop stress;
the latter can remain zero for its `nu=0` material.

These are CPU adapter tests and writer collective emulation, not a full GEOS
build, GPU execution or a real distributed solver test. Run the included
`examples/axisymmetricCZ` on the rebuilt target executable, then compare:

1. FLIP and FMPM elastic and debond variants; no theta velocity, equal/opposite
   axial reactions, intact oscillation versus cohesive release.
2. Mesh/time-step refinement, particularly interface-edge area normalization.
   CZ stiffness can constrain stability beyond the bulk-material CFL limit;
   reduce the timestep when increasing interface stiffness.
3. One and multiple real MPI ranks, particle migration, and restart versus an
   uninterrupted run. Existing reference mappings/areas/F_CZ/history fields
   are reused; this patch does not claim a locally executed restart test.
4. Delayed activation and F resets in a GEOS run, checking no instantaneous
   jump in current CZ area, force or cached displacement.

The manual LaTeX and affected generated event description are updated. A new
manual PDF is not included.
