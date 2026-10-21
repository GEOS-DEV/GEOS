# Unified particle geometry: schema 2

This patch replaces the previous optional native-glyph layer. Every particle
has one material shape and two minimal full-edge matrices. Constitutive F,
material geometry, and the temporary integration domain have separate roles.
Pre-schema-2 GEOS particle files, XML output switches, and restart files are not supported.
Regenerate particles and XML using the updated PFW. Old Python `pfw_input` files
remain accepted through normalization in PFW.

## Persistent data and reconstruction

| State | Storage in this patch | Purpose |
|---|---|---|
| `particleGeometryShape` | integer per particle | 0 box, 1 triangular wedge, 2 tetrahedron |
| `particleGeometryReference` (G0) | 3 x 3 doubles | Full reference edge generators |
| `particleGeometryCurrent` (G) | 3 x 3 doubles | Full current material edge generators |
| Center, mass, physical volume, F | Existing solver fields | Geometry does not duplicate these |
| Integration geometry, vertices, half edges | Stack temporaries | Reconstructed when needed |

Generators are rows in C++ and particle files. Thus `G = G0 F^T` when the
constitutive and material frames coincide. In column-vector notation the same
equation is `G = F G0`. That equality is deliberately not assumed after FOnly.

| Shape | Full generators | Vertices relative to center c | Initial measure |
|---|---|---|---|
| 3D box | g1, g2, g3 | c +/- g1/2 +/- g2/2 +/- g3/2 | det(G) |
| 3D tet | v1-v0, v2-v0, v3-v0 | v0=c-(g1+g2+g3)/4; vi=v0+gi | det(G)/6 |
| Thin-3D wedge | Two triangle edges, thickness vector g3 | triangle centroid offsets, then +/-g3/2 | det(G)/2 |
| Planned 2D parallelogram | Two 2-component edges | c +/- g1/2 +/- g2/2 | det(G) |
| Planned 2D triangle | Two 2-component edges | v0=c-(g1+g2)/3; vi=v0+gi | det(G)/2 |

The matrices cost 144 bytes per 3D particle, plus the shape tag. Boxes already
needed two 3x3 domain matrices, so their matrix count is unchanged. Native
simplices previously also carried a 3x3 glyph-reference matrix: removing it
saves 72 bytes per particle. Future native 2D storage will use two 2x2 matrices
(64 bytes total); this patch's running solver remains 3D/thin-3D.

## Physical update and reset

Advance G with the same physical incremental deformation used to advance F:
`G_next = G_current DeltaF^T`. Do not rebuild G from F and G0 every step.
Rigid rotation, particle-transform events, and overlap volume corrections also
advance current material geometry. Thermodynamic startup can establish G0 from
the prescribed initial F and the supplied current G; it is initialization, not
a numerical reset. Subdivision creates smaller G0/G for each child.

`geometryReset` is independent of `domainResetType`:

| Setting | F | Current G | Reference G0 |
|---|---|---|---|
| `FOnly` (default) | Set to selected F_new | Preserve | Preserve |
| `FAndGeometry` | Set to selected F_new | G0 F_new^T | Preserve |

The policy applies to automatic resets and the PolymerHeal/CrystalHeal reset
paths. Surface/material reference data and cohesive reference F retain their
existing frame updates; they are not the fixed geometry reference G0. A reset
does not change physical particle mass or volume. Geometry measure may differ
from physical volume after a deliberate geometry reset, especially if healing
has changed the constitutive reference volume. Subsequent convection preserves
that distinction. FOnly does not cure material-domain distortion; use the
geometry-reset policy or geometry-aware splitting when that is the intention.

The VP-Hencky F-selection algorithm applies to boxes. Simplex F selection uses
isotropic-polar reset; no box feasibility cap is imposed on a simplex.

## Interpolation and CPDI scaling

SinglePoint and SinglePointBSpline sample at the center while retaining the
same material geometry for output, splitting and CZ anchors. CPDI requires a
box. On every mapping call, construct a temporary full-edge integration matrix
from G, apply the selected projection, and use local half edges `r=Ghat/2` in
the existing eight-corner CPDI weights and gradient formula. Cached mapping,
on-demand mapping, neighbor support and CPDI boundary checks use this same
projection. Cached shape functions remain cached; Ghat is not a particle field.

Homel clipping and VP-Hencky projection retain their existing geometry-dependent
algorithms. They are generally anisotropic and cannot be represented by one
scalar scale factor. VP-Hencky uses the full-edge feasibility bound
`measure <= (D/sqrt(d))**d`; the half-edge mathematical convention divides this
measure by `2**d`. Integration scaling never changes G, G0, F, mass or physical
volume. A projection failure retains its input and reports an invalid-domain
diagnostic. The existing explicit surface-data disabling control remains a
separate side effect. Output always reconstructs unscaled material geometry.

CPTI is a planned interpolation kernel, not enabled by this patch. For an affine
tet its domain weights will be 1/4 and barycentric gradients will be
`b0=-G_column^-T [1,1,1]`, `bi=G_column^-T ei`. In 2D use 1/3 and two edges.
Initial simplex scaling should be uniform and use actual vertices/radial bounds;
it must not reuse the box VP-Hencky measure cap. Arbitrary warped CPDI2 cells,
native 2D mesh topology/mapping, and axisymmetric mechanics remain future work.

## Splitting and cohesive zones

Runtime splitting uses unscaled current G. Boxes retain the generator-mask
split policy. Tets bisect their longest current edge; wedges bisect a triangle
edge or, in 3D, their thickness. Plane strain preserves wedge thickness.
Native children are reevaluated and the largest violating child is split until
all meet the diameter limit or the pass reaches four leaves in plane strain /
eight in 3D. A single tet bisection need not reduce maximum diameter. The log
reports selected parents that still have oversized children; later steps can
continue refinement.

A temporary child map (A,q,w) gives row-stored G_child=A G and G0_child=A G0,
and current offset G^T q. The same map acts on both frames after either reset
policy. No new persistent geometry or file columns are required. Fraction w
partitions mass and physical current/reference volumes even if geometry measure
differs from physical volume. Adaptive child fractions can differ within a family.

Intensive state, specific internal energy, F and material history are copied
using particle-sized wrappers, including the `meltFlag` field. The optional
velocity-gradient correction conserves linear momentum through the weighted
centroid identity; kinetic energy is recomputed and can increase. With correction
off, parent velocity and point-particle kinetic energy are preserved. There is no
energy renormalization or general angular-momentum conservation claim.

`particleSurfaceTraction` contributes an integrated force in grid assembly:
current and reference values are partitioned by w. Surface offsets preserve the
parent's approximate absolute surface point; constitutive-reference offsets use
F^-1 times the current child offset, not the fixed G0 offset. This does not infer
exact child boundary-facet ownership. Cohesive and damaged-cohesive particles,
and particles with active cohesive state, remain excluded. Mapping and neighbor
caches are invalidated; CPDI diagnostics are refreshed before subsequent reset
selection. See `PARTICLE_SPLITTING.md` for controls and verification.

CZ TypeA caches the supplied surface offset. TypeB now intersects the ray from
the particle center along its current surface normal with the actual unscaled
box, wedge or tet at CZ initialization. Both convect that cached anchor using
relative cohesive F; Nodal is unchanged. This changes TypeB's anchor selection
for a deformed particle relative to the old reference-box ray. The anchor is a
particle-boundary estimate, not explicit ownership of an interface facet.
Refinement does not turn internal simplex faces into cohesive interfaces.

## PFW and file contract

Preferred input:

```python
pfw['particleGeometry'] = 'fromRefinement'  # default: box
pfw['surfaceRefinementMethod'] = 'simplex'
pfw['surfaceRefinementLevel'] = 2
pfw['geometryReset'] = 'FOnly'              # or FAndGeometry
```

PFW simplex cuts still require object `particleType=1` (SinglePointBSpline).
Uncut/fallback pieces remain boxes, producing hex bulk plus surface tets in 3D
and boxes plus wedges in the current plane-strain mode. GEOS also accepts point
linear interpolation with supplied canonical simplex geometry; PFW's cut
algorithm retains its B-spline restriction. CPDI with a simplex is rejected.

Old `particleReferenceGeometry='equivalentBox'|'fromRefinement'` maps to the new
selection unless `particleGeometry` is explicitly supplied. Old
`particleGeometryFormat='legacy'|'unified'` names are accepted but both emit schema
2. `RVector` in `particleFileFields` is consumed as an obsolete geometry request.
Old `particleDomainOutput` and `plotUnscaledParticles` inputs are accepted by PFW;
they are not emitted to XML. Output always shows actual material geometry.

Required geometry columns, on every row and for every interpolation type:

```text
GeometryVersion GeometryShape ParticleVolume
GeometryXX GeometryXY GeometryXZ
GeometryYX GeometryYY GeometryYZ
GeometryZX GeometryZY GeometryZZ
```

Version is exactly 2. The reader also requires ID and PositionX/Y/Z, rejects
unknown/partial old geometry schemas, nonfinite or inverted generators, and
initial measure/volume mismatch beyond relative tolerance 1e-9. Initial wedges
must be xy triangles extruded in z. Both geometry matrices are restartable;
restart version 2 is required. There is no old-restart conversion path.

The XML addition is `geometryReset="FOnly"` on `SolidMechanics_MPM`. The old
ParticleRegion output policy is removed. See the complete runnable input in
`examples/unifiedGeometry/pfw_input_unified_geometry.py`.

## Verification

From this directory:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python tests/run_reference_geometry_cpp_tests.py
python tests/run_cpdi_geometry_cpp_tests.py
python tests/run_writer_regression.py --communicator emulated --reference-geometry fromRefinement --output /tmp/pfw-native
python tests/run_reference_geometry_cpp_tests.py --writer-results /tmp/pfw-native
```

The standalone tests compile production geometry and mapping bodies with minimal
CPU adapters. Writer collectives are explicitly emulated unless real MPI is
selected. These checks do not replace a complete GEOS build, MPI migration,
restart round trip, GPU execution or a coupled CZ mechanics regression.
