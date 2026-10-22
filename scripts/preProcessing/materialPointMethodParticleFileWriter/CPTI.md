# CPTI support

This incremental update targets the successfully built unified-geometry and
geometry-aware splitting source. `particleType=3` / XML `CPTI` now selects
affine triangle interpolation in native 2D and tetrahedral interpolation in 3D.
CPDI2 remains unsupported. Plane strain now uses a native XY mesh, not a
thin-3D specialization; see TRUE_2D_PLAN.md.

## Stored data and interpolation

No new particle fields, file columns, or restart schema are introduced. CPTI
requires shape 4 (2D) or 2 (3D), with full current/reference edge generators G/G0 stored by
rows. The 3D center plus three edges reconstructs four vertices:

    x0 = center - (g1 + g2 + g3)/4
    xa = x0 + ga, a=1,2,3

For E=[g1 g2 g3] in column notation, barycentric gradients are
`ba = E^(-T) ea`, `b0 = -(b1+b2+b3)`. The grid weights and gradients are:

    S_I      = (N_I(x0) + N_I(x1) + N_I(x2) + N_I(x3))/4
    grad S_I = sum_a N_I(xa) ba

Here N is the trilinear grid basis evaluated at each tetrahedron vertex.
Both expressions use the same integration geometry. The determinant is geometric;
the gradient is not divided by independently stored physical particle volume.
The 32 raw entries (four vertices times eight grid nodes) are coalesced by the
existing effective-map code. Cached mapping, on-demand mapping, and nodal
neighborhoods use the same CPTI kernel. Outflow checks use the same projected
vertices. Invalid/inverted/nonfinite geometry or invalid stencils return a fully
cleared mapping row before accessing the grid map; the existing bad-map policy
then handles active particles.

In 2D, two edges reconstruct v0=c-(g1+g2)/3 and two other vertices.
Weights are 1/3, barycentric gradients use the 2x2 reciprocal edge basis, and
three bilinear corner stencils give 12 raw entries. Physical volume is area*t;
t does not enter the interpolation weights or gradient coefficients.

The 3D formulation is CPDI-Tet4; see the manual's CPTI section and
[Noh et al., section 2.3.3](https://doi.org/10.1007/s00466-026-02804-8).
The present implementation convects independent particle domains using the
existing physical increment. It does not introduce a shared tetrahedral mesh
or independently advected vertex state.

## Integration scaling, resets and splitting

`cpdiDomainScaling=0` uses current material G unchanged. When enabled, CPTI
uniformly scales its integration domain about its centroid to an actual
maximum-edge diameter limit:

| `cpdiDomainScalingType` | CPTI diameter limit | CPTI algorithm |
|---|---|---|
| `homel` | `0.99998*min(active grid spacings)` | Uniform scale `s=min(1,limit/diameter)` |
| `vpHencky` | `0.99999*neighborRadius` | Same uniform scale |

These selections retain their existing box algorithms for CPDI. **CPTI does
not use the box VP-Hencky volume-preserving projection or box feasibility cap.**
Its integration measure changes by s^d in dimension d. Material G/G0, F, mass and physical volume
are untouched. Glyphs and splitting continue to use unscaled material geometry.
The existing scaled flag/reset selection and explicit surface-data disabling
control include CPTI. Scaling does not repair tetrahedral aspect ratio.

`geometryReset='FOnly'` preserves G/G0. `'FAndGeometry'` assigns G=G0*Fnew^T.
As before, a simplex requesting VP-Hencky **F selection** uses the isotropic-polar
reset fallback. Mapping always reads the resulting current G.

Explicit `subdivideParticles=1` and `oversizedParticleTreatment='split'` use
native triangle/tetrahedron longest-edge subdivision and preserve particle type 3.
Automatic `subdivideParticles=-1` remains specific to eligible CPDI boxes and
does not select CPTI simplices. The oversized policy is independent of that automatic
setting. Existing per-pass limits, rank-count fallback, conservation rules and
cohesive exclusions remain as documented in `PARTICLE_SPLITTING.md`.

## PFW: all tetrahedra or mixed bulk/surface interpolation

All-tetrahedron objects:

```python
pfw['planeStrain'] = 0
pfw['particleGeometry'] = 'fromRefinement'
# On each desired object:
# particleType=3
```

In 3D each Cartesian piece becomes six equal-volume tetrahedra around a consistent
lower-to-upper box diagonal. Fields are sampled at the final centroids, and
material/contact labels are checked against the parent piece. Surface offsets
are evaluated relative to each final centroid. Integrated traction is apportioned
by one sixth. In native 2D, two equal-area triangles replace each Cartesian
piece and each gets half its integrated traction; select formulation=PlaneStrain
and one planeStrainThickness. Native cut triangles/tetrahedra are emitted directly. The ppc setting controls
the seed lattice; six CPTI particles replace each otherwise Cartesian seed piece.

For CPDI boxes with CPTI surface tets, the named preset is:

```python
pfw['planeStrain'] = 0
pfw['surfaceRefinementMethod'] = 'hexTet'
pfw['surfaceRefinementLevel'] = 2
# Objects use particleType=2 for CPDI hex bulk.
```

See [HEX_BULK_TET_SURFACE_PLAN.md](HEX_BULK_TET_SURFACE_PLAN.md) for validation,
layout diagnostics and the complete preset example. The equivalent explicit
configuration remains valid:

```python
pfw['planeStrain'] = 0
pfw['particleGeometry'] = 'fromRefinement'
pfw['surfaceRefinementMethod'] = 'simplex'
pfw['surfaceRefinementLevel'] = 2
pfw['simplexParticleType'] = 3
# Objects use particleType=2 for CPDI bulk.
```

`simplexParticleType` is PFW-only and normally defaults to None (inherit the
object's type); the `hexTet` preset resolves an unset value to 3. Values 0, 1 and 3 are accepted. The override affects actual cut simplices;
uncut and Cartesian-fallback pieces retain their object type. CPTI requires
native geometry of the matching dimension. PFW declares every required material/type pair in the XML,
including on ranks without those particles. No solver XML option is added.

Existing Python inputs retain their defaults and aliases. Geometry schema 2
remains unchanged. New CPTI runs require regenerating particle/XML files with
this writer. Existing schema-2 point/CPDI cases do not require a format conversion.

Examples:

- `examples/cpti/pfw_input_cpti.py`: all-CPTI sphere.
- `examples/cpti/pfw_input_cpti_mixed.py`: CPDI bulk plus CPTI cut tets.

The native mixed example is examples/true2D/pfw_input_true_2d.py.
These examples disable job submission. Use the normal userDefs configuration and PFW
invocation. The refinement JSON distinguishes `emitted_refinement_pieces` from
`emitted_particles`, because CPTI tessellation may turn one piece into six rows.

## Verification and remaining integration checks

From this directory:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/run_cpti_cpp_tests.py
python3 tests/run_cpdi_geometry_cpp_tests.py
python3 tests/run_subdivision_cpp_tests.py
python3 tests/run_cpti_writer_tests.py --output /tmp/geos-cpti-writer-check
# Optional address/undefined sanitizers; no leak-check claim:
ASAN_OPTIONS=detect_leaks=0 python3 tests/run_cpti_cpp_tests.py --sanitize
```

The CPTI C++ runner compiles actual production function bodies and geometry
headers with small host adapters. It checks 4,300 mapping cases against independent
Gaussian-elimination/trilinear references: partition of unity, affine-field
reproduction, gradient moments, duplicate-node accumulation, raw/view output
agreement, rotations, scaling, split children, both resets, invalid inputs and
boundary roundoff. Python tests check six-tet volume/centroid conservation,
shared faces and independent sampled ownership. The splitter runner includes
CPTI explicit/automatic/oversized cases and transfer conservation. Old CPDI
mapping cases are retained as regression checks.

The full writer runner verifies all-tet Cartesian/simplex and mixed CPDI/CPTI
generation on 1/3/7 **emulated** ranks, including zero-row ranks, with identical
particle bytes, positive volumes, valid shape/type pairs and XML declarations.
These are not actual MPI executions.

The full GEOS dependency tree/build is unavailable in this workspace. The patch
still needs the target GEOS build and mechanical smoke tests: uniform translation,
affine deformation, a mixed CPDI/CPTI run, splitting during motion, a checkpoint
round trip and multi-rank migration. GPU and coupled contact/CZ integration are
also unverified. No new cohesive subdivision, facet ownership or shared-vertex
tet-mesh continuity is claimed.
