# Hex bulk with tetrahedral surface refinement

This incremental patch targets the CPTI-patched branch. The CPTI update already
supplies the mixed CPDI-box/CPTI-tet generation and runtime support. This patch
makes that combination a named, checked PFW workflow and adds emitted-layout
diagnostics, a complete example, regression tests and current documentation.
It does not change GEOS C++ source, the cut algorithm or the particle schema.

## Implementation plan and resulting flow

1. **Resolve the requested layout before generation.** The new
   `surfaceRefinementMethod='hexTet'` selects the existing simplex cutter,
   `particleGeometry='fromRefinement'`, and cut-only `simplexParticleType=3`.
   Validate 3D and reject contradictory explicit geometry or all-tet bulk types.
2. **Keep the existing bounded two-pass detection.** Sample coarse centers for
   occupied/empty and region-label changes and surface flags. Use the three-slice
   cache and one-row halos. Reuse staged bulk records outside the selected band.
3. **Refine candidate parents and cut only intersected subdomains.** Apply the
   surface factor and common neighboring material factor. Keep uncut pieces as
   boxes; partition each accepted local plane cut into native tetrahedra.
   Material interfaces retain both sides, free surfaces omit the void side.
4. **Preserve field and geometry contracts.** Existing final-centroid sampling,
   labels, volume integration and surface offsets remain in use. Use schema 2
   full-edge generators. Declare CPDI/CPTI material/type pairs in generated XML.
5. **Report what was actually emitted.** Count final serialized records by
   shape and interpolation, with physical-volume sums and minimum/maximum
   particle volume. Gather only these small rank summaries, including empty
   ranks; retain existing candidate, cut and fallback counters.
6. **Verify the workflow and update documentation.** Exercise free surfaces,
   interfaces, curvature and fallback, rank-independent output and old inputs.
   Include the preset example, generated option tables and rebuilt manual.

The patch implements all six items. Solver-side mixture mapping, physical
geometry evolution, output reconstruction and native splitting come from the
unified-geometry, splitting and CPTI baseline.

## Recommended input

In a complete PFW input:

```python
pfw['planeStrain'] = 0
pfw['surfaceRefinementMethod'] = 'hexTet'
pfw['surfaceRefinementLevel'] = 2
# Set particleType=2 on the geometry objects for CPDI bulk.
```

No separate geometry/type override is needed. An existing explicit
`particleGeometry='fromRefinement'` remains valid. Remove an explicit
`particleGeometry='box'` or legacy `particleReferenceGeometry='equivalentBox'`,
or replace it with native geometry. As before, the preferred `particleGeometry`
key takes precedence over the legacy alias when both are supplied.

The equivalent explicit configuration remains supported and produces identical
particle records:

```python
pfw['planeStrain'] = 0
pfw['surfaceRefinementMethod'] = 'simplex'
pfw['particleGeometry'] = 'fromRefinement'
pfw['simplexParticleType'] = 3
pfw['surfaceRefinementLevel'] = 2
# Geometry objects: particleType=2.
```

The preset accepts object types 0/1 for point-interpolated hex bulk and explicit
`simplexParticleType=0` or `1` for point-interpolated tets in comparisons. Its
normal combination is type 2 / type 3. Type 3 objects are rejected in this mode
because they tetrahedralize all Cartesian pieces. Use the existing ordinary
`simplex`/`cartesian` workflows for all-CPTI generation instead. No plane-strain
CPTI specialization is introduced; true 2D remains a separate task.

`surfaceRefinementLevel` retains its default of 1 and is a factor per active
direction, not recursive depth. `hexTet` cuts even at level 1. A level greater
than 1 also creates smaller hexes in uncut portions of the candidate band.
Material-indexed `particleRefinement` continues to apply in the bulk and uses
the maximum neighboring material factor for a common interface partition.

## Geometry, interpolation and XML

| Piece | Geometry / shape ID | Default interpolation | Stored generators |
|---|---|---|---|
| Interior outside surface band | Hex / 0 | CPDI / 2 | Three full box edges |
| Uncut subdomain in the band | Smaller hex / 0 | CPDI / 2 | Three full box edges |
| Accepted cut material piece | Tet / 2 | CPTI / 3 | Three edges from vertex 0 |
| Whole-subdomain Cartesian fallback | Hex / 0 | CPDI / 2 | Three full box edges |
| Void | No particle | None | None |

The particle centroid, shape, physical volume and nine generator components
reconstruct the geometry. Runtime G/G0 storage, scaling, F-only/F-and-geometry
resets and splitting are unchanged. Current material geometry supplies output
and splitting; temporary integration scaling never overwrites it.

Generated XML uses the already implemented CPDI and CPTI particle types and
matching subregions for each declared material. PFW may declare a CPTI subregion
that ends up empty, for example when every attempted cut falls back. The preset
and layout diagnostics add no solver XML fields and no particle-file columns.

## Diagnostics and fallback

The existing `mpmParticleFile_<case>_refinement.json` gains:

- `requested_method`, `resolved_geometry`, `resolved_simplex_particle_type`:
  distinguish the user's preset from resolved `method='simplex'`.
- `particle_layout.total_count` and `.total_volume`: actual serialized rows
  and their physical volume, checked against the existing emitted count.
- `particle_layout.by_type_and_shape`: one entry per actual pair, containing
  `particle_type`, `geometry_shape`, `shape_name`, `count`, `volume`,
  `min_particle_volume` and `max_particle_volume`.

These diagnostics apply to all two-pass refinement modes, not just hexTet.
They do not parse the completed file or gather per-particle arrays. Per-rank
compensated sums reduce accumulation error; aggregate volume can differ by
summation roundoff with rank count. Particle bytes and integer counts remain
rank-independent for deterministic geometry and property functions.

The original `fallback_*` counters count failed subdomains, not emitted hexes.
A failed subdomain is classified at its midpoint, so it may emit a hex or no
particle. The default fallback is `cartesian`; choose
`surfaceRefinementFallback='error'` to stop at the first unresolved cut.
No individual sliver is discarded by an arbitrary volume threshold, and no
mass/volume renormalization conceals clipping failures. Tiny tets can have
consequences for downstream mechanics; this patch reports their minimum volume
but does not change the accepted partition.

## Verification and target checks

Run from the writer directory:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/run_hex_tet_writer_tests.py --output /tmp/geos-hex-tet-check
python3 tests/run_cpti_writer_tests.py --output /tmp/geos-cpti-check
```

The new full-writer runner covers unchanged bulk, planar free surfaces, complementary material
interfaces, a curved free surface and missing-plane fallback at 1/3/7 emulated
ranks. It checks exact planar volume, positive geometry measures, recovered tet
vertices, correct type/shape pairs, summary/file agreement, valid XML and
byte-identical particles across ranks. It compares all five scenes with the
explicit configuration and tests rejection of incompatible inputs and strict
fallback. These processes emulate collectives; they do not establish real MPI
transport, solver mechanics or performance.

The complete example is
`examples/hexBulkTetSurface/pfw_input_hex_bulk_tet_surface.py`.
Its job submission is disabled. Run it with the usual PFW userDefs configuration;
regenerate particle and XML inputs to select the new layout.

After the CPTI build completes, use that executable for a mixed-domain mechanical
smoke run, including translation/affine deformation, contact, runtime splitting
if enabled, restart and migration. This incremental patch changes only PFW and
its documentation/tests, so it needs no additional C++ rebuild. It does not
replace those target integration checks.

## Scope for later work

This is particle initialization, not a conforming shared-vertex tet mesh. Planes
are independently fitted per subdomain; global surface continuity and mesh
quality are not guaranteed. Features missed by coarse occupancy/flags remain
resolution dependent. Multi-region junctions or inconsistent local planes use
the documented fallback. No general runtime conversion of a newly exposed hex
surface into tets, cohesive facet ownership/subdivision, true 2D or axisymmetry
is introduced. Those tasks should build on this mixed-layout contract.

Native 2D now has the analogous quad/triangle workflow via explicit simplex,
fromRefinement and simplexParticleType=3 settings. The hexTet preset remains
3D only. See TRUE_2D_PLAN.md and PFW_2D_OBJECTS.md.
