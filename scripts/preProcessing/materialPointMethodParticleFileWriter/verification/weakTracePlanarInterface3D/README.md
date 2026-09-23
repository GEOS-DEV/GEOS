# 3D Weak Trace Planar Interface

This folder is a smoke test for the full-3D prescribed-surface weak trace/C0
vertex-mortar path.  It uses two elastic half-space blocks inside one Cartesian
box, split by an analytic plane normal to the x-loading direction.  The x extent
includes material buffer cells between the split plane and the moving boundary
faces so weak-tied interface nodes do not also see prescribed boundary
conditions.  The y/z extents include empty side-buffer cells outside the
material sample so interface-mapped nodes at the free side edges are not also
domain-boundary nodes.  The `traceContactGroups` variant gives the two sides
separate contact groups and writes particle-carried
`WeakDiscontinuity` surface normals and positions.

The default run uses the level-4 simplex-refined interface on `2x2x2` ranks,
which keeps the refined particle count near 30-50k per rank:

```bash
cd /g/g19/homel1/geos/scripts/preProcessing/materialPointMethodParticleFileWriter/verification/weakTracePlanarInterface3D
./runTest --force --jobs 1 --case-id weakTracePlanarInterface3D_simplexL4_native
```

Use `--variant traceContactGroups` to run only the C0 mortar/contact-groups
case, or omit `--variant` to run it with the single-field baseline.

This benchmark defaults to the two-pass PFW surface refiner with
`surfaceRefinementMethod="simplex"` and `surfaceRefinementLevel=4`, plus
`particleReferenceGeometry="fromRefinement"` and
`particleDomainOutput="native"` so VisIt can show the native refined particle
domains.  Override those with `WEAK_TRACE_SURFACE_REFINEMENT_METHOD`,
`WEAK_TRACE_SURFACE_REFINEMENT_LEVEL`,
`WEAK_TRACE_PARTICLE_REFERENCE_GEOMETRY`, and
`WEAK_TRACE_PARTICLE_DOMAIN_OUTPUT`.  The buffer defaults to four x cells per
side and can be changed with `WEAK_TRACE_BOUNDARY_BUFFER_CELLS`.  The empty
side buffer defaults to four y/z cells per side and can be changed with
`WEAK_TRACE_SIDE_BUFFER_CELLS`.
`WEAK_TRACE_INTERFACE_OFFSET` controls the plane location; the default keeps it
slightly off the central grid plane.

When simplex refinement is enabled, both variants expose the analytic interface
to PFW through surface flags so the single-field baseline can use the same
refined particle geometry.  The single-field variant still leaves contact and
weak-interface projection disabled.
