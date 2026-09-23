# 3D Weak Trace Rotated Bar

This folder is a rotated full-3D companion to `weakTracePlanarInterface3D`.
The material sample is a square prism running along the long diagonal of the
Cartesian box, with an internal material interface normal to that diagonal.
The grid box is deliberately larger than the prism cross-section so the central
mortar interface and its free side edges are several cells away from the
domain-boundary nodes.

The default run uses level-4 simplex surface refinement, native particle
domains, and `2x2x2` ranks.  The `traceContactGroups` variant assigns the two
materials to separate contact groups and enables the prescribed-surface C0
vertex-mortar weak interface.  The `singleField` variant keeps the same refined
particle geometry but disables contact and weak-interface projection.

```bash
cd /g/g19/homel1/geos/scripts/preProcessing/materialPointMethodParticleFileWriter/verification/weakTraceRotatedBar3D
./runTest --force --jobs 1 --case-id weakTraceRotatedBar3D_simplexL4_native_diagPad
```

Use `--variant traceContactGroups` to run only the mortar case, or omit
`--variant` to run it with the single-field baseline.  Useful overrides are
`WEAK_TRACE_BAR_HALF_WIDTH`, `WEAK_TRACE_BOX_HALF_EXTENT`,
`WEAK_TRACE_SIDE_BUFFER_CELLS`, `WEAK_TRACE_SURFACE_REFINEMENT_LEVEL`, and
`WEAK_TRACE_FINAL_STRAIN`.
