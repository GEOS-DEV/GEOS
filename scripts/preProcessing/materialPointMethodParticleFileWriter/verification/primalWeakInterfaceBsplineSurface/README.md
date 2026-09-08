# Primal weak-interface B-spline surface-reconstruction checks

This convenience runner dispatches the three serial CPU verification variants
used by the particle-surface reconstructed primal weak-interface prototype:

- `weakTraceRotatedBar/primalEnrichment`
- `weakTraceDynamicWave/primalEnrichment`
- `weakTraceCircularInclusionRadial/primalEnrichment`

All command-line options are forwarded to the existing case runners. For
example:

```bash
./runTest -y --walltime 00:20:00 --post-walltime 00:10:00
```

The solver reconstructs a current, continuous two-dimensional polyline from
`SurfaceFlag::WeakDiscontinuity` particle surface points and normals at the
beginning of every explicit step. Adjacent interface cells share canonical
edge-crossing points. The mass-orthogonal primal ridge uses the unsigned
distance and one-sided gradient to that reconstructed polyline.

This is an experimental serial, CPU-only, plane-strain, PIC,
`SinglePointBSpline` implementation. Contact, DFG, CZ, and the dual
weak-interface projection must remain disabled in these variants.
