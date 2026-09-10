# 3D Weak Trace Planar Interface

This folder is a smoke test for the full-3D prescribed-surface weak trace/C0
vertex-mortar path.  It uses two elastic half-space blocks inside one Cartesian
box, split by an oblique analytic plane.  The `traceContactGroups` variant gives
the two sides separate contact groups and writes particle-carried
`WeakDiscontinuity` surface normals and positions.

The first intended run is the trace-only pdebug smoke test:

```bash
cd /g/g19/homel1/geos/scripts/preProcessing/materialPointMethodParticleFileWriter/verification/weakTracePlanarInterface3D
./runTest --force --jobs 1 --variant traceContactGroups --no-visit
```

For a baseline comparison, add `--variant singleField` or omit `--variant`.
