# Weak Trace Large-Deformation Compaction

This benchmark is a qualitative large-deformation stress test for the
duplicated-field C0 vertex-mortar weak-interface path.  A copper matrix with two
aluminum inclusions is compressed between tool-steel platens in plane strain.
The lateral faces are free so the compact can flow into the side void as the
platen separation is reduced.

The `traceContactGroups` variant uses contact group 1 for the copper matrix and
contact group 0 for the tool-steel platens and aluminum inclusions, with the C0
vertex-mortar tie applied to `weakInterfaceTracePairs = [[0, 1]]`.  The
`singleField` variant keeps the same weak surface markers and surface-refined
particle sampling, but collapses all particles to contact group 0 and disables
contact and trace projection.

Default controls:

- `WEAK_TRACE_VARIANT=traceContactGroups`
- `WEAK_TRACE_PARTICLE_TYPE=1` for `SinglePointBSpline`
- `WEAK_TRACE_PPC=4`
- `WEAK_TRACE_SURFACE_REFINEMENT_LEVEL=2`
- `WEAK_TRACE_FINAL_Y_STRETCH=0.62`
- one MPI rank, `planeStrain=1`, `nK=3`

Example:

```bash
cd ~/geos/scripts/preProcessing/materialPointMethodParticleFileWriter/verification/weakTraceLargeDeformationCompaction
WEAK_TRACE_PPC=4 WEAK_TRACE_SURFACE_REFINEMENT_LEVEL=2 ./runTest --force --jobs 1 --walltime 00:30:00 --post-walltime 00:10:00
```

Single-field comparison only:

```bash
cd ~/geos/scripts/preProcessing/materialPointMethodParticleFileWriter/verification/weakTraceLargeDeformationCompaction
WEAK_TRACE_PPC=4 WEAK_TRACE_SURFACE_REFINEMENT_LEVEL=2 ./runTest --variant singleField --force --jobs 1 --case-id weakTraceLargeDeformationCompaction_singleField --walltime 00:30:00 --post-walltime 00:10:00
```
