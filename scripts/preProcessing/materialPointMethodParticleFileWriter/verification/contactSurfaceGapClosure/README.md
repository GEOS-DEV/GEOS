# Contact surface/gap-closure verification

This folder owns a compact folder-scoped verification test for curved MPM contact surfaces.
The default `runTest` dispatches three representative variants; use `--all-variants` to run
the full 12-case matrix.

```bash
./runTest --case-id contactSurfaceGapClosure --force --no-visit
./runTest --case-id contactSurfaceGapClosure --force --all-variants
```

The PFW input exposes only a small set of test-defining environment variables: `CONTACT_GAP_REFINE`,
`CONTACT_GAP_CPP`, `CONTACT_GAP_PPC`, the contact field/surface/gap-correction variant values
set by `runTest`, and the geometric gap/compression parameters documented in the input.
