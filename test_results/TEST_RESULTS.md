# Test results

## Unit and geometry tests

`python -m unittest discover -s tests -p 'test_surface_refinement.py' -v`

All **18 tests passed**. Coverage includes:

- 200 random plane cuts in 2D/3D with positive volumes, parent-volume closure,
  and first-moment closure; additional cuts through vertices, edges, faces, and
  thin slices; translation, anisotropic dimensions, and plane-strain thickness.
- Empty-center recovery in both methods, including an object whose bounding
  interval excludes the coarse center but intersects the parent domain.
- Unflagged A/B adjacency detection, staged bulk-record reuse, and no extra
  geometry sampling in an unflagged uniform bulk region.
- Planar complementary two-material partitions and unequal material refinement
  factors, curved free surfaces/inclusions, strict fallback, and explicit
  non-detection of an isolated feature missed by every coarse center.
- Integer ownership, expanded object-factory bounds, periodic marker images,
  independent rank slabs, and ranks with zero owned rows.

## Complete-writer multi-process regression

All **24 configurations** passed. Each row below was tested with 1, 3, and 7
writer processes. Particle files were byte-identical to the one-process result.

| Fixture | Method | Particles |
|---|---|---:|
| 2D planar A/B | Cartesian | 80 |
| 2D planar A/B | Simplex | 113 |
| 2D separated objects | Cartesian | 42 |
| 2D separated objects | Simplex | 94 |
| 3D planar A/B | Cartesian | 480 |
| 3D planar A/B | Simplex | 1,382 |
| 3D separated objects | Cartesian | 252 |
| 3D separated objects | Simplex | 1,296 |

Checks included finite particle data, strictly positive volume reconstructed
from R-vectors, consecutive global IDs, total volume reduced over all processes,
all material/type blocks in rank 0's XML, and absence of writer-only options in
solver XML. The separated-object fixtures use a slice-filtered object factory;
seven processes on five coarse x rows exercise empty writer ranks.

**Transport limitation:** these runs used a deliberately labeled file-backed
test communicator, not mpi4py/MPICH/OpenMPI. They validate the complete writer's
control flow and decomposition consistency with independent Python processes,
but cannot validate MPI installation, actual collectives, cluster filesystems,
or speedup. `tests/run_writer_regression.py --communicator mpi` is supplied for
that next check on the user's machine. The production writer still uses its
original mpi4py communicator, not the test helper.

The first long regression invocation reached the execution-time limit after
15 successful cases. The remaining 9 cases were completed in two targeted
invocations. All 24 output sets were subsequently rechecked together; the
combined machine-readable results are in `writer_regression_results.json`.

## Supplied workshop-input PFW smoke test

The original `weakTracePlanarInterface3D` input was run through the original and
updated writers with `CPP=6`, `PPC=1`, and the prescribed trace-contact variant.
This is a reduced-resolution **preprocessor** smoke test, not a recreation of
the workshop mechanics results. The one-process launcher supplied a test-only
MPI API shim because the runtime was unavailable.

| Configuration | Particles | Initialized volume |
|---|---:|---:|
| Original, surface refinement off | 64 | 0.42 |
| Updated, surface refinement off | 64 | 0.42 |
| Updated, Cartesian level 2 | 512 | 0.42 |
| Updated, simplex level 2 | 1,502 | 0.42 |

Volumes agree with the filled box to floating-point accumulation accuracy. The
two no-refinement particle files were byte-identical. The simplex case reported
90 successful cut subdomains and 1,080 simplex-derived particles, with no
fallback in this planar case. Remaining particles were ordinary uncut children.

## Packaging checks

- New/modified Python files compile.
- `git diff --check` passes.
- Forward `git apply --check` and `git apply` pass on the exact provided baseline.
- All nine installed file contents match the packaged replacements after apply.
- `pfw_geometryObjects.py` is unchanged.

No C++ solver build, GEOS time stepping, mechanical verification, or physical
performance claim is included in these results.
