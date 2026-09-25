# Boundary-mortar performance update (2026-09-25)

## Scope and baseline

Incremental update against the exact files in `2609251035_geos_mpm_bundle.zip`,
review commit `5f2aee3aee772a22dcfb64c7f79b56a046aaf03a` with the supplied worktree
overlay. Apply to that current tree; do not reapply the earlier PFW, native
geometry, boundary-mortar, or fallback fixes.

This step accelerates the existing CPU boundary-compatible projection. It does
not change particle quadrature, P2G/G2P, stress updates, surface fitting, mortar
integration, fallback selection, boundary dependencies, or the accepted physical
residual tests. Each Cartesian component remains one complete projection over
all accepted interfaces. Per-interface/block scheduling is deliberately deferred.

Standalone kernel and explicitly emulated scheduling tests were performed.
A complete GEOS build, MPI transport, timestep/mechanics comparisons, GPU runs,
and restart/migration tests were NOT performed in this environment. The helper
and limited CMake build tests are not substitutes for that integration validation.

## Runtime controls

For an existing boundary-compatible input, use:

```python
pfw["weakInterfaceTraceBoundaryCoupling"] = 1  # already used by the test
pfw["weakInterfaceTraceQRBackend"] = "auto"
pfw["weakInterfaceTraceComponentRanks"] = 3
pfw["weakInterfaceTraceFastMarkerSearch"] = 1
pfw["weakInterfaceTracePerformanceLog"] = 1    # short profiling runs; default 0
```

These four new keys are forwarded to the solver XML. Their defaults in the new
solver are `auto`, `3`, `1`, and `0`. They affect only the boundary-coupled branch.
Existing particle files can be reused: add/change these attributes on the
SolidMechanicsMPM solver in the XML, or regenerate the XML with the updated PFW.
No new particle fields, changed volumes, or new geometry queries are required.

For a reference comparison using the new executable:

```python
pfw["weakInterfaceTraceQRBackend"] = "reference"
pfw["weakInterfaceTraceComponentRanks"] = 1
pfw["weakInterfaceTraceFastMarkerSearch"] = 0
pfw["weakInterfaceTracePerformanceLog"] = 1
```

Keep rank tolerance, jump tolerances, particle file, grid, boundary gap, loading,
and transfer method unchanged during performance comparisons.

### QR backends

- `reference`: unchanged original hand-written pivoted Householder kernel.
- `native`: dependency-free faster kernel using safeguarded partial-column-norm
  updates. Norms are refreshed when cancellation is significant and before a
  rank decision can discard remaining directions. Primal/dual recovery and all
  original-row acceptance checks remain in place.
- `lapack`: LP64 `DGEQP3`, triangular solves, and `DORMQR`. No explicit normal
  equations or explicit Q. Requires the optional build flag below.
- `auto`: LAPACK if compiled, otherwise native. This is a build-based choice,
  not a runtime autotuner. Benchmark `native` and `lapack` on the actual host;
  a reference/system LAPACK need not be faster than the native kernel.

There is no automatic fallback to GS, no relaxed constraint tolerance, and no
added compliance. Near a numerical-rank threshold, different pivot schemes can
choose different directions; inspect rank and the ORIGINAL accepted-row
residuals. Multipliers can differ for redundant rows; compare physical velocity,
impulses, work, and reactions rather than demanding bitwise dual equality.

The numerical factors are rebuilt for the current operator every step. The
patch does not freeze geometry or reuse stale mass-dependent factors.

## Build

The native/reference backends build with the existing minimal-MPM profile and
need no added numerical dependency. Reconfigure the existing build so CMake
includes the new QR translation unit, then rebuild the CPU executable.

For LAPACK, append to the normal GEOS CMake configuration:

```text
-DGEOS_ENABLE_MPM_MORTAR_LAPACK=ON
```

The new option requests BLAS/LAPACK through the existing GEOS TPL setup even when
`GEOS_ENABLE_DENSE_LINEARALGEBRA=OFF`. It does NOT turn on the implicit linear
algebra stack or change the minimal-MPM profile's other choices. Provide the
usual GEOS host-config BLAS/LAPACK library settings. The backend requires LP64
(32-bit Fortran integers); ILP64 libraries with identical symbol names are not
supported. CMake rejects `BLA_SIZEOF_INTEGER=8` when specified, but the caller
must also select compatible libraries. Underscored and the existing GEOS
no-underscore Fortran naming conventions are handled.

Selecting `lapack` in XML without that build support is an explicit error.
Start with one BLAS thread per MPI process to avoid oversubscription:

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
```

The code does not change global BLAS thread settings at runtime.

## Concurrent MPI component scheduling

The existing gathered canonical-key state is retained on all ranks. Once it is
ready, each worker computes its assigned complete component problems BEFORE
any result broadcast. This is not a root change inside a serial solve/broadcast
loop.

| MPI ranks | requested component ranks | Assignment in 3D |
|---|---:|---|
| 1 | 3 | Rank 0 computes x, y, z |
| 2 | 3 | Rank 0 computes x,z; rank 1 computes y |
| 3 or more | 3 | Ranks 0,1,2 compute x,y,z concurrently |
| Any | 1 | Rank 0 computes x,y,z (reference scheduling) |

All ranks still perform their usual spatially partitioned GEOS work. On four
or eight ranks, only this component stage uses at most three workers. The
number of component workers is independent of the Cartesian partition shape:
there is no need to decompose the domain into a long 3x1x1 layout.

Preserve the physical grid and choose a valid roughly isotropic decomposition.
For about 161k total particles, four ranks average about 40k particles/rank, but
surface refinement and the octant geometry can make local loads uneven. Check
max/min owned particles and interface work, not just the average. Match the
existing -x/-y/-z partition controls and the generator's decomposition settings.

All worker success/failure statuses are exchanged before ANY projection results
are committed. Failure on one worker produces the same error on every rank.
Nonworker ranks and ranks without local interface still participate. Receive
allocation failure is checked collectively. Each solve occurs exactly once;
results, scalar diagnostics, and reaction data are then broadcast from its owner.

The dense matrix cap is per component. Additional storage is three sets of
Result arrays per rank; dense factors remain on their workers and are released
when each solve returns. Input geometry/state allgathers are NOT removed.

No MPI communication is issued from a threaded QR kernel. Multiple BLAS threads
on workers, if selected externally, only perform numerical linear algebra.

## Surface search

The old query concatenates every visited bin, then sorts and uniquifies all
marker indices. The new helper first deduplicates the at-most-27 canonical bin
keys (important with periodic wrapping), gathers each bin once, and skips marker
sorting when the gathered index ranges are already ordered. This covers the
common single-bin case. Interleaved streams retain optimized `std::sort`:
a tested heap-based merge was slower and was not retained.

The returned candidate vector is sorted/unique and identical to the old one.
The subsequent distance calculations, tie comparator, nearest count, fit radius,
normal averaging, phi values, and geometry assembly are unchanged. No surface
marker thinning or search-radius reduction is used.

## Timings and diagnostics

Enable `weakInterfaceTracePerformanceLog=1` AND the existing MPMSubroutines log
category/solver log level. Turn timing back off after short profiling runs; its
extra timing reductions and log I/O are not free.

The existing `C0 boundary-compatible QR component=` lines retain their fields
and additionally report `worker` and `backend`. Rank zero can log each component
because scalar result metadata is now broadcast too.

Additional per-component times:
- setup: dependency reduction, mass/constraint setup, and dense packing;
- qr: the numerical backend (factor, triangular solves, and reflector application);
- recovery: final physical corrections, impulses/reactions, and residual checks.

Additional rank-max wall times:
- total projection;
- geometry+rows (including fitting);
- fitting (a nested subset, NOT an additional independent stage);
- node-state gathering/deduplication;
- common setup;
- component solves plus distribution;
- commit.

Candidate visits and maximum candidates/query are also reported. Maxima from
different ranks need not sum to the maximum total. The solves+distribution time
includes waiting for the slowest worker; compare it with the individual qr
measurements. Do not infer the whole timestep bottleneck from matrix dimensions.

## Verification supplied

From the GEOS root:

```bash
T=src/coreComponents/physicsSolvers/solidMechanics/tests/mortarPerformance
python "$T/run_tests.py" --output /tmp/mortar-perf-native
python "$T/run_tests.py" --output /tmp/mortar-perf-lapack --lapack --benchmark
python "$T/run_tests.py" --output /tmp/mortar-perf-sanitize --lapack --sanitize
```

The existing boundary-mortar checks are rerun through every backend. Independent
SVD fixtures include tall/wide and rank-deficient systems. The scheduler tests
use an explicitly labeled thread-based collective emulator, including ranks
with no assigned components, empty systems, worker failure, and a concurrency
barrier that detects a secretly serialized scheduling loop.

Actual MPI transport test (on a host with mpicxx/mpiexec):

```bash
python "$T/run_tests.py" --output /tmp/mortar-perf-mpi --mpi --lapack
```

That tests the production scheduler with real MPI at 1/2/3/4/8 processes. It is
still not a full GEOS partition/migration test. If the launcher is unavailable,
the runner records NOT RUN; it never calls an emulator a real MPI test.
For a Slurm-only environment, compile `testMPI.cpp` with the QR source using the
site MPI compiler and launch the resulting executable with `srun -n N`.

Standalone CMake is also supported in that test folder:

```bash
cmake -S "$T" -B /tmp/mortar-cmake -DCMAKE_BUILD_TYPE=Release -DMPM_TEST_LAPACK=ON
cmake --build /tmp/mortar-cmake -j 2
ctest --test-dir /tmp/mortar-cmake --output-on-failure
```

Full PFW pass-through and particle-byte comparison:

```bash
python "$T/testPFWOptions.py" --output /tmp/mortar-pfw-options \
  --baseline-pfw /path/to/unpatched/scripts/preProcessing/materialPointMethodParticleFileWriter
```

Its writer processes use explicitly emulated collectives. The existing 26 PFW
unit checks remain unchanged.

## GEOS acceptance before using results

1. Rebuild the full CPU executable from the patched tree. Verify the resolved
   backend and worker assignments in the log.
2. On the SAME particle file/grid, compare reference/one worker/old search with
   native/one worker/new search over a short matched physical interval. Compare
   row counts, fallback, rank, jump and BC residuals, momentum/work/reactions,
   and material stress histories. Separate init/output from per-step timings.
3. Repeat the fast path with 2/4/8 spatial ranks where the existing grid permits.
   Test an interface spanning a rank boundary and a rank with no local interface.
   Verify owner masses and donor/buffer values are counted once.
4. Compare restart versus uninterrupted execution and check the final state
   consumed by G2P/XPIC, not only the solver's internal residual.
5. Benchmark warmed-up 50-100 step windows with controlled output cadence and
   thread counts. No speedup target is an acceptance substitute for correctness.

Do not change the gradient, geometry, tolerances, or support just to improve the
performance number. Existing contact, higher-order transfer, and CPU/GPU
limitations of boundary-compatible mortar remain unchanged.

## Next iteration: several interfaces

This patch does not split by sphere/interface ID. Three jobs per disconnected
interface can provide useful coarse-grained parallelism, but only when their
rows are independent AFTER boundary/fallback reduction. Geometrically separate
surfaces can share a free matrix-side grid unknown; such entities must remain
in the same algebraic block. The current multiple-sheet support/fallback policy
is unchanged. That future scheduler can reuse the compute / collective-status /
result-distribution separation introduced here.

## Numerical references

The optional backend uses Netlib LAPACK's DGEQP3 (pivoted QR, Level-3 BLAS),
DTRTRS, and DORMQR. The native norm-refresh criterion follows the safeguarded
partial-column-norm strategy documented in DLAQP2 / LAPACK Working Note 176.

- https://netlib.org/lapack/explore-html/d0/dea/group__geqp3_gae96659507d789266660af28c617ef87f.html
- https://netlib.org/lapack/explore-html/dc/db8/group__laqp2_gaf8ebda4d584de0767a5563dc3fd62fb6.html

These references describe the linear algebra, not validation of this GEOS patch.
