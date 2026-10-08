# PFW two-pass surface refinement update

Baseline: the files actually supplied in `2609161200_geos_mpm_bundle.zip`.
The bundle identifies review commit `271e3edc9553080f388395d86df7b0c02f68a4ee`
and includes a worktree overlay. `BASELINE_AND_FILES.json` records the exact
archive/file hashes; the patch targets those file contents, not an assumed
clean checkout of that commit.

## Install

From the GEOS repository root:

```bash
git apply --check /path/to/pfw_two_pass_surface_refinement.patch
git apply /path/to/pfw_two_pass_surface_refinement.patch
```

Alternatively, copy the contents of `updated_files/` over the corresponding
GEOS repository paths. Use one installation method, not both. The complete
updated `particleFileWriter.py` and both new helper modules must be installed
together. Existing runClean staging copies `pfw_*.py`, so it also picks up the
new modules. No geometry-object or C++ solver file is changed.

## Use

```python
pfw["surfaceRefinementMethod"] = "cartesian"
pfw["surfaceRefinementLevel"] = 4
```

Or use `surfaceRefinementMethod="simplex"` with `particleType=1` on objects
whose cut domains will emit simplex-centroid particles. The local simplices
supply position and volume; equal-volume R-vector boxes preserve compatibility
with the existing particle reader. This does not add triangle/tet reference
geometry to the solver or turn CPDI into a simplex method.

Both methods use a band selected from ordinary coarse-center occupancy,
surface flags, and region-label changes. Empty-center neighboring domains can
be refined. Interior records are reused, and selected coarse records are
replaced. Three cached coarse x slices bound the preprocessing storage.

Rank-sliced generation uses integer parent-row ownership, a read-only one-row
halo, enlarged `make_objects` bounds, periodic marker shifts, global
material/type metadata, and globally reduced particle volume. Each parent is
emitted by exactly one rank. Geometry and state getters must be deterministic
and object factories must preserve complete geometry and priority order.

See `updated_files/scripts/preProcessing/materialPointMethodParticleFileWriter/
SURFACE_REFINEMENT.md` for configuration, implementation details, factory
requirements, diagnostics, limitations, and executable tests.

## Validation performed

- All 18 geometry, detection, ownership, periodicity, and conservation unit tests passed.
- All 24 complete-writer regression configurations passed. Merged particle files
  were byte-identical for 1, 3, and 7 independently executing processes, including
  empty ranks and a rank-filtered object factory.
- The complete-writer multi-process tests used an explicitly emulated file-backed
  collective API because an MPI runtime was unavailable. **They are not actual
  MPI transport or cluster-performance tests.** A real-MPI regression command is
  included in the guide.
- The supplied 3D workshop input also completed PFW generation with both methods
  at reduced resolution. No-refinement output matched the original writer byte
  for byte in that smoke case.
- Python compilation, patch whitespace checks, and forward `git apply --check`
  plus `git apply` against the exact baseline passed.

No GEOS mechanics simulation, solver build, MPI cluster benchmark, or stress
accuracy comparison was run. The simplex approximation is local: it does not
guarantee a globally C0 surface or solve the existing mortar/fallback issues.
Unresolved local cuts have counted Cartesian fallback (or explicit strict-mode
failure), and subparticle features entirely missed by the coarse pass can remain
undetected. See `test_results/TEST_RESULTS.md` and the accompanying logs.
