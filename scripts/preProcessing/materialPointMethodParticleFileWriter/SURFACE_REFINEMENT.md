# Two-pass surface refinement in PFW

## Configuration

```python
# Existing geometry objects need no changes.
pfw["surfaceRefinementLevel"] = 4
pfw["surfaceRefinementMethod"] = "cartesian"  # Default method.
```

Or, for boundary-conforming centroid/volume initialization of single-point
B-spline particles:

```python
pfw["surfaceRefinementLevel"] = 4
pfw["surfaceRefinementMethod"] = "simplex"
# Set particleType=1 on the corresponding geometry objects.
```

`"simplex-conforming"` is also accepted as an alias for `"simplex"`.
These are writer options, not GEOS solver XML attributes. With the default
`method="cartesian", level=1`, the original no-surface-refinement generation
path is retained. Explicit `method="simplex", level=1` enables plane cutting
without additional Cartesian subdivision of the candidate parents.

The level is a subdivision **factor per active direction**, not a recursion
count. Plane strain subdivides x/y and retains the original z thickness.
An object's existing `particleRefinement` factor is also applied. At a
candidate domain, the maximum material refinement factor in its immediate
neighborhood supplies one common partition for both material sides.
Consequently, level 4 does not mean the same particle count as in the older
surface-flagged-center-only implementation: more parents can be selected, and
cut polygons/polyhedra can produce several simplex particles per subdomain.

Two optional safeguards are available:

```python
pfw["surfaceRefinementFallback"] = "cartesian"  # Default; alternatives: "error".
pfw["surfaceRefinementPlaneTolerance"] = 0.05  # Default; fraction of child diagonal.
```

The tolerance controls allowed disagreement between a local tangent plane and
sampled corner labels. It is not a geometric-error guarantee or a conservation
correction. Strict `fallback="error"` is useful when checking whether a
benchmark actually uses the intended simplex cuts throughout its surface band.

## What the two passes do

The lattice here is the **coarse particle reference-domain lattice**, not the
MPM background grid. Pass 1 queries ordinary coarse centers, recording
occupancy, object/material/contact-group labels, and surface flags. An absent
entry in a slice represents an empty sampled location. Existing particle
properties are staged as ordinary output text for owned occupied parents.

Pass 2 selects a parent when any of its immediate neighbors is surface flagged
or has a different region label. This includes the flagged parent itself, both
sides of an occupied/empty boundary, and both sides of an A/B interface.
The neighborhood includes diagonal neighbors: 8 in 2D or 26 in 3D.
No corner sampling is performed throughout the bulk merely to detect this band.

The passes are pipelined through three x-slices rather than keeping or rereading
a complete provisional particle file. Unchanged bulk parents reuse their staged
record without a second membership/property query. Candidate records are
replaced, never appended alongside their children. Empty-center parents in the
band can now contribute particles. Output IDs are assigned in parent order.

### Cartesian method

Subdivide the selected parent and classify each child center using the existing
object priority. Emit ordinary equal-sized children where material is present.
Unlike the previous method, a parent need not have an interior center to be
visited. This remains midpoint/stair-step volume quadrature, not geometric
clipping. Children in the interior away from the band retain the ordinary
initialization behavior, including material-wise particle refinement.

### Simplex method

Each selected parent is subdivided at the requested level. On each subdomain:

1. Query its center and corners, only within this selected band. Geometry
   queries at corners are shifted one floating-point step toward the subdomain
   interior to respect objects with half-open bounds; clipping uses the exact
   corner coordinates.
2. Obtain a usable surface normal and position offset from a local flagged
   sample, or a neighboring coarse surface marker. Reuse staged surface data
   where available. Re-query at the child scale rather than merely subdividing
   a single parent tangent plane.
3. Use one shared plane for both sides of the cut. Intersect edges analytically,
   clip the square/cube, and triangulate/tetrahedralize the convex pieces.
4. Assign the sides using object/material/contact-group labels. Retain material
   pieces and omit void. Check the simplex partition's volume closure before
   removing the void side; never rescale volumes to conceal a failed cut.
5. Place a particle at each simplex centroid and assign its physical volume.
   Recompute surface-position offsets relative to that particle's centroid.

The cutter operates in a unit square/cube for numerical scaling. Its generated
simplex volumes and first moments reproduce the **adopted plane-cut region**,
not the exact curved object. Positive-volume pieces are retained; numerical
coincidences at vertices/edges are handled without an arbitrary physical
sliver-volume threshold.

## Current GEOS file-format compatibility

This update changes PFW only. It does not add a C++ `particleReferenceGeometry`
field, triangle/tet connectivity, CPDI2 support, or a new particle interpolation.
The cut simplices are temporary initialization geometry.

For `particleType=1` (SinglePointBSpline), the existing reader obtains volume
from the determinant of the R-vectors. The writer therefore emits an
**equal-volume axis-aligned box** at each simplex centroid. Its aspect ratio
matches the local sampling box; in plane strain its z thickness is unchanged.
Thus the existing reader receives the intended physical volume without new
columns or a density correction. Visualization through the existing R-vectors
shows these equivalent boxes, not the temporary simplices.

Emitting an actual simplex-cut particle with another particle type produces a
clear error. Use Cartesian refinement for CPDI/CPTI. Merely giving a CPDI
particle an equal-volume box is not the same as supplying a conforming CPDI
integration domain and is deliberately not advertised as such.

## Parallel PFW and rank-based object slicing

Ownership remains a half-open interval of coarse x rows:

```
start = ceil(number_of_rows * rank / number_of_ranks)
stop  = ceil(number_of_rows * (rank + 1) / number_of_ranks)
```

The implementation uses integer arithmetic for those ceilings. Only owned
parents emit particles. All children inherit their parent's writer ownership;
they are not reassigned according to an equivalent box's extent.

Each rank queries one additional coarse x row on either side as a read-only
halo. The halo supplies occupancy/flags across the rank boundary but creates
no particle records or output. Recomputing this thin halo avoids communicating
Python geometry objects, full particle state, or a global occupancy array.
Small volume/diagnostic reductions and material/type metadata collectives remain.
Ranks with zero owned rows participate in the collective stages normally.

For `make_objects(rankxmin, rankxmax)`, PFW now passes bounds derived from the
actual integer-owned rows, expanded by one coarse row when refinement is active.
This fixes a mismatch between fractional rank ranges and integer row ownership.
The factory must return all **complete objects whose bounding boxes intersect
these bounds**, in the same relative priority order on every rank. It must not
clip or redefine those objects to the requested rank interval. The geometry,
material IDs, and property functions must be deterministic functions of global
position; process-order-dependent random draws are not rank-count invariant.

Refined child queries filter objects using the parent's complete x extent,
not just the parent's center. This is necessary for empty-center cap recovery.
The writer also unions material/type metadata across ranks before rank 0 writes
XML, and reduces particle volume across all rank-owned outputs rather than
printing rank 0's local volume.

Periodic neighbors wrap their lattice indices. Neighboring surface planes are
translated into the local periodic image before use. With periodic x and a
slice-limited factory, the writer passes the **full x interval** to that factory:
a single `(xmin,xmax)` interface cannot express a disjoint opposite-edge halo
safely while preserving object priority. This is a documented object-memory
tradeoff, not a global particle/occupancy gather.

The original code's serial-only Flux writer path is unchanged. The MPI/sliced
path discussed here is the existing non-Flux PFW path.

## Limitations and fallback

- Features with no coarse occupied sample, label change, or nearby surface flag
  can remain undetected. This is a coarse-resolution requirement, not a certified
  search for arbitrary subparticle features.
- One local plane may not describe a corner, thin layer, or multi-region
  junction. Missing/invalid surface data, more than two sampled regions, large
  side-label inconsistencies, invalid partitions, or a simplex centroid in the
  wrong actual region trigger **whole-subdomain** fallback. Individual
  simplices are never filtered out of an otherwise accepted material partition.
- Default fallback is ordinary Cartesian midpoint sampling at that subdomain's
  requested resolution. It can retain the corresponding staircase error.
  Strict mode stops instead. Fallback reasons and counts are reported.
- The same plane partitions A and B inside one subdomain, so those two local
  volume partitions are complementary. Independently chosen planes in adjacent
  subdomains are **not guaranteed to form a globally C0 surface** or to coincide
  with GEOS's reconstructed mortar. This change does not fix the solver's
  topology, global-boundary, or post-projection fallback issues.
- Conservation is relative to the adopted local geometry, not the previous
  staircase file. Recovered material and clipped-away overfill can change total
  initialized mass at fixed physical density. No cell-filling correction is used.
- Additional simplex particles cost preprocessing and subsequent solver time.
  No claim about improved stress accuracy is made without GEOS verification.

## Diagnostics and tests

The writer prints and saves `mpmParticleFile_<case>_refinement.json`, including
candidate domains, recovered empty parents, simplex particles, fallback counts,
geometry-query counts, total particle volume, and per-rank peak slice-cache
sizes. Query counters count the writer's call sites, not additional searches
performed internally by an object's methods. Peaks are maxima across ranks;
counts are sums. More ranks can increase halo query counts without changing
particles.

Run geometry/ownership unit tests without MPI:

```bash
python -m unittest discover -s tests -p 'test_surface_refinement.py' -v
```

Run complete writer regression with actual MPI on an appropriate machine:

```bash
python tests/run_writer_regression.py --output /tmp/pfw-regression \
    --communicator mpi --launcher mpiexec
# On an allocated Slurm node, substitute --launcher srun.
```

The regression compares merged particle files byte-for-byte across 1, 3, and
7 ranks, including empty-rank cases, a material interface, separated free-surface
objects, and a slice-filtered object factory. It validates positive R-vector
volumes, global volume reporting, IDs, and material/type XML coverage.

When an MPI runtime is unavailable, `--communicator emulated` explicitly opts
into a **test-only file-backed collective API**. This exercises the entire
writer in independent processes but is not evidence of actual MPI transport or
cluster performance. The test worker is never imported by production PFW.
