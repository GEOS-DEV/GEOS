# Native particle reference geometry for MPM output

## Scope

Native triangular-wedge/tetrahedral output preserves the cut-simplex initialization
geometry while keeping SinglePointBSpline mechanics based on the existing
equivalent boxes. It does not introduce CPTI, CPDI2, or a connected element mesh.

## Enable

In the PFW input:

```python
pfw["surfaceRefinementMethod"] = "simplex"
pfw["surfaceRefinementLevel"] = 4
pfw["particleReferenceGeometry"] = "fromRefinement"
pfw["particleDomainOutput"] = "native"       # default; optional
```

Geometry objects producing cut particles must use `particleType=1`
(`SinglePointBSpline`). A cut triangle is extruded through the existing plane-strain
reference thickness. A 3D cut simplex is a tetrahedron. Uncut and fallback pieces
remain boxes. Cartesian refinement also works with metadata enabled, but produces
box geometry. Selecting `fromRefinement` does not convert every bulk particle to a
simplex or add any object/surface queries.

The default `particleReferenceGeometry="equivalentBox"` writes the old particle
file with no added columns. Rerun PFW to obtain native geometry: the discarded
simplex shapes cannot be recovered uniquely from an old equivalent-box file.

To inspect the original numerical domains while retaining the new metadata:

```python
pfw["particleDomainOutput"] = "equivalentBox"
```

The corresponding XML option belongs on each `ParticleRegion`, NOT on the solver:

```xml
<ParticleRegion ... particleDomainOutput="equivalentBox"/>
```

PFW writes this attribute when native metadata is enabled. Old XML/files are valid;
new files with geometry columns require the updated reader. The output mesh name
and one-value-per-particle field association are unchanged.

## Build compatibility

Native geometry files require the patched reader and output writers. Build GEOS
with the desired Silo/VTK backend. Full solver compilation and runtime output
validation require the dependency submodules absent from the supplied archive.

## Data representation

The new particle-file columns are appended after all existing columns:

```
ParticleReferenceGeometry
GeometryRVectorXX GeometryRVectorXY GeometryRVectorXZ
GeometryRVectorYX GeometryRVectorYY GeometryRVectorYZ
GeometryRVectorZX GeometryRVectorZY GeometryRVectorZZ
```

All ten columns must be present together. The importer validates finite data,
geometry tags, positive initial orientation, allowed interpolation type, and
agreement with the INITIAL mechanical volume (relative tolerance 1e-9). It never
rescales geometry, volume, mass, or density to repair disagreement.

| Tag | Geometry | Three ROW vectors | Signed reference volume |
|---|---|---|---|
| 0 | Parallelepiped | Existing half edges | 8 det(G) |
| 1 | TriangleWedge | Two full triangle edges and full z extrusion | det(G)/2 |
| 2 | Tetrahedron | Three full edges from one vertex | det(G)/6 |

The particle position is the centroid. For a tetrahedron, vertex zero is the
center minus one quarter of the vector sum. For a wedge, the first bottom vertex
is the center minus one third of the first two vectors and half the third vector.
Plane-strain wedges require in-plane triangle edges and z extrusion; tetrahedra
require a 3D solver configuration. Unsupported combinations fail explicitly.

The existing `particleRVectors` and `particleReferenceRVectors` retain their
original equivalent-box meaning. They remain authoritative for legacy volume
initialization, B-spline support, contact, CZ TypeB, scaling and other mechanics.
In particular, `numberOfVerticesPerParticle()` is NOT changed: its existing value
also sizes interpolation mappings and is not a visualization-corner count.

New state is registered as:
- `particleReferenceGeometry`: one integer per particle, also available as a plot field.
- `particleReferenceGeometryVectors`: nine doubles per particle, restart data but not plotted by default.

The optional dense arrays are allocated for all particles in geometry-enabled
subregions, including ordinary boxes. Cost is 72 bytes plus an integer per particle
per allocated memory copy, plus allocator overhead. Old files allocate neither array.
For tag 0, output deliberately uses the live legacy box R-vectors; its supplemental
input vectors are not a second box-deformation/scaling policy.

## Motion and reference changes

For simplex particles the current geometry vectors are `G_current = G_ref F^T`
in row-vector storage, centered on the current particle position. There is no new
current-geometry tensor, per-step corner integration, or mechanical equation.

A numerical reference change from F_old to F_new updates ONLY the new generators:

```
G_ref_new = G_ref_old F_old^T F_new^(-T)
```

This preserves current native vertices through the supplied solver's general
F reset, initial temperature/pressure reference change, polymer healing, and
crystal-healing volumetric reference change. Singular/nonfinite rebase operations
fail rather than silently destroying the shape. Physical deformation updates and
`TransformParticles` rotations are not canceled by rebasing.

The native shape can remain distorted while an equivalent numerical box is reset
or regularized. That distinction is intentional. It is not forcibly made to have
the mechanical volume after a reference change. The INITIAL volume equality check
is not a constraint on every later output snapshot.

Runtime subdivision is guarded: when a simplex-tagged particle would actually be
split, execution reports that native simplex subdivision is not implemented. It
does not clone the original simplex geometry into every child. Changing only `particleDomainOutput="equivalentBox"` does not
remove the guard. To permit runtime box splitting, regenerate PFW files with
`particleReferenceGeometry="equivalentBox"`. Offline PFW
refinement remains supported, as does the existing runtime subdivision of boxes.

## Output and ordering

Each particle is exactly one output zone: hex, prism/wedge, or tetrahedron. Particle
and field ordering remains subregion order followed by particle order. Silo uses
consecutive runs of equal zone shape, allowing repeated runs, rather than sorting
particles by shape. Ghost flags remain determined by the particle owner rank.
Vertices stay disconnected between particles; no global ownership or node welding
is introduced. Both native points and existing F/R-vector fields are synchronized
to host for output.

Geometry vertex order is independent of backend order. VTK linear wedges through
9.6.2 require the first triangle normal to point AWAY from the other triangle;
positive reference generators point toward it, so wedge connectivity reverses
both triangles. Silo has a further prism/tetrahedron permutation. The helper makes
these conversions explicit, and the standalone C++ geometry tests exercise the vertex mappings.
A future VTK release changing its wedge orientation convention would need the
small backend permutation updated; it must not change reference geometry or mass.

This displays the local cut domains, not a guaranteed global conforming mesh.
Neighboring local cuts may already differ, and independent particle motion can
produce gaps/overlaps. Visualization does not hide or repair them.

## MPI and PFW rank slicing

No new partition or halo algorithm is added to the two-pass writer. Geometry
columns travel on the SAME record as their particle; there is no sidecar keyed by
rank-local IDs. Merged IDs and existing global type/material metadata are unchanged.
Every rank uses the same column schema, including ranks with no owned coarse rows.

The GEOS importer enables optional arrays from the file HEADER on every block,
including empty local blocks. Subregion copying explicitly transfers the arrays;
they are parent-sized registered wrappers, using existing packing, erasure and
restart machinery. Output vertices have no cross-rank connectivity and need no
additional grid halos. Full GEOS migration/restart and device execution remain
integration checks, not claimed results of the standalone tests.

## Tests

From the repository root:

```bash
PFW=scripts/preProcessing/materialPointMethodParticleFileWriter
python -m unittest discover -s "$PFW/tests" -p 'test_*.py' -v
python "$PFW/tests/run_reference_geometry_cpp_tests.py"
python "$PFW/tests/run_writer_regression.py" --output /tmp/pfw-native-check   --reference-geometry fromRefinement --communicator mpi
```

Use `--launcher srun` on an appropriate allocation. Without MPI, explicitly choose
`--communicator emulated`. Emulated tests launch independent processes with a
file-backed collective API; they are NOT MPI transport/performance tests.

The C++ math driver compiles the production header with a minimal DataTypes stub;
it is NOT a complete GEOS compilation. See documentation/AUDIT.md for the scope of the source and documentation checks. Full Silo/VTK runtime output is not validated by the standalone tests.

Before relying on this in production, build the complete CPU solver, inspect a
native output at t=0 and after deformation, compare mechanical fields with metadata
off/on, trigger an F reset, and test MPI migration plus restart across a changed
rank count. Test a CZ case to establish that the untouched mechanical R-vector
path gives unchanged forces. Native geometry itself does not make CZ/contact
simplex-aware. GPU compilation, Silo runtime I/O, and full GEOS tests were not
available in the patch-preparation environment.
