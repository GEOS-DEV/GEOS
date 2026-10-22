# Geometry-aware particle splitting

Splitting operates on unified current/reference geometry. Native plane strain
uses compact parallelogram/triangle domains; 3D retains boxes, tets and wedges.

## Controls

- `subdivideParticles=-1`: automatic box CPDI VP-Hencky policy only.
- `subdivideParticles=0`: disable general splitting; oversized splitting can
  still operate independently.
- `subdivideParticles=1`: existing box criteria plus native triangle/tet/wedge diameter
  threshold `0.99998 * min(active grid spacings)`. Native selection does not
  depend on the box scaling/reset algorithms.
- `oversizedParticleTreatment='split'`: apply the configured full-diameter / grid
  diagonal ratio to boxes, parallelograms, triangles, tets, and wedges, including point/B-spline particles.
- `oversizedParticleResetRankParticleCountThreshold=-1`: disable rank-local
  oversized-reset fallback. A nonnegative threshold uses the existing strict
  greater-than comparison against start-of-step active rank count. Independent
  general splitting is unaffected. FOnly fallback does not shrink material G.
- `applyVelocityGradientWhenSplittingParticles=1`: affine child velocity; 0
  preserves parent velocity. Kinetic energy is recomputed in either case.

No new file columns, XML controls, geometry storage, or interpolation types.
Old pfw_input aliases continue to be normalized. Use the complete
`examples/geometrySplitting/pfw_input_geometry_splitting.py` input.

## Geometry and bounded work

Boxes keep their previous generator-mask selection and child order. A tet
bisects its longest current edge; the two children share its midpoint and
opposite vertices. Wedges split the longest triangle edge or thickness in 3D.
Native 2D triangles bisect their longest XY edge and parallelograms split
generators. Children retain their parent's shape and interpolation type. Native
2D uses 2x2 G/G0; there is no through-thickness geometry to subdivide.

Native plans repeatedly split the largest violating leaf, bounded to four
children per original parent in plane strain or eight in 3D. Child fractions
may be unequal. A single tet bisection need not reduce maximum diameter. The
log reports additional particles, split tet/wedge parents, and selected parents
with children still above the diameter limit. Remaining refinement can occur
on later steps. Box VP-Hencky measure feasibility stays box-specific.

Each temporary map has coefficients A, centroid q, and fraction w. With rows
as generators, G_child=A G and G0_child=A G0. Current and reference offsets are
G^T q and G0^T q. Current geometry chooses the topology, and the same map acts
on both frames, independently of constitutive F. This works after FOnly and
commutes with subsequent physical convection or FAndGeometry reset.

## Transfer and lifecycle

Mass and current/reference physical volumes use w. Geometry measure is not
substituted for physical volume. Intensive fields (including specific internal
energy), F and constitutive history are copied. Particle-sized wrappers include
meltFlag and imported fields without a particle prefix; infrastructure maps are
excluded. One child keeps the parent ID, while additional IDs use global maximum
and rank offsets. No persistent refinement tree or additional geometry fields.

The optional velocity correction is v_child=v_parent+L*dx. The weighted centroid
identity conserves linear momentum. It increases point kinetic energy by
0.5*sum(m_child*|L*dx|^2), and does not generally conserve point angular momentum.
There is no energy renormalization. With correction off, point kinetic energy
is preserved. Each child's kinetic-energy field is recomputed immediately.

SurfaceTraction enters grid assembly as force: both current and reference forces
are partitioned by w. Ordinary surface offsets preserve the existing approximate
absolute surface point; reference offsets use F^-1*dx rather than G0 coordinates.
Nonzero estimated offsets are rebased; zero retains the existing unavailable
sentinel. Surface force is conserved, but no exact exterior-facet ownership or
contact-quadrature invariance is claimed. Cohesive/damaged-cohesive and active-CZ
particles remain excluded. CPTI tetrahedra now use this same native splitter;
interpolation type 3 is inherited by every child. CZ subdivision and axisymmetry remain separate work. See [CPTI.md](CPTI.md).

Splitting occurs in end-of-step cleanup, after that step's direct history output.
It precedes out-of-range deletion, MPI repartitioning and F reset. Current mapping
and neighbor caches are cleared and rebuilt by the next step's normal passes.
CPDI diagnostics are refreshed before reset selection. Plans exist only for
selected parents on the host; all ranks participate in counts and ID allocation.

## Verification

From the PFW directory:

```bash
python3 -m unittest discover -s tests -p 'test_*.py'
python3 tests/run_subdivision_cpp_tests.py
python3 tests/run_subdivision_cpp_tests.py --sanitize
python3 tests/run_reference_geometry_cpp_tests.py
python3 tests/run_cpdi_geometry_cpp_tests.py
```

The subdivision suite uses the production header and extracts the complete
production `subdivideParticles` body plus its selection helpers. Host adapters
provide views, wrapper resizing/copying, tensor operations and collectives;
these do not replace full GEOS/LvArray/MPI verification. The suite checks 2100
random geometry plans, point ownership (no gaps/overlap), orientation, weighted
centroids, reference/current correspondence, both reset policies, physical
convection, all simplex edges, box masks/order and repeated refinement. Solver
cases cover mixed shapes, unequal child fractions, mass/volume/linear momentum,
energy modes, force/offset/history/melt transfer, deleted/ghost/CZ exclusions,
automatic/explicit/oversized policies, CPDI feasibility, rank fallback, empty
ranks and global-ID offsets. Compiler warnings include Wshadow and Wfloat-equal.
If LeakSanitizer cannot inspect /proc in a restricted host, use
`ASAN_OPTIONS=detect_leaks=0` for address/undefined sanitizers; no leak-check claim.

Required integration checks on a configured GEOS build:

1. Run the supplied mixed box/tet example; inspect child glyphs and split logs.
2. Repeat with plane-strain surface wedges; thickness must remain unchanged.
3. Run 1 and multiple MPI ranks, including an initially empty rank and children
   crossing partition/periodic boundaries; compare mass/volume/linear momentum.
4. Write/reload a post-split restart; compare against uninterrupted execution.
5. Run after both reset policies with constitutive history and surface loading.
6. Exercise the same cases on the configured GPU backend.

Different particle quadrature after splitting can change the numerical solution;
conservation at child creation is not a guarantee of identical grid forces or
contact response. Full GEOS, MPI restart/migration, and GPU checks require the
complete build dependencies and are not established by standalone tests.
