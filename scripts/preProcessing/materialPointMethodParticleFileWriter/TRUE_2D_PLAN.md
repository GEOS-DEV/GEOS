# Native Cartesian 2D: implementation and integration plan

This incremental patch targets the branch with unified geometry, geometry-aware
splitting, CPTI, and the hex-bulk/tet-surface workflow already applied. It adds
native XY plane strain to the same timestep solver as Cartesian 3D. Axisymmetry
and plane stress are separate future formulations. Rebuild GEOS and regenerate
particle/XML inputs. Old Python PFW inputs are normalized; old GEOS inputs and
checkpoints are not a compatibility target.

## Implemented architecture

1. **Select the formulation once.** `SolidMechanics_MPM formulation="PlaneStrain"`
   requires an `InternalMesh spatialDimension="2"`. Cartesian3D remains the
   default. The shared timestep dispatches dimension-specialized mapping,
   gather/update, integration-domain, and subdivision operations.
2. **Use a native mesh and communication topology.** Plane strain generates XY
   C2D4 cells, four-node support, two-node faces, and one stored z coordinate
   equal to zero. There are no z cells, z ownership searches, z ghost slabs, or
   z Cartesian neighbors. MPI Cartesian topology has two dimensions. XY
   partitioning, periodicity, and the four in-plane boundaries remain active.
3. **Store compact material geometry.** G and G0 are 2x2 for all native 2D
   particle types, and 3x3 in 3D. Positions and framework vectors retain three
   components with inactive z zero. Constitutive tensors remain 3x3: Fzz=1,
   Lzz=0, and sigma_zz remains available for plane-strain material response.
4. **Use one physical thickness.** A finite positive
   `InternalMesh planeStrainThickness` (default 1) is stored on the domain
   partition and restored on restart. PFW writes this from the matching scalar.
   V=A*t, mass=rho*A*t, and an interface line represents area=length*t.
   In-plane diameter, CPDI scaling, splitting, and particle characteristic
   length use XY geometry; particle fracture length uses sqrt(V/t).
5. **Generate and split native shapes.** Shape 3 is a parallelogram; shape 4
   is a triangle. A point or B-spline can carry either shape. CPDI requires
   parallelograms; CPTI requires triangles. Mixed quad-bulk/triangle-surface
   generation uses the existing explicit simplex/native/cut-type settings.
   The named hexTet preset retains its 3D-only meaning.
6. **Write native output.** Silo writes spatial dimension 2 and quad/triangle
   zones. VTK writes quad/triangle cells with z=0 coordinates. Grid quads have
   consistent perimeter connectivity. The shared VisIt renderer selects a
   fitted View2D for native-2D metadata; the supplied case renderers also handle
   four-value XY extents. No artificial extrusion is generated for display.

## Particle data and mapping

Generators are rows. With column edge matrix E=G^T, vertices are reconstructed
from the centroid; neither vertex arrays nor half-edge arrays are persistent.

| Geometry | Tag | Stored G and G0 | Reconstructed vertices | Geometric measure |
|---|---:|---|---|---|
| 2D parallelogram | 3 | Two full XY edges each | c +/- g1/2 +/- g2/2 | det(G) |
| 2D triangle | 4 | Two full XY edges each | v0=c-(g1+g2)/3; v1=v0+g1; v2=v0+g2 | det(G)/2 |
| 3D box | 0 | Three full XYZ edges each | c +/- g1/2 +/- g2/2 +/- g3/2 | det(G) |
| 3D wedge | 1 | Two base edges and one thickness edge each | Base triangle offsets +/- g3/2 | det(G)/2 |
| 3D tetrahedron | 2 | Three full XYZ edges each | v0=c-(g1+g2+g3)/4; va=v0+ga | det(G)/6 |

The two matrices occupy 64 bytes in 2D, 144 bytes in 3D, plus the tag. Physical
particle volume is a separate field; it equals measure*t initially in 2D.

| Particle type | 2D support | 3D support | Domain restriction |
|---|---:|---:|---|
| SinglePoint (0) | 4 | 8 | Any supported shape of the selected dimension |
| SinglePointBSpline (1) | 16 | 64 | Any supported shape of the selected dimension |
| CPDI (2) | 16 raw | 64 raw | Parallelogram / parallelepiped |
| CPTI (3) | 12 raw | 32 raw | Triangle / tetrahedron |

Raw domain entries include repeated grid nodes; effective mapping coalesces
those entries. Point bases are bilinear/trilinear or tensor-product cubic
B-splines. CPDI averages the 2^d corner basis values and forms gradient
coefficients from the reciprocal full-edge basis. CPTI averages d+1 corner
values, with barycentric gradients b_a=E^(-T)e_a and b_0=-sum(b_a):
`S_I=sum_a N_I(x_a)/(d+1)` and `grad S_I=sum_a N_I(x_a)*b_a`.
Both positions and gradients use the same temporary integration geometry.

Homel and VP-Hencky retain their box projection policies in the active dimension.
CPTI uses uniform diameter clipping: Homel limit 0.99998*min(active h), or
VP-Hencky selection limit 0.99999*neighborRadius. It does not use the box
volume-preserving projection. Integration scaling changes neither stored G/G0
nor physical volume/mass/F. G advances by the physical increment. FOnly preserves
G; FAndGeometry rebuilds G=G0*Fnew^T using the active block. Splitting applies the
same child map to G and G0, conserves physical mass/volume and weighted centroid,
and retains shape/type. Cohesive particles remain excluded from generic splitting.

## Input and file changes

```python
pfw['formulation'] = 'PlaneStrain'
pfw['planeStrainThickness'] = 0.01  # physical length, not a grid cell count
pfw['xmin'], pfw['xmax'] = 0.0, 1.0
pfw['ymin'], pfw['ymax'] = 0.0, 1.0
pfw['nI'], pfw['nJ'] = 42, 42
pfw['periodic'] = [False, False]
pfw['objects'] = [geom.disk('sample', [0.5, 0.5], 0.3, particleType=2)]
# Optional CPDI bulk / CPTI surface triangles:
pfw['particleGeometry'] = 'fromRefinement'
pfw['surfaceRefinementMethod'] = 'simplex'
pfw['surfaceRefinementLevel'] = 2
pfw['simplexParticleType'] = 3
```

The complete, non-submitting example is `examples/true2D/pfw_input_true_2d.py`.
Native XML has xCoords/yCoords and nx/ny, C2D4, spatialDimension=2, and the single
mesh thickness input. It has no zCoords/nz. The solver selects PlaneStrain;
planeStrain is no longer an XML input. Native particle rows use GeometryVersion=2,
shape 3/4, ParticleVolume and GeometryXX/XY/YX/YY. PositionZ is written as zero;
the native reader permits it to be omitted. Restart geometry version is 3,
which distinguishes compact native storage from previous checkpoints.

Old `planeStrain=1` PFW input maps to PlaneStrain. If no explicit thickness is
supplied, its old zmax-zmin supplies t; otherwise new inputs default to t=1.
Old nK/ppcz do not create z layers. zpar must be 1. Old planeStrain inputs suppress an inactive z-periodic flag;
a new native formulation input rejects z periodicity.
Conflicting formulations and nonpositive/nonfinite thickness are rejected.
Old 6-face BC inputs are accepted and inactive z faces are cleared; native
boundaryConditionTypes/padSides may use four entries. Boundary-type tables may
use time plus four entries; old six-face tables are accepted. F tables retain
[time,Fxx,Fyy,Fzz], with Fzz required to equal 1. The optional transverse-velocity
matrix retains its existing 6x2 interface; only in-plane face/components act.

## Objects, strength statistics and output

See `PFW_2D_OBJECTS.md` for the audited example objects, new disk/annulus and
rectangle APIs, the legacy extrusion adapter, and the limits of general 3D
cross-section objects. The 11 plane-strain example inputs now omit artificial
z mesh extents. Brazilian-disk and borehole Weibull fields use dim=2 with
physical sample volumes area*t. Their random realizations change from the old
3D fields; preserving an old realization is not promised.

Built-in native Weibull wrappers apply t to their raw Voronoi/layer areas before
the physical minimum-volume clamp. They retain the sampled random draws and
support idempotent rebinding. A custom field may implement
`setPlaneStrainThickness(t)`; arbitrary user StrengthScale values are not rescaled.
Existing 3D stochastic wrappers retain their 3D meaning. The domain scalar is
available for fracture regularization; this patch changes the particle-sized
lengthScale field to sqrt(V/t), without overwriting a scalar constitutive
lengthScale parameter.

Type-B CZ anchors intersect the actual XY quad/triangle. Native line-interface
nodal weights integrate the Q1 basis along the segment and multiply by t.
Contact/nodal physical-volume normalizations use hx*hy*t, without the old
thin-layer half factors. Force and energy therefore retain physical 3D units;
area/geometry and plot extents do not include an invented z cell. Constitutive
output retains full tensors including sigma_zz.

## Verification and target acceptance

Local tests compile production geometry/mapping headers and selected full
production bodies with small CPU adapters. They cover 2D/3D mapping moments,
compact storage, resets, tiny/large scales, native subdivision conservation,
Q1 line integrals, quad topology, and partition-neighbor enumeration. Writer
checks use 1/3/7 explicitly emulated ranks, including empty ranks. Object tests
vary query z and legacy extrusion depth and check flags, normals, surface
positions, and thickness-dependent Weibull scales. VisIt tests mock the camera
API; they do not run VisIt. See the delivery validation.json for exact results.

The complete GEOS dependency tree is unavailable here. Required target checks:

1. Full GEOS build with project warnings and the desired CPU/GPU configuration.
2. Native quad-grid initialization and translation/affine-deformation cases for
   point, B-spline, CPDI, CPTI and mixed particles, then PIC/FLIP/XPIC/FMPM.
3. Compare forces/mass/energy for two thicknesses, retaining the same in-plane
   geometry. Check plane-strain constitutive stress including sigma_zz.
4. Real MPI periodic/nonperiodic ownership, migration, ghosts and restart
   round trip; compact geometry must survive all of them.
5. Splitting during motion, rigid-body jamming, contact and coupled CZ mechanics,
   followed by Silo/VTK readback and the actual VisIt renderer.

No full GEOS, real MPI, GPU, restart round trip, or coupled contact/CZ run is
claimed by the standalone checks. Axisymmetry, plane stress, CPDI2, shared-facet
cohesive subdivision and runtime conversion of exposed quad/hex surfaces remain
separate work.
