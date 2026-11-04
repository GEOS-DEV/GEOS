# Axisymmetric reconstructed surfaces, nodal area and surface tension

This extends the pressure-regions baseline. Coordinates remain x=r, y=axial,
with no swirl; tensor zz is hoop. All measures include a full revolution.
No particle-file schema, XML spelling or additional particle geometry changes.

## Controls and flow

- `computeParticleSurfaceNormalsAndPositions=1` requests one reconstruction
  before mechanical P2G, then resets to zero. It updates particles with Surface,
  Cohesive or DamagedCohesive flags. Use PFW flags or `surfaceDetection=1`.
- `normalAndPositionMethod='DFGAndVolumeIntegration'` uses a geometric boundary
  normal and volume-matched plane. `LogisticRegression` refines the normal from
  the dominant neighboring field; free surfaces and failed fits use the volume
  normal. Each field retains its own volume-matched offset, including gaps.
- Existing nonzero normals are preserved unless
  `overwriteExistingNormalsAndPositions=1`. Active CZ reference geometry is
  preserved even when overwrite is requested. The normal/offset is transported
  back to reference configuration using the existing F/cofactor transforms.
- A `ReferenceCohesiveZones` event may request `computeNormalsAndPositions=1`.
  Region-level `computeParticleSurfaceNormalsAndPositions=1` schedules a pass
  before that enabled region's reference initialization, using the solver's
  selected method. PFW permits zero initial CZ normals only with a reconstruction
  request; keep ContactGroup, SurfaceFlag, SurfaceNormal, SurfacePosition and
  CZTag columns. Reconstruction must occur by the interface's activation step.
- `computeNodalArea=1` reconstructs current grid planes every step, independently
  of the one-shot particle update. `gridSurfaceArea` is per field, in physical
  area units. The field-0 log sums owner nodes globally (no MPI ghost duplicates).
  A shared interface appears once in each of its participating fields.
- `enableSurfaceTension=1`, `surfaceTensionPairs=[[A,B,gamma],...]` activates
  constant pair coefficients. B=-1 denotes void. Gamma is force/length; phase
  numbers are ContactGroup ids, not material ids. The existing gradient
  tolerance, force sign and capillary timestep limit are retained. Surface
  flags, particle reconstruction and nodal-area output are not prerequisites.

## Geometry and measures

A transient Q1 geometric projection uses current unified quad/triangle geometry
for all mechanical interpolations, including point Q1 and cubic B-splines. It
splits the polygon at grid cells and integrates the ring-weighted Q1 basis
exactly on each piece. This differs from the mechanical CPTI interpolant:
using the latter for color with a Q1 support capacity produces false gradients
near symmetry boundaries even in a fully tiled triangle mesh.

The projection uses the existing temporary integration-domain scaling and axis
fit, with physical particle volume unchanged. A glyph outside the local grid
uses a Q1 point fallback. Large scaling or under-resolved geometry therefore
reduces geometric accuracy. No extra r vectors or per-particle vertices are
stored. One transient nodal scalar per field, `gridSurfaceMaterialVolume`, is
added. Owned particles scatter; shared sums are reduced before reconstruction;
owner planes replace ghost copies before particle gathers.

The boundary integral uses V*(G_r+T,G_y), with T=<N/r>. The radial N/r term
removes the interior circumference gradient. Support is clipped at r=0 and
physical Symmetry faces. Logistic samples use meridional area V/(2*pi*r),
dimensionless coordinates and reflected axis images, avoiding radius and
particle-count bias. The criterion that assigns a particle to a DFG field is
unchanged. The fit is local and has one plane per field; corners, thin features
and multiple competing surfaces need refinement.

For a plane n.(x-xI)=d, solve

    V_I = integral[n.(x-xI)<=d] 2*pi*r*N_I dA

by monotone bisection. The support volume is radius dependent; a plane through
a near-axis node does not generally split it in half. The exact line integral

    A_I = integral[plane] 2*pi*r*N_I ds

uses the same clipped support. The revolution axis has zero surface area.
Nodal planes are converted to particle-relative normal distances before
averaging; averaging node-relative offsets directly introduces a centroid bias.

## Surface stress and curvature

Void color is geometric material volume divided by the corresponding Q1 ring
support volume. Material/material color is VA/(VA+VB), independently on each
existing DFG side. Scalar colors are reflected at symmetry faces for B-spline
buffer-node gathers. The mechanical gradient still uses the particle's actual
mapping; the third slot T is never treated as a color-gradient component.

    tau = gamma * |grad(c)| * (I - n outer n),  n=grad(c)/|grad(c)|
    f_r = -V * (tau_rr*G_r + tau_ry*G_y + tau_tt*T)
    f_y = -V * (tau_ry*G_r + tau_yy*G_y)

The hoop stress supplies circumferential curvature without explicitly dividing
curvature by radius. Full-ring V already supplies 2*pi*r. No z force is created.
The axisymmetric material/void layer is sampled only on the material side, so
its surface stress carries a factor two. Material/material layers sample both
sides and do not receive this factor. This normalizes the resolved flat-layer
limit; finite-radius, cut-cell and diffuse-interface errors remain. The new
normalization is confined to axisymmetry; Cartesian force formulas are retained.

This is the existing diffuse surface-stress formulation with cylindrical work
terms, not an exactly balanced Laplace-pressure discretization. Check radius,
area, residual velocity and pressure under refinement. Constant pair coefficients
are supported; no prescribed wall contact angle or Marangoni law is added.

## Examples and checks

`examples/axisymmetricSurfaces` provides a constrained-length cylinder and a
free sphere with initial Laplace pressure gamma/R and 2*gamma/R respectively.
The default cylinder uses CPDI core/CPTI skin and initially zero surface normals.
The sphere's mixed mode uses the existing curved-surface refinement. CPDI, CPTI
and B-spline comparison options, a runner and an area/particle-volume report are
included. Set AXIS_SURFACE_TENSION=0 for geometry-only diagnostics.

From the PFW directory:

```sh
ASAN_OPTIONS=detect_leaks=0 python3 tests/run_axisymmetric_surfaces_cpp_tests.py
python3 tests/run_axisymmetric_surfaces_writer_tests.py --output /tmp/axis-surfaces
```

Tests compile the production mapping, geometric helpers and surface-stress
assembly excerpt with strict warnings and ASan/UBSan. They cover volume/area
inversion, uniform quad/triangle tilings, symmetry faces and axis, regression
scale/refinement invariance, particle-relative offsets, cylindrical area and
capillary-work refinement for cylinders/spheres. Writer tests execute complete
PFW emission at 1/3 emulated ranks and compare particle files byte for byte.
These do not replace a full GEOS build, dynamics convergence, actual MPI/restart
or GPU execution. Existing two-group CZ restrictions, CZ-tagged B-spline block,
separability criteria, SPH/overlap and generic PFW traction limitations remain.
