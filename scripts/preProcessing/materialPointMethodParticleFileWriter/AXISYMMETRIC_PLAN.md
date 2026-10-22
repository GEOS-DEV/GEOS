# Axisymmetric MPM: implementation and acceptance plan

This patch is the first, single-field stage, against the branch after rigid-body
removal. ParticleColor remains available. The existing timestep dispatches 2D
kernels; there is no second timestep solver. This is not yet a cylindrical
replacement for every plane-strain option.

## Coordinates, geometry and units

Choose `formulation='Axisymmetric'` in PFW. x is radius r, y is axial position,
and the coordinate/velocity z slot is zero (no swirl). The physical radial
bounds must be nonnegative and particle centers strictly positive. The mesh
retains its ordinary exterior buffer cells; a buffer may have negative r,
but physical particle quadrature does not use negative-radius support.
PFW emits `spatialDimension="2" axisymmetric="1"` on InternalMesh. The domain
partition retains that flag on restart and the solver checks agreement.

Store the existing 2x2 G and G0 for parallelograms or triangles. Their center
is the meridional area centroid. Do not add vertices, hoop generators, or an
extra particle radius field. F, L and stress remain three-dimensional in the
(r, axial, theta) basis: zz means hoop, not a third geometric coordinate.

All extensive measures are for a complete ring: V=2*pi*r*A and m=rho*V.
There is no axisymmetric thickness input. Domain volume is
pi*(rmax^2-rmin^2)*(ymax-ymin). Particle-sized crack-band length uses
sqrt(V/(2*pi*r)); explicit scalar constitutive lengths remain user choices.
The physical volume remains independently integrated, as in the baseline;
F-only resets do not redefine it from the current glyph geometry.

## Implemented in this patch

1. Read/write and validate full-ring particle volumes using the unchanged
   schema-2 native quad/triangle record. PFW keeps old plane-strain and 3D inputs.
   Both ordinary and two-pass surface-refinement emission use each final
   child's own centroid radius, including CPTI triangle tessellation.
2. Radius-weighted Q1-point, CPDI and CPTI operators in the common mapping
   routine, covering cached and on-the-fly mappings. Gradient entry 2 stores
   T=<N/r>, not a z derivative. Coalescing retains it.
3. Radial hoop force and hoop strain rate in FLIP, PIC, XPIC and FMPM; hoop
   strain enters F, physical volume and stress power. Constitutive kernels
   already accept complete 3x3 tensors. Use constitutive symmetries/material
   orientations compatible with no swirl (no r-theta or axial-theta shear).
4. Symmetry at physical r=0. Radial periodicity is rejected. A CPDI/CPTI
   integration domain is uniformly shrunk about its center when necessary to
   avoid negative r; G, G0, mass and volume remain unchanged. A center reaching
   r<=0 is an error, not a clamped trajectory.
5. Ring-weighted subdivision: child fraction is areaFraction*rChild/rParent.
   Mass, current/reference volumes and integrated traction inheritance use
   that fraction. Optional affine child velocities subtract their weighted
   mean offset to conserve momentum. As in the existing splitter, resolving
   an affine velocity gradient can add resolved kinetic energy; it is not a
   total-energy-conserving remap. Use the existing flag to disable it if needed.
   A physical domain that would split into nonpositive-radius children fails.
6. F-only and F-plus-geometry resets keep their existing meanings. Automatic
   planar reset projection uses J/F_theta, preserves the hoop stretch and
   reconstructs only the 2x2 geometry when requested.
7. Existing Silo/VTK native quad/triangle output remains a meridional section.
   Volumes, masses, body forces and reactions have full-ring units. Display
   stress zz as hoop stress. Reaction-history length_z is zero and the x/y
   reactions are full-ring forces; Cartesian force-to-pressure conversions
   must not be reused. The patch does not generate a revolved 3D mesh.

For area A and centroid radius r, the operators are
S_i=integral(r*N_i dA)/(A*r), G_i=integral(r*grad(N_i) dA)/(A*r),
and T_i=integral(N_i dA)/(A*r). In CPDI/CPTI these integrate the material-domain
interpolant of the grid basis. Thus f_r=-V*(sigma_rr*G_r+sigma_ry*G_y+sigma_tt*T),
f_y=-V*(sigma_ry*G_r+sigma_yy*G_y), and L_tt=sum(v_ir*T_i).

For CPDI corner signs s,t and full edges g1,g2, the scalar corner coefficient
is 1/4+(s*g1r+t*g2r)/(24*r). The reciprocal-edge gradient coefficient is
s*b1*(1/2+t*g2r/(12*r))+t*b2*(1/2+s*g1r/(12*r)); the hoop coefficient is 1/(4*r).
For CPTI, scalar corner coefficient is 1/3+(rCorner-r)/(12*r), gradients are
constant barycentric gradients, and hoop coefficient is 1/(3*r).
Shape-weighted nodal position is the volume average, not generally the area
centroid. Tests therefore check cylindrical moments rather than imposing the
Cartesian centroid identity on S.

The cylindrical weak form follows Nairn and Guilkey (2015),
*Axisymmetric form of the generalized interpolation material point method*,
DOI: https://doi.org/10.1002/nme.4792 . Full-ring measures multiply the
per-radian convention by 2*pi consistently.

## Explicit first-stage limits and next patches

- **Contact and CZ:** one contact group (0), no DFG or cohesive particles.
  Next implement meridional-area contact-normal reconstruction, radius-weighted
  interface integrals, axis nodal support volumes and cohesive tractions/work;
  then test annular compression and separation. ParticleColor is independent.
- **Surface tension and tractions:** surface tension, reconstructed surface
  positions/nodal-area output, borehole pressure, spatially restricted confining
  pressure and PFW surface traction generation are rejected. Uniform
  `ConfiningPressure uniformPressure="1"` is supported by the jacketed-sample
  follow-up patch, using background stress and the existing hoop operator. Next integrate 2*pi*r on each
  physical surface segment and add circumferential curvature. Direct GEOS
  particle tractions, when supplied, are already integrated full-ring forces;
  no automatic per-area or per-radian conversion is performed.
- **Weibull:** built-in area-only native Weibull objects are rejected.
  Statistical volumes require 2*pi*integral(r dA), evaluated over each actual
  statistical cell/layer before minimum-volume clamping. Custom wrappers may
  implement setAxisymmetric(); explicit StrengthScale values are unchanged.
  Weibull and localized hoop-driven fracture are deferred: circumferential
  localization breaks the assumed axisymmetry and needs a separate treatment.
- **Other controls:** no swirl, B-spline axis extension, SPH/overlap correction,
  crack-tip detection, profileHistory slice-volume normalization, prescribed BC/F tables or stress-controlled domain motion.
  Boundaries are Outflow or Symmetry. Existing point/CPDI/CPTI body forces,
  initial velocity and compatible stress/thermal initialization remain usable.
  Most geometry-changing/healing/CZ runtime events are rejected; allowed events
  are BodyForceUpdate, ResetDeformationGradient, InitializeStress,
  SetInitialTemperatureAndPressure, Anneal, MaterialSwap, TemperatureProfile,
  and uniform ConfiningPressure.
- **Validation next:** full CPU build, actual MPI migration and split/restart,
  GPU execution, and Silo/VTK readback. Compare cylindrical elastic wave and
  thick-cylinder pressure solutions after traction integration is added.

## Example and checks

`examples/axisymmetricCylinder` creates a radius-1, length-1 cylinder with CPDI
core, CPTI outer annulus and a small radial contraction velocity. It is a smoke
case, not an analytical dynamics benchmark. Its initial volume is pi and its
mass is 1000*pi. Supply ordinary userDefs and run `./runProblem --generate-only`
from a new run directory, or set GEOS_EXECUTABLE and run without that option.
The runner does not submit a scheduler job. `checkParticles.py` verifies the
initial geometry/volume; use the existing plot/history tools for dynamics.

Run from the PFW directory:

```sh
python tests/run_axisymmetric_cpp_tests.py
python tests/run_axisymmetric_writer_tests.py --output /tmp/axisymmetric-pfw
```

Local results: 1,800 production mapping cases compiled with strict warnings,
address/undefined-behavior sanitizers; cylindrical partition/gradient moments,
hoop strain and stress work, uniform-stress radial equilibrium, near-axis
restriction, geometry immutability and actual quad/triangle split plans passed.
The existing production CPDI tests passed 3,000 2D/3D cases. Leak detection
was disabled locally because this container prevents LeakSanitizer process
inspection; address and undefined-behavior checks remained enabled.

Six complete axisymmetric writer runs (ordinary cylinder and curved-surface
refinement, each at 1/3/7 emulated ranks) produced byte-identical files across
rank counts. Four native-2D/3D writer regressions also passed. These are emulated
collectives, not real MPI tests. No full GEOS build, dynamics run, restart or
GPU run was performed here. The target build and smoke run remain acceptance
gates before treating this stage as production validated.

## Uniform-confinement test follow-up

`examples/axisymmetricJacketedSample` models the illustrated sample, steel end
plates and a closed rubber jacket in a single bonded velocity field. It adds
`uniformPressure=1` to ConfiningPressure (no box/sphere exclusion), retains the
pressure state across restart, and corrects axisymmetric box-history density
normalization to use ring volume. The same setup with equal material properties
provides a homogeneous hydrostatic reference. See its README for dimensions,
model assumptions, commands and local validation limits.
