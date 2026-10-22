# Axisymmetric MPM: implementation and acceptance plan

The implementation includes native axisymmetric mechanics, uniform confinement,
ordinary contact, and prescribed two-group cohesive interfaces. ParticleColor remains
available. The existing timestep dispatches 2D kernels; there is no second
timestep solver. Finite-thickness polymer CZ and damage-field partitioning remain separate follow-ups.

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

## Implemented mechanics

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

## Current scope and next patches

- **Contact:** ordinary contact between nonnegative ContactGroup indices is
  supported with `enableContact=1`. The common P2G normal uses
  `V*(G_r+T, G_y, 0)` to remove the radial interior bias of the ring-volume
  gradient. Existing mass-based normal/friction forces and FMPM impulses need
  no additional circumference multiplier. Explicit PFW surface positions or
  the existing approximate centre-of-mass gap criterion can be used.
  Solver-reconstructed surface positions and overlap correction remain disabled.
  ParticleColor is independent of contact groups.
- **CZ:** prescribed interfaces between groups 0 and 1 support Uncoupled,
  Coupled and Bicrystal traction-separation laws. Ring-weighted line integrals
  and nodal support volumes initialize area; the full relative-F cofactor
  transports it, including hoop stretch and reset invariance. See
  `AXISYMMETRIC_CZ.md` and `examples/axisymmetricCZ` for input, tests and limits.
  Finite-thickness polymer CZ, DFG, multiple field pairs and active-CZ
  subdivision remain future work.
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
  Most geometry-changing/healing runtime events are rejected; allowed events
  are BodyForceUpdate, ResetDeformationGradient, InitializeStress,
  SetInitialTemperatureAndPressure, Anneal, MaterialSwap, TemperatureProfile,
  ReferenceCohesiveZones with explicit normals/positions, and uniform ConfiningPressure.
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

## Ordinary-contact follow-up

`examples/axisymmetricContact` contains two 10-mm-radius, 8-mm-long elastic
cylinders approaching at +/-1 m/s across a 2-mm gap. Each body mixes CPDI core
particles with CPTI skin particles. Separate groups permit separation/rebound;
the within-body particle-type seam stays in one field. Explicit interface
positions gate contact; the corrected volume-gradient estimate supplies normals.
The example defaults to FLIP; `AXIS_CONTACT_UPDATE=FMPM` selects FMPM(2).

The source change is shared by the ordinary and FMPM contact paths. It does
not change the existing normal-selection or friction laws. At a straight
interface, weighting identical unit normals by circumference cannot change
their direction; the new T term instead corrects the cylindrical boundary
identity used to reconstruct normals from bulk particle data.

Local checks compile actual production mapping, normal assembly and both
pairwise laws with strict warnings and address/undefined-behavior sanitizers.
Mixed CPDI/CPTI tiles test axial, cylindrical and inclined interfaces, interior
normal cancellation and axis behavior. 540 force/impulse comparisons cover
radius/mass-ratio/friction variations, supplemented by separation, gap gating,
equal-and-opposite impulses and contact kinetic-energy checks. Four complete
writer runs (FLIP/FMPM, 1/3 emulated ranks) check contact groups, mixed geometry,
ring measures, XML options and byte-identical particle emission.

Run from the PFW directory:

```sh
python3 tests/run_axisymmetric_contact_cpp_tests.py
python3 tests/run_axisymmetric_contact_writer_tests.py --output /tmp/axis-contact-writer
```

The user reported a successful target build and jacketed-sample run before
this contact follow-up. The contact patch still needs a full target build,
GEOS dynamics, actual MPI migration/restart and GPU validation. Its example
postprocessor produces diagnostics, not an established convergence benchmark.
The patch updates manual LaTeX; it does not include a regenerated manual PDF.

Uniform confinement is not selective to the outside of a sealed jacket when
internal interfaces separate. Such an assembly needs exterior-only pressure
loading; the existing bonded jacketed example is left unchanged.

## Prescribed-CZ follow-up

The current patch adds full-ring cohesive areas, corrects reference-relative
area transport and partial-node debonding, and preserves CZ surface data under
quadrature scaling and reset. These shared reference/damage corrections also
apply to Cartesian CZ runs. The user reported that ordinary axisymmetric
contact built and ran before this patch. The new CZ examples still require a
target build and dynamics/MPI/restart checks; see `AXISYMMETRIC_CZ.md`.
