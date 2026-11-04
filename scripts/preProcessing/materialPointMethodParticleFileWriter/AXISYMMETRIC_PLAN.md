# Axisymmetric MPM: implementation and acceptance plan

The implementation includes native axisymmetric mechanics, spatial pressure regions,
ordinary multi-group contact, DFG contact-surface insertion, and prescribed
two-group cohesive interfaces. ParticleColor remains
available. The existing timestep dispatches 2D kernels; there is no second
timestep solver. Finite-thickness polymer CZ and general multi-field CZ remain
separate work.

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
2. Radius-weighted Q1-point, cubic B-spline, CPDI and CPTI operators in the common mapping
   routine, covering cached and on-the-fly mappings. Gradient entry 2 stores
   T=<N/r>, not a z derivative. Coalescing retains it.
3. Radial hoop force and hoop strain rate in FLIP, PIC, XPIC and FMPM; hoop
   strain enters F, physical volume and stress power. Constitutive kernels
   already accept complete 3x3 tensors. Use constitutive symmetries/material
   orientations compatible with no swirl (no r-theta or axial-theta shear).
4. Symmetry at physical r=0. Radial periodicity is rejected. A CPDI/CPTI
   integration domain is uniformly shrunk about its center when necessary to
   avoid negative r; G, G0, mass and volume remain unchanged. A center reaching
   r<=0 is an error, not a clamped trajectory. The same temporary fit is used
   by out-of-range cleanup, so valid near-axis particles are not deleted based
   on their un-fitted material corners. Grid endpoints are accepted within
   floating-point tolerance; genuine outflow escapes still receive delete flags.
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

## Axis boundary and surface-detection corrections

`particleType=1` is accepted for axisymmetric bulk/ordinary-contact particles.
At the physical axis, the cubic mapper replaces the negative-radius control
value by `u[-1]=2*u[0]-u[1]` before particle/grid transfers. For `s=r/h` in the
first radial cell the physical weights are `(1-s+s^3/6, s-s^3/3, s^3/6)` at
`r=(0,h,2h)`. They are nonnegative, linearly complete and C2 at the transition
to the interior cardinal basis. Both mapping paths use the same derivatives
and `T=N/r` hoop operator. Direct polynomial evaluation retains tiny positive
radii without a numerical radius floor. No additional particle state is stored.

This closure acts only at the physical revolution axis; Cartesian/F-table,
annular, other physical grid faces and MPI partition faces keep their original
basis. It differs from post-mapping continuation of constrained velocity
components: mass and force assembly use the new basis, too. The one-sided
second radial derivative is zero at the axis; zero axial-velocity slope is
not strongly imposed. Check curved-profile accuracy by refinement. The grid-driving follow-up enables moving/contact/F-table walls using this
existing basis and constrained-velocity continuation on the other faces.
B-spline CZ-tagged particles remain rejected until the Q1 cohesive
area/reference-support operator is generalized. Untagged bulk can be B-spline.

`surfaceDetection=1` now measures the meridional area-weighted centroid,
including axis images within the unchanged neighbor support. It keeps the
existing threshold and surface-flag convention. Uniform material no longer
has a radius-induced centroid shift and the axis is not a free surface.
True axial/cylindrical boundaries remain detectable. This is separate from
reconstruction of contact surface positions/normals, now supported as described in
`AXISYMMETRIC_SURFACES.md`.
Previously assigned flags are retained, including on restart.

`tests/run_axisymmetric_boundary_cpp_tests.py` compiles production mapping,
surface-neighbor, cleanup and event-handler excerpts with strict warnings and
address/undefined-behavior sanitizers. It checks the old false-surface case,
axis/free-surface discrimination, valid near-axis cleanup, true escapes,
B-spline moments/derivatives/hoop work down to r/h=1e-20, hydrostatic balance,
consistent overlapping local stencils and the friction-swap allow-list.
These tests do not run full GEOS, actual MPI/restart or GPU execution.

## Current scope and next patches

- **Contact:** ordinary contact between nonnegative ContactGroup indices is
  supported with `enableContact=1`. The common P2G normal uses
  `V*(G_r+T, G_y, 0)` to remove the radial interior bias of the ring-volume
  gradient. Existing mass-based normal/friction forces and FMPM impulses need
  no additional circumference multiplier. Explicit PFW surface positions or
  the existing approximate centre-of-mass gap criterion can be used.
  Solver-reconstructed surface positions are supported; overlap correction remains disabled.
  ParticleColor is independent of contact groups.
- **DFG:** `damageFieldPartitioning=1` is accepted with ordinary contact.
  Reconstruct the damage gradient using meridional area `V/(2*pi*r)` and an
  even extension across the axis. Existing neighbor selection, field mapping,
  separability criteria and their defaults are retained. Physical grid
  measures keep full-ring weights. See `AXISYMMETRIC_DFG.md` and
  `examples/axisymmetricDFG`. This enables contact-surface insertion, not new CZ bonds.
- **CZ:** prescribed interfaces between groups 0 and 1 support Uncoupled,
  Coupled and Bicrystal traction-separation laws. Ring-weighted line integrals
  and nodal support volumes initialize area; the full relative-F cofactor
  transports it, including hoop stretch and reset invariance. See
  `AXISYMMETRIC_CZ.md` and `examples/axisymmetricCZ` for input, tests and limits.
  Combining CZ with DFG, general multi-field CZ, finite-thickness polymer CZ
  and active-CZ subdivision remain separate work.
- **Pressure and tractions:** BoreholePressure selects r < r0 independently
  of axial y. ConfiningPressure supports a pressure-free (r,y) box, an
  origin-centered sphere exclusion r*r+y*y <= R*R, or uniformPressure=1.
  All use background stress including hoop pressure. Surface reconstruction,
  nodal areas and diffuse surface tension are supported; see `AXISYMMETRIC_SURFACES.md`.
  Generic PFW surface-traction generation remains rejected. Direct GEOS
  particle tractions, when supplied, are already integrated full-ring forces;
  no automatic per-area or per-radian conversion is performed.
- **Weibull:** built-in area-only native Weibull objects are rejected.
  Statistical volumes require 2*pi*integral(r dA), evaluated over each actual
  statistical cell/layer before minimum-volume clamping. Custom wrappers may
  implement setAxisymmetric(); explicit StrengthScale values are unchanged.
  Weibull and localized hoop-driven fracture are deferred: circumferential
  localization breaks the assumed axisymmetry and needs a separate treatment.
- **Other controls:** no swirl, SPH/overlap correction,
  crack-tip detection, profileHistory slice-volume normalization, or periodic
  affine-bulk `prescribedFTable` driving. Boundary F tables, radial/axial stress
  control, BC tables and Outflow/Symmetry/Moving/Contact walls are supported.
  Contact walls remain unavailable with FMPM. Existing point/CPDI/CPTI body forces,
  initial velocity and compatible stress/thermal initialization remain usable.
  Most geometry-changing/healing runtime events are rejected; allowed events
  are BodyForceUpdate, ResetDeformationGradient, InitializeStress,
  SetInitialTemperatureAndPressure, Anneal, MaterialSwap, TemperatureProfile,
  ReferenceCohesiveZones with explicit or reconstructed normals/positions, ConfiningPressure, BoreholePressure,
  and FrictionCoefficientSwap (its existing scalar/table update and completion).
- **Validation next:** full CPU build, actual MPI migration and split/restart,
  GPU execution, and Silo/VTK readback. Compare cylindrical elastic wave and
  thick-cylinder pressure solutions using the pressure-region example below.

## Example and checks

`examples/axisymmetricCylinder` creates a radius-1, length-1 cylinder with CPDI
core, CPTI outer annulus and a small radial contraction velocity. It is a smoke
case, not an analytical dynamics benchmark. Its initial volume is pi and its
mass is 1000*pi. Supply ordinary userDefs and run `./runProblem --generate-only`
from a new run directory, or set GEOS_EXECUTABLE and run without that option.
Set `AXIS_CYLINDER_INTERPOLATION=bspline` for the B-spline variant, and
`AXIS_CYLINDER_UPDATE=FMPM` for FMPM(2). The runner does not submit a scheduler job. `checkParticles.py` verifies the
initial geometry/volume; use the existing plot/history tools for dynamics.

Run from the PFW directory:

```sh
python tests/run_axisymmetric_cpp_tests.py
python tests/run_axisymmetric_boundary_cpp_tests.py
python tests/run_axisymmetric_writer_tests.py --output /tmp/axisymmetric-pfw
```

Local results: 2,400 production mapping cases compiled with strict warnings,
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

## Spatial pressure regions

`examples/axisymmetricPressureRegions` supplies a hollow cylinder with
independent inner/outer pressure ramps and a sphere-confinement variant.
Pressure masks are evaluated on current grid coordinates with fixed physical
cutoffs. They do not follow moving material surfaces. No surface flags or
particle traction records are required.

For material a < r < b choose a < r0 < r1 < b, with cutoff margins appropriate
to the interpolation support and expected deformation. BoreholePressure uses
`boreholeRadius=r0`. The outer ConfiningPressure box has maximum radial entry
r1, minimum radial entry below every grid node, and y bounds beyond every grid
node (including buffers). Box XML still has three entries; z is ignored.
The pressure-free box interior then leaves only r > r1 loaded. A sphere
exclusion uses r*r+y*y > R*R and revolves about y. Equalities are not selected.
Pressure values are positive in compression; all three background normal
stresses, including hoop, are -p. Full-ring particle volumes already provide
the circumference factor.

Confining pressure takes precedence if masks overlap. The solver retains one
borehole and one confining region; successive events of the same kind replace
the shared state. Axisymmetric pressure handlers set the exact endPressure on
the step reaching the ramp end, then hold it; a zero-duration event sets its
endPressure. The general event scheduler and Cartesian time sampling are
unchanged. Both pressure states are registered for restart, so completed
events need not run again to restore a held load.

Pressure acts on any boundary in selected nodal support, not on a named
external surface. Radial bands can also load end-face bands or internal
surfaces. Check those intersections when using contact, cracks or a deforming
jacket. Generic PFW per-area traction integration remains separate work.

`tests/run_axisymmetric_pressure_cpp_tests.py` compiles production selectors
and force assembly and checks uniform/box/sphere/borehole selection, overlap
precedence, moved-grid reevaluation, mixed CPDI/CPTI inner/outer cylindrical
resultants, hoop balance, pressure work and Cartesian regressions.
`tests/run_axisymmetric_pressure_event_cpp_tests.py` compiles pressure event
validation and handlers, checks ramp endpoints/completed gates, and statically
checks restart registrations. `tests/run_axisymmetric_pressure_regions_writer_tests.py
--output NEW_DIR` generates both examples and the uniform-pressure regression
at 1/3 emulated ranks. These are not a full solver dynamics or restart test.

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

## Axisymmetric DFG/contact follow-up

The focused DFG patch removes axisymmetric rejection in the solver and PFW,
then supplies a meridional normalized-kernel damage gradient with axis parity.
It does not alter the common separability law or defaults, neighbor-selection
policy, ordinary pairwise contact laws, or CZ machinery. Ordinary contact
already loops over all occupied velocity-field pairs; examples now exercise
three groups, optionally with two DFG slots per group.

The new example has `intact`, `damage`, `contact` and `combined` variants,
with FLIP/FMPM selections and CPDI/CPTI particles. Damage is prescribed initially
to isolate the partition/insertion path; no new constitutive damage law is
introduced. See its README for commands and limits. The new CPU and writer
checks are listed in `AXISYMMETRIC_DFG.md`; full target dynamics, MPI and restart
validation are still needed for this follow-up.

## Radial/axial grid-driving follow-up

- `prescribedBoundaryFTable=1`: rows `[time,Frr,Fyy,1]`. Geometry uses exact
  finite ratios; wall velocities retain endpoint Fdot/F. `domainF` and `domainL`
  remain solver state, not new PFW/XML inputs. Use the F table to prescribe them.
- `stressControl=[1,0,0]` targets the full-domain-volume average of sigma_rr,
  including void. `[0,1,0]` targets axial stress; `[1,1,0]` controls both.
  This is not wall pressure or radial-plus-hoop generalized stress. Existing
  legacy PID, adaptive-seek (0/1 only) and one-sided limits (2) are retained.
- Grid zz stays inactive (F=1, L=0); particle hoop rate is still v_r/r.
  Full-ring domain volume scales with Frr squared times Fyy. Both radial faces
  of an annulus share one scale about r=0; axial scaling is about y=0.
- Moving faces impose velocity bilaterally; Contact faces are unilateral,
  with inelastic normal response and impulse-limited meridional Coulomb
  friction. No additional 2*pi*r factor enters the velocity correction.
  Contact without explicit surface positions retains the center-of-volume
  gap estimate. Contact walls require FLIP/PIC/XPIC; Moving also permits FMPM.
- A BC schedule cannot release/move r=0. Driven directions need a Moving or
  Contact face; fixed Symmetry away from the origin and axial periodic driving
  are rejected. Swirl and uniform nonzero radial cap velocities at r=0 are
  rejected. Other deformation/healing event restrictions remain unchanged.
- `examples/axisymmetricGridDriving` includes stretch/stress/mixed modes,
  `runProblem`, and a post-processor which converts material-average stress
  to the whole-domain normalization used by the controller.

Validation: `tests/run_axisymmetric_driving_cpp_tests.py` compiles production
interpolation, grid maps, wall targets, volume/stress metrics and both controller
bodies with strict warnings and ASAN/UBSAN. CPU array, reduction, MPI and material
lookup adapters are used; this is not a full GEOS build.
`tests/run_axisymmetric_driving_writer_tests.py --output NEW_DIR` runs eight
full-writer cases using one/three emulated ranks and validates accepted/rejected
inputs. Real MPI migration, restart and loading convergence still need GEOS runs.
