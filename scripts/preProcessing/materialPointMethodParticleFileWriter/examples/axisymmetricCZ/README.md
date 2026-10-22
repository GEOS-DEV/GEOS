# Axisymmetric cohesive-interface smoke test

Two touching elastic cylinders, each radius 10 mm and length 8 mm, initially
move apart at 0.1 m/s. The interface at y=0 is bonded by an Uncoupled CZ with
normal/shear stiffness 1e8 Pa/m. Density is 1000 kg/m^3, E=10 MPa, nu=0.
Each body has a CPDI core and CPTI outer layer, with 768 particles in total.
The two groups are 0/1 and CZTag is 1. The axis uses radial symmetry; the
remaining domain boundaries are outflow, sufficiently far from the sample.

The input explicitly sets `useEvents=1`. GEOS defaults this switch to zero:
without it, `ReferenceCohesiveZones` is parsed but never runs, so the initially
separating cylinders remain unbonded. `SurfaceFlag=3` and `CZTag=1` identify
interface candidates; the runtime `particleCohesiveZoneFlag` is set by CZ
activation. The generated file has 48 interface candidates per body. PFW and
the example checker now reject an omitted or disabled event switch.

If correcting an input generated before this fix, set `useEvents="1"` on
`SolidMechanics_MPM` in the generated XML and start again at time zero. Do not
resume an unbonded restart past the event's activation time. This input fix
alone needs no further C++ rebuild once the axisymmetric CZ patch is installed.

The subsequent reference-event source fix does require rebuilding GEOS. Each
`activateBond` event must initialize the reference only once, even when the
timestep grows from its startup value. Previously, a 1 ns initial step
followed by a 2.4975 microsecond step triggered activation at both cycles 0
and 1. In the failing run, the logged area sum fell from 6.28319e-4 m^2 to
3.14159e-11 m^2 on the second initialization. This overwrote the bonded
reference; it was not failure at the specified opening limit. The new
internal activation flag is saved across restart. Rerun from time zero
after this fix; restarting the damaged run will not restore its reference.

The default `elastic` variant sets a very large failure displacement, giving
a compliant intact interface with a 1 m failure displacement. `debond` uses a 5 micrometre normal opening
limit. The latter should release under the tensile pulse. This is a target-run
smoke test, not a previously validated quantitative dynamics benchmark.

## Generate and run

After applying the patch, rebuild GEOS. Run from a **fresh run directory**:

```bash
example=/usr/WS1/homel1/geos/scripts/preProcessing/materialPointMethodParticleFileWriter/examples/axisymmetricCZ
"$example/runProblem" --generate-only
```

Supply the usual `userDefs_<username>.py` with the rebuilt executable's
`geosPath`, either in the run directory or the PFW root. The runner invokes
PFW and checks the resulting particles/XML. Submit the generated GEOS batch
script through your site's scheduler; the runner does not submit it.

After that run completes:

```bash
python3 "$example/postProcess.py"
```

Alternatively, on an appropriate compute allocation, set `GEOS_EXECUTABLE`
and run without `--generate-only` for a direct single-rank run and plots.
`GEOS_EXECUTABLE` affects this direct launch only; the generated batch script
uses `geosPath` from userDefs.

Use separate fresh directories for variants:

```bash
AXIS_CZ_CASE=debond "$example/runProblem" --generate-only
AXIS_CZ_UPDATE=FMPM "$example/runProblem" --generate-only
AXIS_CZ_DISPLACEMENT=TypeB "$example/runProblem" --generate-only
```

Accepted choices: elastic/debond; FLIP/FMPM; Nodal/TypeA/TypeB. Controls can be
combined. The job ends at 2 ms and requests a restart at 1 ms. Resume using your
normal GEOS restart invocation and compare with an uninterrupted run; restart
execution is not part of the sandbox validation.

## Inspect

- `initialChecks.json`: finite 2D geometry, both particle types in each body,
  full-ring mass/volume, equal/opposite initial momentum, surface flags/tags,
  opposing axial normals, interface positions at y=0, `useEvents=1` and
  canonical XML enums. Candidate counts are input checks, not confirmation
  that the solver has activated the bond.
- `czHistories.png`: four tracer histories of axial velocity, tensile stress,
  particle-centre position and radial velocity. `czDiagnostics.json` reports
  their final velocities and peak axial tension. The script checks finite
  histories, completed time coverage and zero theta velocity, but does not
  assert an unverified dynamics accuracy tolerance.
- Grid `gridCohesiveForce`, `gridCohesiveArea`, `gridDisplacement`,
  `gridCohesiveNode`, `gridContactForce`, and particle cohesive/surface flags:
  inspect intact load transfer and release. Grid CZ data is attached to the
  cached reference grid, not to an explicitly advected interface mesh.

Plot `particleStressYY` for axial tension. `particleStressZZ` is hoop stress
in axisymmetry; this example uses `nu=0`, so zero hoop stress alone does not
indicate a missing bond. Look for the initialization log's nonzero
`summed reference areas (both sides)` and active cohesive flags after the
first timestep, followed by CZ forces and axial stress as separation develops.
For this single-event example, there should be only one
`Enabling cohesive zone region, bond` message and one reference initialization
in a fresh run. A second initialization at cycle 1 indicates that the
reference-event source fix is not in the executable being run.

The physical interface area is pi*(0.01 m)^2. The logged sum includes **both
sides**, so compare half that sum with pi*R^2 and refine the outer edge before
setting an accuracy tolerance. The local-line/volume-normalization estimate
can differ from the finite disk's exact area at its rim.

Tracer-centre separation includes deformation of the neighbouring particles;
it is not the cohesive opening. For displacement diagnostics use the difference
of the two cached `gridDisplacement` fields projected on the interface normal.
With `preventCZInterpenetration=1`, ordinary contact supplies compression and
CZ supplies tensile/shear forces; after full debonding only ordinary contact
remains. Partial failure retains intact cached-node forces.

The PFW documentation `AXISYMMETRIC_CZ.md` specifies supported models,
reference-area transport, reset behavior, restrictions and target tests.
Finite-thickness polymer CZ models and DFG are not enabled by this example.
