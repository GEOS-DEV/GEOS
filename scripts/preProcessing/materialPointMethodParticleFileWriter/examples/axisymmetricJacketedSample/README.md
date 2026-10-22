# Axisymmetric jacketed-sample confinement test

A small model of the supplied sketch: a cylindrical specimen between steel end
plates, enclosed by a rubber jacket with end caps, immersed in a virtual
pressure bath. Apply this patch **after GEOS-mpm-axisymmetric-stage1.patch**, on
the branch with rigid-body MPM removed, then rebuild GEOS. This patch adds the
example and the small uniform-pressure extension it needs.

## Model and defaults

All quantities use SI units. Coordinates are x=r, y=axial, geometric z=0. The
full axial height is modeled. Only r=0 is a symmetry boundary; the other grid
boundaries are outflow and lie outside the initially occupied material.

| Part/control | Default |
|---|---|
| Sample | Radius 10 mm, height 20 mm; effective density 800 kg/m3, E=50 MPa, nu=0.20 |
| Steel plates | Two 2 mm discs, radius 10 mm; density 7850 kg/m3, E=200 GPa, nu=0.30 |
| Rubber jacket | 1 mm side wall and two 1 mm caps; density 1100 kg/m3, E=5 MPa, nu=0.45 |
| Outer dimensions | Radius 11 mm, height 26 mm |
| Pressure | Cosine ramp from 0 to 10 kPa over 1 ms, followed by a 1 ms hold |
| Grid | 1 mm cells; 2 particles per cell direction; 1,328 initial particles |
| Particle domains | CPDI quads in sample/plates; CPTI triangles in jacket/caps |
| Interfaces | Bonded through shared velocity field/contact group 0 |
| Time integration | ExplicitDynamic + PIC; CFL factor 0.25 |

These properties are illustrative. The sample is a homogeneous *effective
porous elastic continuum*: there are no resolved pores, pore pressure, Darcy
flow, compaction law, plasticity or damage. Rubber is a nearly incompressible
elastic surrogate appropriate to this small-load test, not a calibrated
finite-strain rubber model. Steel is deformable and its high wave speed limits
the explicit timestep. PIC introduces numerical damping; inspect the histories
before interpreting the result as quasistatic.

The fluid is represented only by prescribed external pressure. There are no
fluid particles, fluid inertia or fluid/solid flow coupling. All interfaces
are bonded: no sliding, opening, friction or CZ. This is the supported simple
starting point before axisymmetric contact is added. No Weibull or
symmetry-breaking hoop fracture is included.

## Run

From the GEOS checkout root, after applying the patch and rebuilding:

```bash
case_dir="$PWD/scripts/preProcessing/materialPointMethodParticleFileWriter/examples/axisymmetricJacketedSample"
export GEOS_EXECUTABLE=/absolute/path/to/rebuilt/geosx
mkdir -p runs/jacketed
cd runs/jacketed
"$case_dir/runProblem" --generate-only
```

The runner copies your `userDefs_<username>.py` from the PFW directory if it is
present. Otherwise place your configured file in the run directory first.
`PFW_PYTHON` can select the Python interpreter with your normal PFW dependencies.
The generated `initialGeometry.png` shows the actual particle domains.

Run the already-generated input inside an allocation as required by your site:

```bash
"$GEOS_EXECUTABLE" -i mpm_axisymmetricJacketedSample.xml > geos.log 2>&1
python3 "$case_dir/postProcess.py"
```

Alternatively, call `"$case_dir/runProblem"` without `--generate-only` from a
fresh run directory to generate, check, run and plot in sequence. The runner
uses one GEOS rank, does not submit a scheduler job, and refuses to overwrite
existing case files. `case_dir` is an absolute path and survives `cd`.

For a homogeneous reference run in a different fresh directory:

```bash
mkdir ../homogeneous
cd ../homogeneous
JACKET_CASE=homogeneous "$case_dir/runProblem"
```

This keeps the same geometry/material labels and mixed CPDI/CPTI domains but
assigns the sample's density/E/nu to every region. It removes stiffness contrast
and tests uniform confinement across the particle-type interfaces.

Edit `caseSettings.py`, or use these environment overrides **during generation**:
`JACKET_CASE=jacketed|homogeneous`, `JACKET_REFINE=2`, `JACKET_PRESSURE=10000`,
`JACKET_RAMP_TIME=0.001`, `JACKET_END_TIME=0.002`. The runner records the resolved
settings in `case.json`; postprocessing reads that file, not the current
environment. Increase end time if the ramp time is increased. Refinement doubles
the resolution in each direction and gives four particles across jacket thickness.

## Loading and source changes

The example uses:

```xml
<ConfiningPressure name="pressureRamp" startTime="0" endTime="0.001"
  startPressure="0" endPressure="10000" interpType="1" uniformPressure="1"/>
```

`uniformPressure=1` is a new optional event attribute (default 0). It sets the
background stress to `[-p,-p,-p,0,0,0]` everywhere, without a box or sphere
exclusion. The pressure enters the existing axisymmetric weak force:

- radial: `-V * [(sigma_rr+p)*G_r + sigma_ry*G_y + (sigma_tt+p)*T]`;
- axial: `-V * [sigma_ry*G_r + (sigma_yy+p)*G_y]`.

Here V is the full-ring particle volume and T is the existing hoop operator.
No extra factor of 2*pi, thickness or estimated surface area is applied. For a
conforming assembly, internal-interface contributions cancel and the outer
boundary receives pressure. Thus the sealed outer jacket/caps are loaded
without applying pressure at each material seam. A pressure field has no
knowledge of cavity connectivity; do not use this test to represent isolated
unpressurized internal cavities by simply adding voids.

The existing Cartesian box/sphere modes retain their behavior. Axisymmetry
still rejects those spatial exclusions, borehole pressure, nonuniform surface
tractions, contact/CZ and other unsupported stage-1 features. PFW checks the
pressure event before writing files. Pressure enable/mode/stress/exclusion
state now persists on restart, including after a ramp has completed. Restart
I/O still needs verification on the full GEOS build.

A related fix normalizes axisymmetric box-history density by the selected ring
volume. `Sxx/Syy/Szz` are volume-weighted radial/axial/hoop material stresses.
The box history's `F00/F11/F22` describe the **background domain**, so they are
not used as specimen strains here. The sample tracer's density gives volume
strain `rho_initial/rho - 1` under mass conservation.

## Outputs and expected behavior

- `initialChecks.json`: per-material counts, full-ring volumes and masses.
- `initialGeometry.png`: actual quad/triangle domains with pressure directions.
- `boxAverageHistory.csv`: the sample interior (r<9 mm, |y|<9 mm).
- `jacketTracer_000000.csv` through `000004.csv`: sample center, sample near upper
  plate, upper steel plate, side jacket, upper jacket cap.
- `histories.png`: sample stresses, tracer axial stresses, sample volume strain,
  and radial displacement relative to each tracer's first record.
- `historyChecks.json`: final diagnostics and homogeneous-reference checks.
- Standard Silo and restart output from GEOS. Glyphs remain an r-y section;
  stressZZ is hoop stress, not out-of-plane plane-strain stress.

Expect inward radial motion and axial compression with top/bottom symmetry.
Bonding and stiffness contrast produce end effects in the heterogeneous case;
its sample stress need not be exactly hydrostatic. There is no analytical
pass/fail check for that case.

For the homogeneous variant, after settling, the small-strain reference is
`sigma_rr = sigma_yy = sigma_tt = -p` and `dV/V = -p/K`, where
`K=E/[3(1-2*nu)]`. At the defaults, K=27.7778 MPa, each normal strain is
approximately -0.00012 and volume strain approximately -0.00036 (-0.036%).
The postprocessor checks the final 20% of the hold: each mean stress within 5%
of p, stress oscillation range below 10% of p, and volume strain within 20% of
-p/K. These are initial smoke-test tolerances, not established convergence
bounds. A failure needs inspection, a longer ramp/hold and a resolution study.

## Local verification and limits

The full writer generated both variants on 1 and 3 **emulated** ranks, with
byte-identical particle files across counts. Every initial domain had positive
area/radius, the expected material placement, and `V=2*pi*r*A`; per-material
volumes matched the cylinder/annulus formulas. The default total volume is
9.883450488e-6 m3; the heterogeneous total mass is 0.01746914011 kg.

A standalone C++ test compiles the production mapping, background-pressure and
2D force bodies with strict warnings, ASan and UBSan. It verifies pressure
resultants, cancellation at bonded CPDI/CPTI interfaces and the axis, hydrostatic
balance and pressure work, plus unchanged Cartesian box/sphere selection.

From the PFW directory:

```bash
python3 tests/run_axisymmetric_pressure_cpp_tests.py
python3 tests/run_jacketed_sample_writer_tests.py --output /tmp/new-jacketed-checks
```

Local LeakSanitizer process inspection is unavailable in the container, so the
C++ test was run with `ASAN_OPTIONS=detect_leaks=0`; ASan/UBSan remained enabled.
Postprocessing was exercised with synthetic histories only. No full GEOS build,
dynamics, actual MPI, restart or GPU run was performed here. Those remain
acceptance checks on the target machine. The manual's LaTeX sources are updated;
this small example patch does not regenerate the compiled manual PDF.
