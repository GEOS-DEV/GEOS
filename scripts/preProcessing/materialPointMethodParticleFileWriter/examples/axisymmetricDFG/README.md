# Axisymmetric contact and DFG isolation cases

A true r-y mesh with CPDI bulk and CPTI outer skin. All cylinders have radius
6 mm. Material is ElasticIsotropic, density 1000 kg/m^3, E=10 MPa, Poisson ratio
zero. Grid spacing is 1 mm, with two particles per cell direction before CPTI
tessellation. The axis is radial symmetry; other boundaries are remote outflow.
There are no cohesive regions/events or explicit particle surface normals.

| AXIS_DFG_CASE | Groups / allocated fields | Initial condition | Purpose |
| --- | --- | --- | --- |
| intact | 1 / 2 | Continuous 12-mm cylinder, zero damage; halves move axially at +/-0.1 m/s | Intact comparison |
| damage (default) | 1 / 2 | Same cylinder with d=1 in a 1-mm central band | DFG opening from a prescribed damage field |
| contact | 3 / 3 | Outer cylinders approach at +/-0.5 m/s; central body initially stationary | Ordinary contact with more than two groups |
| combined | 3 / 6 | Three bodies, central damage band and initial central velocities +/-0.1 m/s | DFG and ordinary multi-group contact together |

The three-body setup has axial intervals [-8,-3], [-2,2], [3,8] mm for groups
0,1,2. It initially has 1-mm gaps between bodies. Run time is 1 ms for the
single-body cases and 4 ms for the three-body cases. Initial velocities are
assigned once; no tensile traction or maintained end velocity is applied.

Damage is an initial field, not a new failure model. ElasticIsotropic does not
generate additional damage. The existing DFG law decides opening/contact;
this example does not change that law or its defaults. Its input explicitly
sets minimum mapped damage and quality thresholds to 0.5, disables single-field
fraction suppression and fully-damaged F resets, and uses neighborRadius=1.5 mm.
These settings isolate DFG partitioning for this smoke case. They are not a
material calibration. Surface reconstruction and overlap correction remain off.

## Run

Use your usual PFW userDefs configuration. From a fresh run directory, with
`GEOS_ROOT` pointing at your source checkout:

```sh
export AXIS_DFG_CASE=damage        # intact, damage, contact, combined
export AXIS_DFG_UPDATE=FLIP        # FLIP or FMPM (order 2)
"$GEOS_ROOT/scripts/preProcessing/materialPointMethodParticleFileWriter/examples/axisymmetricDFG/runProblem" --generate-only
```

Generation writes XML, particles, case.json and initialChecks.json. Submit the
XML with your normal scheduler/run command. Or set `GEOS_EXECUTABLE` to the
rebuilt executable and omit `--generate-only` for a direct single-process run.
Use a separate fresh directory for each case/update combination.

After a scheduler run:

```sh
python3 "$GEOS_ROOT/scripts/preProcessing/materialPointMethodParticleFileWriter/examples/axisymmetricDFG/postProcess.py"
```

The postprocessor plots axial/radial velocity, axial stress (tension positive)
and tracer positions. It checks finite complete histories and no hoop velocity;
it does not assert an analytical dynamics tolerance. Compare intact/damage
near-midplane motion and the contact/combined response to the approaching
outer bodies. Grid fields include damage, damage gradient and contact force;
particle output retains damage and contact-group data.

## Verification boundary

The included particle checker tests groups, damage, ring volume, geometry,
zero explicit surface flags/normals, initial momentum and solver XML. The
standalone C++ adapter can map the generated particles to verify central-band
separation under the existing criterion (commands in ../../AXISYMMETRIC_DFG.md).
This is still a smoke example requiring a target GEOS run, not an established
fracture benchmark. Multi-field CZ and separability-law development are separate.
