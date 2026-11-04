# Axisymmetric grid driving

A unit-radius, unit-height elastic cylinder in radial/axial grid walls. The
mesh and particles are 2D; mass, volume and reactions represent full rings.
This example exercises the new controls; a full GEOS dynamics run is still
required to assess controller tuning and mesh/time-step convergence.

From a new run directory with the usual PFW `userDefs` configuration:

```bash
export GEOS_EXECUTABLE=/path/to/geosx
/path/to/examples/axisymmetricGridDriving/runProblem
AXIS_GRID_MODE=stress /path/to/examples/axisymmetricGridDriving/runProblem
AXIS_GRID_MODE=mixed /path/to/examples/axisymmetricGridDriving/runProblem
```

Use a separate run directory for each mode. `--generate-only` stops after PFW.
The default `stretch` mode ramps Frr from 1 to .99 while Fyy=1. `stress` ramps
mean sigma_rr from zero to -10 kPa with Fyy=1. `mixed` applies the same radial
stress target while prescribing Fyy from 1 to .995. All ramps finish at .2 s,
then hold to .3 s. The normal stress sign is positive in tension.

`AXIS_GRID_WALL=contact` makes only the outer radial wall unilateral. Its
friction coefficient defaults to zero; set `boundaryFaceFrictionCoefficients`
to add Coulomb friction. Without explicit surface data the contact gate uses
the mapped center of volume, so surface/gap resolution should be checked.
`AXIS_GRID_UPDATE=FMPM` is supported with Moving walls; contact walls require
FLIP/PIC/XPIC. The axis and bottom platen remain Symmetry.

The outer wall is a cylinder. Its reaction is the sum of radial components
around the ring, not a Cartesian vector resultant. Its current area is
2*pi*R*H. Axial platen area is pi*R^2. The stress controller uses the integral
of sigma_rr over material divided by the whole cylindrical domain volume,
not reaction/area. These measures need not agree in a nonuniform dynamic test.
`postProcess.py` reports the controller-compatible average from the history
files and accounts for their different normalizations.

For a small homogeneous strain with Fyy=1 and E=1 MPa, nu=.25, linear elasticity
gives sigma_rr=sigma_hoop approximately 1.6 MPa*(Frr-1). Thus the default .99
radial stretch gives about -16 kPa; a -10 kPa stress target corresponds to
Frr about .99375. These are small-strain estimates, not asserted run results.
