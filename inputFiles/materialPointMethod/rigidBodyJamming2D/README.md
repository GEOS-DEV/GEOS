# Rigid-body MPM MPI/periodic jamming verification

This case exercises the deterministic MPI rigid-body implementation and then
switches to ordinary continuum MPM.

## Discretization and bodies

- plane strain;
- `2 x 3 x 1` MPI decomposition;
- periodic x direction;
- packing between opposing moving `Contact` walls in y;
- 60 physical x cells and 28 physical y cells, giving `DY = 1.5*DX`;
- two particles per x cell and three per y cell, giving square particles;
- `SinglePointBSpline` particles with CPDI scaling disabled;
- six rigid bodies (three disks and three blocks), all in catch-all
  `ContactGroup=0` and distinguished only by `ParticleColor`;
- one disk crosses the x-periodic seam, so its particles appear near both
  `x=-0.5` and `x=+0.5` while rank zero solves it as one unwrapped body.

The rigid event reserves only two nodal fields. The deterministic first color at
a node uses field 0; all remaining colors use the weighted overflow field. Thus
the grid allocation is independent of the six colors and matches a later
single-contact-group, two-DFG-field continuum solve.

## Loading sequence

1. A y-only F-table moves the top and bottom contact walls inward while
   `RigidBodyMPM` is active.
2. The rigid event ends at its timeout, when opposed y-face reactions reach
   `maxJammingStress`, or after the normalized contact penetration limit is
   exceeded on two consecutive steps.
3. A dependent `DeformationUpdate` starts from the actual rigid-event completion
   time and ramps `Syy` to `-140 MPa`, twice the `70 MPa` copper yield stress.
   x remains periodic and is not stress controlled.

## Run

```bash
python3 generateRigidBodyJamming2DParticles.py
srun -n 6 geosx -i mpm_rigidBodyJamming2D.xml -x 2 -y 3 -z 1
python3 checkRigidBodyJamming2D.py --require-run-output
```

`run_rigidBodyJamming2D.sh` performs the same sequence and accepts
`GEOSX_EXECUTABLE` and `MPI_LAUNCHER` environment overrides.
