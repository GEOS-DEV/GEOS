# rigidBodyJamming2D

PFW-first MPI/periodic rigid-body MPM jamming example. It creates six initially
separated 2D copper bodies, assigns one `ParticleColor` per body, and keeps every
body in the same catch-all `ContactGroup=0`. The rigid stage uses only two nodal
fields, so the number of global colors does not increase the grid-field storage.
The last field is the weighted overflow field for the rare node touched by more
than two colors.

The verification geometry and discretization exercise the MPI/periodic path:

- plane strain;
- `2 x 3 x 1` MPI partitions;
- periodic x direction;
- opposing moving contact walls and packing in y;
- one disk split across the x-periodic seam but represented by one color/body;
- `SinglePointBSpline` particles with CPDI scaling disabled;
- rectangular initial cells with `DY = 1.5*DX`;
- PFW `ppcx=2`, `ppcy=3`, which gives equal x/y particle spacing and square
  particle domains.

After the rigid event reaches `maxJammingStress` or another rigid safeguard, a
dependent `DeformationUpdate` starts from the actual packing-completion time and
ramps the y stress to `2 x` the copper yield strength.

Run in the examples workflow:

```bash
cd ~/geos/mpm/examples/rigidBodyJamming2D
./runProblem -y
```

Or copy `pfw_input_rigidBodyJamming2D.py` and adapt it for the normal PFW
`runClean_<user>.sh rigidBodyJamming2D` workflow.

Generated suite products include:

- `rigidBodyJamming2D_color_four_times.png`: `ParticleColor` snapshots at four
  output states;
- `rigidBodyJamming2D_reactions_y.png`: y-direction reaction history;
- `rigidBodyJamming2D_reactions_y.csv`: plotted reaction data.

After PFW creates the particle file, run the geometry/discretization checks with:

```bash
python3 checkRigidBodyJamming2D.py mpmParticleFile_rigidBodyJamming2D
```
