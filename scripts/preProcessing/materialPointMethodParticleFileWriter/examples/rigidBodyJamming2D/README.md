# rigidBodyJamming2D

PFW-first rigid-body MPM jamming example.  It creates six initially separated
2D copper objects, assigns one `ParticleColor` per object while keeping a single
`ContactGroup`, runs a `RigidBodyMPM` compaction event, then switches to
continuum MPM and drives y-compression to `2 x` the copper yield strength.

Run in the examples workflow:

```bash
cd ~/geos/mpm/examples/rigidBodyJamming2D
./runProblem -y
```

Or copy `pfw_input_rigidBodyJamming2D.py` and adapt it for the normal PFW
`runClean_<user>.sh rigidBodyJamming2D` workflow.

Generated suite products include:

- `rigidBodyJamming2D_color_four_times.png`: ParticleColor snapshots at initial,
  quarter, three-quarter, and final output states.
- `rigidBodyJamming2D_reactions_y.png`: compression-direction reaction history.
- `rigidBodyJamming2D_reactions_y.csv`: plotted reaction data.

After PFW creates the particle file, the geometry and particle-column checks can
be run with:

```bash
python3 checkRigidBodyJamming2D.py mpmParticleFile_rigidBodyJamming2D
```
