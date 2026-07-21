# Rigid-body MPM jamming verification

This verification case exercises the color-partitioned rigid-body MPM event and
then switches back to the continuum MPM update.

Run sequence:

```bash
python3 generateRigidBodyJamming2DParticles.py
geosx -i mpm_rigidBodyJamming2D.xml
python3 checkRigidBodyJamming2D.py
```

The geometry contains six initially non-overlapping 2D objects: three disks and
three blocks.  All objects use one copper material model and one contact group,
but each object has a distinct `ParticleColor`.  The rigid event therefore must
form one rigid body per color and must assign local nodal velocity fields by
color, not by contact group.

The particle generator writes exact analytic surface normals and surface-position
vectors for particles flagged as analytic surface particles.  Interior particles
carry zero normals and positions so they do not contribute spurious contact
geometry.

The loading sequence is:

1. F-table driven boundary compaction while `RigidBodyMPM` is active.
2. Completion at the rigid event end time or when the maximum resultant body
   force reaches the configured jammed-state threshold.
3. A dependent `DeformationUpdate` switches the solver to continuum stress
   control and ramps `sigma_x = sigma_y = -140 MPa`, which is twice the 70 MPa
   copper yield stress used in the XML.

The checker always validates the generated geometry and surface data.  After a
GEOS run it also checks for `rigidBodyHistory.csv` and box-average stress output
when those files are present.
