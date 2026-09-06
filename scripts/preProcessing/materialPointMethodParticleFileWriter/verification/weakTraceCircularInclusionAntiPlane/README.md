# weakTraceCircularInclusionAntiPlane

Circular inclusion benchmark for explicit-surface weak trace development under
anti-plane shear.

The folder dispatches three variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact groups tied by the weak trace projection.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.

The top and bottom grid boundaries receive opposite out-of-plane transverse
velocities.  This case uses a thin 3D slab, rather than GEOS plane-strain mode,
so the out-of-plane velocity component is active.  The circular
inclusion/matrix interface is the only weak-trace or cohesive-zone surface.
Post-processing compares axis tracers with the bonded anti-plane circular
inclusion solution and adds VisIt stress/interface diagnostics.
