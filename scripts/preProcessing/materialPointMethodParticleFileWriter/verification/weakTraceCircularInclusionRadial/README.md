# weakTraceCircularInclusionRadial

Centered circular inclusion benchmark for explicit-surface weak trace
development under biaxial radial-type loading.

The folder dispatches three variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact groups tied by the weak trace projection.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.

The inclusion/matrix interface is the only weak-trace or cohesive-zone surface.
The first version uses a circular finite specimen loaded by the borehole-style
virtual-fluid `ConfiningPressure` event so the outer curved boundary is loaded
by traction rather than by an affine square-box stretch.  Post-processing
compares axis tracers with the bonded composite-cylinder plane-strain solution
and adds VisIt stress/interface diagnostics.
