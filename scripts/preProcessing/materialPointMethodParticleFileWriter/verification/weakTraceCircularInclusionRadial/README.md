# weakTraceCircularInclusionRadial

Centered circular inclusion benchmark for explicit-surface weak trace
development under biaxial radial-type loading.

The folder dispatches four variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact-group velocity fields tied by the reconstructed C0 P1 vertex mortar.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.
- `primalEnrichment` - the earlier single-field primal ridge prototype retained for comparison.

The inclusion/matrix interface is the only weak-trace or cohesive-zone surface.
The first version uses a circular finite specimen loaded by the borehole-style
virtual-fluid `ConfiningPressure` event so the outer curved boundary is loaded
by traction rather than by an affine square-box stretch.  Post-processing
compares axis tracers with the bonded composite-cylinder plane-strain solution
and adds VisIt stress/interface diagnostics.

The C0 mortar variant defaults to cubic single-point B-spline particles.  Set
`WEAK_TRACE_PARTICLE_TYPE=2` to exercise the same interface reconstruction and
mortar projection with CPDI and the underlying Q1 grid trace.  Particle mass,
volume, density, and EOS state are not altered by the mortar treatment.  The
current PFW `WEAK_TRACE_PPC` control provides uniform particle refinement; a
future interface-local refinement option can be added independently if the
surface-marker resolution proves limiting.
