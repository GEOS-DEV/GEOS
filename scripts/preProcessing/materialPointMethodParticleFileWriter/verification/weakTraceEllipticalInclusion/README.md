# weakTraceEllipticalInclusion

Eshelby-type elliptical inclusion benchmark for explicit-surface weak trace
development.

The folder dispatches three variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact groups tied by the weak trace projection.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.

The inclusion is an axis-aligned ellipse embedded in a square matrix and loaded
by a slow remote uniaxial displacement ramp.  The inclusion/matrix interface is
the only weak-trace or cohesive-zone surface.
