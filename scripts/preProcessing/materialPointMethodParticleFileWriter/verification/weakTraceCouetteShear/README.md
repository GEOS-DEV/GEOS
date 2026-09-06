# weakTraceCouetteShear

No-slip parallel-plate Couette shear benchmark for explicit-surface weak trace
development.

The folder dispatches three variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact groups tied by the weak trace projection.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.

The domain is periodic in the flow direction.  The top and bottom boundaries
receive opposite tangential velocities, and the horizontal material interface is
the only weak-trace or cohesive-zone surface.
