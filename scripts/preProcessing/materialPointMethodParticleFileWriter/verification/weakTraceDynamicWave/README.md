# weakTraceDynamicWave

Quasi-one-dimensional dynamic wave transmission benchmark for explicit-surface
weak trace development.

The folder dispatches three variants:

- `singleField` - ordinary single-field mixed-material mapping.
- `traceContactGroups` - two contact groups tied by the weak trace projection.
- `falseElasticCZ` - a stiff elastic cohesive-zone tie used as a reference.

The left boundary is driven by an F-table approximation to a small Gaussian
velocity pulse.  The right computational boundary is at `x = 0`, so the global
F-table velocity field prescribes zero normal velocity there while imposing the
requested velocity at the left boundary.  The material interface is inside the
bar and is the only surface flagged for weak-trace or cohesive-zone treatment.
Profile-history output writes x-line `velocityX` and `stressXX` CSV files so the
post-process can compare numerical and analytical wave profiles at several
times before and after the pulse reaches the interface.
