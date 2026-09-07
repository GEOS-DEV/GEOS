# weakTraceRotatedBar

Rotated two-material MPM verification for local mixed-cell stress error.  The
same diagonal plane-strain tensile coupon is dispatched through several
variants: a single-field baseline, a false-elastic cohesive-zone reference, and
prescribed weak-interface trace projection between contact groups.  DFG-surface
and adaptive stress-residual variants are reserved as placeholders.

Only the internal material interface is surface flagged in the CZ and trace
variants.  Exterior and boundary-adjacent faces are intentionally unflagged.
Use `./runTest --geometry-only --variant traceContactGroups` to inspect the
PFW-generated surface flags, normals, and surface-position vectors before
running GEOS.

For the CZ and weak-trace variants, the case maps prescribed interface normals
with `explicitSurfaceNormalInfluence=25.0` by default.  Override with
`WEAK_TRACE_EXPLICIT_NORMAL_INFLUENCE=<value>` when comparing normal-mapping
sensitivity.

The current experimental weak-interface implementation uses one canonical
constraint per cut cell (plus nodal constraints for conforming or sliver cuts)
and solves the coupled mass-metric projection with PCG.  The default
`weakInterfaceTraceProjectionIterations=100` is therefore the maximum PCG
iteration count, not a number of repeated Jacobi projection passes.

This prototype intentionally disables displacement/gap feedback:
`weakInterfaceTraceGapStabilization=0` and
`weakInterfaceTraceMaxGapCorrectionVelocity=0` are required.  Use
`WEAK_TRACE_PROJECTION_SCALE=1.0` for the full coupled projection and inspect
the solver log for PCG convergence, post-projection velocity jump, and kinetic
energy change.
