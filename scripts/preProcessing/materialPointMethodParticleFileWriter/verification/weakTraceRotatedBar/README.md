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

The weak-trace variant defaults to a small gap-feedback test,
`weakInterfaceTraceGapStabilization=0.01`, with the accumulated-gap correction
velocity capped by `weakInterfaceTraceMaxGapCorrectionVelocity=1.0e-5`.
Override `WEAK_TRACE_GAP_STABILIZATION=0.0` for the stable velocity-only
projection baseline or `0.1` to reproduce the overcorrected failed case.
Override `WEAK_TRACE_MAX_GAP_CORRECTION_VELOCITY=<value>` to test the cap.
