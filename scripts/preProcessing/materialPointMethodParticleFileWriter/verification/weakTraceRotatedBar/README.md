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
