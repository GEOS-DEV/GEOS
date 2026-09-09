# weakTraceRotatedBar

Rotated two-material MPM verification for local mixed-cell stress error. The
same diagonal plane-strain tensile coupon is dispatched through several
variants: a single-field baseline, a false-elastic cohesive-zone reference, and
prescribed weak-interface projection between contact-group fields. DFG-surface
and adaptive stress-residual variants remain placeholders.

Only the internal material interface is surface flagged in the CZ and mortar
variants. Exterior and boundary-adjacent faces are intentionally unflagged.
Use `./runTest --geometry-only --variant traceContactGroups` to inspect the
PFW-generated surface flags, outward normals, and surface-position vectors
before running GEOS.

For the CZ and mortar variants, the case maps prescribed interface normals with
`explicitSurfaceNormalInfluence=25.0` by default. Override it with
`WEAK_TRACE_EXPLICIT_NORMAL_INFLUENCE=<value>` when comparing normal-mapping
sensitivity.

The current experimental implementation reconstructs a connected C0,
piecewise-linear interface from particle-carried oriented surface samples. It
places one P1 mortar-hat constraint at each reconstructed interface vertex,
integrates the B-spline or Q1 trace on the incident segments with two-point
Gauss quadrature, and enforces the overlapping constraints with alternating
forward/reverse local mass-weighted projection sweeps. It does not assemble a
global constraint matrix and does not alter particle mass, volume, density, or
EOS state.

The defaults are `weakInterfaceTraceProjectionIterations=30` and
`weakInterfaceTraceProjectionScale=0.6`. The iteration count is the number of
local projected Gauss--Seidel sweeps, not a Krylov iteration count. The
prototype intentionally disables displacement/gap feedback:
`weakInterfaceTraceGapStabilization=0` and
`weakInterfaceTraceMaxGapCorrectionVelocity=0` are required. Inspect the solver
log for reconstruction and constraint counts, the relative mortar residual,
post-projection velocity jump, kinetic-energy change, and momentum balance.

Cubic single-point B-spline particles are the recommended default. Linear CPDI
can use the same reconstructed mortar with the underlying Q1 grid trace.
Linear single-point MPM is not supported because an interface tie does not cure
its cell-crossing instability.
