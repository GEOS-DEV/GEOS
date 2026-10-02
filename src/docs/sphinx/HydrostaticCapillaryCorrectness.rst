Hydrostatic capillary pressure correctness
=========================================

This correction does not change the default world-z coordinate convention of
``HydrostaticEquilibrium``. It repairs two independently incorrect capillary
initialization operations in that existing mode and in any later coordinate
extensions:

* For gas/oil, GEOS stores capillary pressure on the gas slot with the sign
  implied by ``p_gas = P_oil - Pc_gas``. The old initializer populated the oil
  slot instead, causing an incorrect saturation inversion.
* The solver stores a primary reference pressure, not always the pressure of
  the contact-ordered phase. After capillary inversion, it must reconstruct
  ``P = p_present_phase + Pc_present_phase``. A geometric phase may be absent
  before a capillary entry pressure is reached, so the reference is chosen
  from a present phase with actual positive relative permeability.

The active law is evaluated again after inversion. Pressure differences for
all phases with strictly positive saturation and actual positive relative permeability must agree with the integrated phase-pressure
curves, within ``equilibrationTolerance`` plus a pressure-scaled floating-point
roundoff allowance. Nonconverged inversion, invalid reconstructed states, and
incompatible mobile-phase endpoints are errors. Contact pressure matching is
also tightened to the absolute pressure tolerance for active capillarity.
Previously accepted but physically incorrect capillary cases intentionally
change; unaffected no-capillary cases retain their prior behavior.

The actual non-hysteretic relative-permeability model is evaluated, including
three-phase blending. Reported residual minima and small positive saturation
cutoffs are not reliable mobility proxies. A phase excluded as absent or
immobile must satisfy the appropriate one-sided endpoint complementarity.
Only out-of-domain arithmetic roundoff is clipped; small positive phases are
not snapped to zero. History-dependent relative permeability is unsupported.

Acceptance is repeated after the normal component-density/composition/EOS,
saturation, relative-permeability and capillary update, using actual initialized
fields. Temporary phase-pressure records are non-input, non-output and
non-restart data, and are removed afterward. This detects endpoint mobility
introduced by a one-ULP roundtrip or component-density chopping. The authored
``allowLocalCompDensityChopping`` setting is never silently changed. If it
prevents the strict per-mobile-phase pressure contract, the error reports the
measured residual and the setting. Exact-absence fixtures explicitly use zero.
A very small resulting flux may be negligible yet fail this strict contract;
no flux-weighted tolerance or saturation threshold is substituted.

The validation update preserves all unselected rows of solver and constitutive
fields exactly. Its host-resident kernel and explicit memory migration do not
constitute GPU numerical verification.

Datum pressure retains its existing meaning: it anchors the phase chosen by
contact ordering at the datum coordinate (water below the lower contact, oil
between contacts, gas above the upper contact). With capillarity this need
not equal the stored primary pressure. If that contact-ordered phase is absent
in an entry-pressure transition, the datum is a continued phase-pressure
reference, rather than a measurement of a present phase.

Verification and scope
----------------------

Run the independent regression against the corrected executable and an
unmodified pinned baseline::

  python3 scripts/testHydrostaticCapillary.py \
    --geos /path/to/corrected/geosx-wrapper \
    --baseline /path/to/original/geosx-wrapper \
    --output /new/capillary/evidence/directory

The tests reproduce the old gas/oil and three-phase errors, verify analytical
capillary pressures and saturations, entry-pressure clipping, rejection of
incompatible mobile-phase states, and a real closed-domain solver time step.
They compute TPFA mobile-phase potential/flux residuals using actual output
pressure, capillary pressure, density and mobility fields. Physical rock
compressibility fixes the closed-domain pressure nullspace; the tests do not
introduce arbitrary pressure-gauge forcing. The rest-state fixture permits
convergence at Newton iteration zero and verifies actual outputs at both
0 and 1 seconds; it does not force a Newton correction of an already converged
state. Default-sensitive chopping rejection and actual Table/Baker/Stone2
mobility endpoints are covered separately.

These tests use pressure-independent immiscible fluids. This focused repair
alone does **not** establish a self-consistent pressure-dependent EOS/capillary
initialization: the inherited marcher evaluates EOS properties at its chosen
phase pressure, whereas the flow solver evaluates them at its primary
reference pressure. A separate numerical extension is required to remove
that discrepancy. General compositional flash closure and spatially varying
rock-dependent capillary laws must not be inferred to be validated here.
