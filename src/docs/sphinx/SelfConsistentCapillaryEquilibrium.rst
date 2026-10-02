Self-consistent fixed-composition capillary equilibrium
======================================================

This extension adds a coupled hydrostatic initializer for ``DeadOilFluid`` and
``InvariantImmiscibleFluid`` with spatially uniform ``TableCapillaryPressure``.
It supports the gas/oil, oil/water, gas/water and gas/oil/water phase arrangements,
including permutations of the declared phase order, when the authored contact
and capillary constraints admit a physically compatible state.

One initializer's selected union must use one fluid/capillary/relative-
permeability constitutive-model tuple, checked across all MPI ranks in both
coordinate modes. Empty target sets do not constrain that union. Different
model names are conservatively rejected even if their authored parameters
happen to be equivalent; lateral changes in entry pressure or mobility cannot
silently masquerade as one spatially uniform equilibrium.

No coordinate default changes. ``coordinateSystem="elevation"`` remains the
legacy world-z convention. ``gravityAligned`` retains its potential-distance
coordinate and unchanged physical gravity vector. The corrected coupled route
is used in both modes for these supported fluid/capillary combinations. The
no-capillary route is unchanged. Legacy pressure-dependent capillary answers
that used an inconsistent EOS pressure intentionally change.

One pressure for capillarity and EOS
-----------------------------------

GEOS evaluates every phase's EOS properties at the stored primary reference
pressure :math:`P`, and represents phase pressure as
:math:`p_\alpha=P-P_{c,\alpha}`. Evaluating the EOS at a geometrically selected
phase pressure instead produces inconsistent pressure-dependent densities and
component masses, even when the old hydrostatic iteration reports convergence.

The new serial boundary-value solve performs these operations at **every**
hydrostatic fixed-point evaluation:

1. Invert the integrated phase-pressure differences with the actual capillary
   model to obtain saturations, including exact endpoint clipping and a reduced
   gas/water inversion when the three-phase oil saturation vanishes. A bounded
   safeguarded Newton solve uses pressure residuals rather than a fixed endpoint
   snap band or a saturation-step stopping rule
2. Evaluate the forward capillary law at those saturations
3. Evaluate the actual relative-permeability law and select the largest present
   phase with strictly positive relative permeability, then reconstruct :math:`P=p_r+P_{c,r}`
4. Evaluate phase densities and fixed phase compositions at exactly that
   :math:`P` and the authored temperature
5. Solve the trapezoidal hydrostatic equations
   :math:`dp_\alpha/d\zeta=-\rho_\alpha(P,T)|\mathbf g|`

The independent unknown phase pressures at the datum are adjusted by a coupled
shooting solve to satisfy **all** zero-capillary contact pressure equalities.
They are not corrected sequentially with density-independent pressure offsets.
The contact-ordered datum phase pressure stays fixed throughout the solve.
For the legacy coordinate, the signed world-z gravity component is used instead
of :math:`-|\mathbf g|`.

At an actual target cell, the integrated phase-pressure curves are interpolated
and the capillary law is inverted again. The EOS is reevaluated at that cell's
reconstructed primary pressure. Component amounts are formed from its saturation,
its mass/molar phase density and its phase component fractions; interpolated
neighbor-node EOS densities are never substituted for the cell's EOS state.
The solver's configured ``useMass`` basis is respected before building the EOS
wrapper.

Datum and composition semantics
-------------------------------

``datumPressure`` preserves its existing meaning. It anchors the contact-ordered
phase: water below the lower contact, oil between the contacts, and gas above
the upper contact, restricted appropriately for a two-phase model. It does not
universally anchor the stored primary pressure. A phase may be absent before a
positive capillary entry pressure is reached; its datum phase pressure is then
a continued reference curve. The present phase still determines stored primary
pressure through the actual capillary law.

The authored composition tables are validated and supplied to the fixed-
composition EOS. Final overall component fractions are reconstructed to realize
the capillary saturations. They are not promised to equal the input table values.
This is unambiguous for the two supported immiscible fluid families, whose phase
compositions are fixed. General compositional EOS are outside this extension:
a prescribed overall composition, capillary constraints and contact constraints
can be incompatible, and an exact-preservation versus reconstruction contract
requires separate design. This is not a claim that general compositional flash
is inherently underdetermined.

Validity, tolerances and failures
--------------------------------

The final table and cell states must have finite bounded saturations summing to
one, a positive finite primary pressure, positive finite mass/molar densities
and viscosities, and finite normalized phase compositions. Every phase with strictly positive saturation and actual positive relative
permeability must reproduce the integrated phase-pressure differences through
the forward capillary law, matching the solver's phase-presence and mobility
predicates. Absent or actually immobile phases do not incorrectly force the
stored reference pressure to follow their continued phase-pressure curves.
Neither a small saturation cutoff nor the reported residual saturation is a
reliable mobility proxy for all three-phase interpolation laws. Only
out-of-domain saturation roundoff at machine precision is clipped; a small
positive saturation is never zeroed. No saturation error is converted into a
relaxed pressure allowance.

Acceptance is checked again after the normal solver component-density,
composition, EOS, saturation, relative-permeability and capillary update, using
the actual initialized fields. Temporary reference records are removed after
this check and are not input, visualization or restart data. This catches a
residual-saturation endpoint that becomes mobile after even a one-ULP floating-
point roundtrip, as well as phases introduced by component-density chopping.

The authored ``allowLocalCompDensityChopping`` setting is preserved. Exact
phase-absence fixtures use ``allowLocalCompDensityChopping="0"``. If a positive
component-density floor creates a mobile phase whose phase-pressure residual
exceeds the contract, initialization fails with the measured residual and the
setting name. A tiny resulting flux may be numerically negligible and still
fail this strict phase-pressure contract; the implementation neither declares
all tiny phases unphysical nor silently introduces a mobility/saturation cutoff
or a flux-weighted replacement tolerance.

``equilibrationTolerance`` remains an absolute pressure tolerance. The coupled
contact solve and local trapezoid/capillary inversion residuals use that tolerance, with only
pressure-scaled floating-point roundoff allowances for subtraction of large
pressures. ``maxNumberOfEquilibrationIterations`` bounds both the local EOS
fixed-point iterations and the contact iterations. Contact backtracking is
bounded to twenty halvings, as is admissibility backtracking for EOS updates. Nonconvergence, singular contact constraints,
incompatible mobile endpoints and inadmissible EOS states fail explicitly.
There is no automatic tolerance relaxation. An author may need more iterations
or a finer integration table for a difficult admissible case.

For DeadOil, validation inspects the authoritative initialized PVT objects,
including file-backed tables. Every bulk pressure interval must have finite
positive endpoint properties and a nonincreasing formation-volume factor. This
conservative check prevents a coarse integration step from skipping an invalid
unsampled segment. Pressure must be covered by the actual PVT table. An unused
first interval is exempt, preserving the existing surface-condition convention,
but is checked if used. This positional exemption does not independently prove
that arbitrary first-interval data represent a physical surface state.
Water compressibility is nonnegative and water/surface properties must be
finite and physically admissible. No admissible PVT data are invented by
extrapolation or clipping.

In gravity-aligned mode, active capillarity with other fluid families or with
J-functions/spatially or pressure-dependent rock coupling is explicitly
unsupported. History-dependent relative-permeability models are also outside
this initial revision; the non-hysteretic relative-permeability law is evaluated
explicitly, including its three-phase blending. ``initialPhaseName`` with active
capillarity remains unsupported;
use the contact formulation. No chemical or thermal-equilibrium guarantee is
made. The serial path explicitly migrates inputs and target fields to host
memory and marks written fields there. This does not imply GPU numerical
validation; the regression scope is CPU/MPI.

Numerical verification
----------------------

Run the real-solver regression against the preceding capillary-correctness
checkpoint, which contains the common-pressure EOS counterexample::

  python3 scripts/testSelfConsistentCapillary.py \
    --geos /path/to/coupled/geosx-wrapper \
    --baseline /path/to/preceding/capillary-corrected/geosx-wrapper \
    --mpiexec /path/to/mpirun \
    --output /new/self-consistent/evidence/directory

The harness requires Python VTK and writes XML, native meshes, actual VTK output,
solver logs and a numerical evidence summary. It tests pressure-dependent EOS
closure, two/three-phase arrangements and phase permutations, mass/molar basis,
positive entry and disappearance/reference switching, zero gravity, upward and
oblique covariance, the legacy corrected path, native selected element sets
with an empty MPI rank, explicit invalid-input/nonconvergence failures, and a
real closed-domain flow time step.

The EOS closure test is separate from discretization accuracy. The integrated
continuum state need not be an exactly discrete TPFA rest state on a coarse
mesh when density varies with pressure. The tests calculate TPFA phase flux
from actual primary pressure, capillary pressure, mean phase density and upwind
mobility fields. They separately refine the hydrostatic integration table and
the spatial mesh, rather than relabeling truncation error as solver roundoff or
weakening the pressure tolerance to hide it.

Measured CPU regression results
-------------------------------

The focused suite contains 17 groups and 67 actual solver invocations, including
MPI and deliberate rejection cases. In the pressure-dependent gas/oil negative
control, the preceding implementation produced a maximum phase mass flux per
area of :math:`2.2444\times10^{-8}`. The coupled implementation reduces this to
:math:`1.1509\times10^{-12}` on the same spatial mesh, with a closed-form
primary-pressure error of :math:`1.5181\times10^{-6}` Pa. A real one-second flow
step changes pressure by at most :math:`9.89\times10^{-7}` Pa.

In the smooth, stronger-compressibility refinement fixture, halving spatial
spacing changes maximum TPFA phase mass flux per area from
:math:`1.1431\times10^{-10}` to :math:`3.1703\times10^{-11}` to
:math:`8.3279\times10^{-12}`. Water-containing arrangements also include a
primary-pressure derivative change at phase appearance: their interface flux
converges approximately first order, while the corresponding phase-potential
residual converges approximately second order. These finite-mesh errors are
reported separately from local constitutive and hydrostatic convergence.

The selected-set test compares every emitted field outside the selected native
set against an initializer-free control, including an outside pure-phase state
sensitive to component-density chopping. Those output fields are bit-identical.
The verification-only normal state update preserves unselected solver and
constitutive rows using the existing host indexed pack/unpack API.
