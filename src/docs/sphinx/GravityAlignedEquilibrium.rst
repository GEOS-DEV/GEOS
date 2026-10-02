Gravity-aligned hydrostatic initialization
=========================================

``HydrostaticEquilibrium`` accepts the optional attribute
``coordinateSystem="gravityAligned"``. The default, ``elevation``, retains
legacy world-z coordinates and the legacy restriction to vertical gravity.
Opting in changes the initializer's coordinate system, **not** the gravity
vector used by the solver.

Coordinate and physical contract
-------------------------------

For a finite, nonzero constant gravity vector :math:`\mathbf g`, define

.. math::

   \mathbf u=-\mathbf g/|\mathbf g|,\qquad
   \zeta=\mathbf u\cdot\mathbf x,\qquad
   \frac{dp_\alpha}{d\zeta}=-\rho_\alpha|\mathbf g|.

``datumElevation``, every ``phaseContacts`` entry,
``elevationIncrementInHydrostaticPressureTable``, and all coordinates of
``temperatureVsElevationTableName`` and
``componentFractionVsElevationTableNames`` use this same scalar potential
*distance*. It has length units (metres in the usual GEOS SI convention), not
energy units. A world-space datum point is encoded as its projection
:math:`\mathbf u\cdot\mathbf x_D`; no independent datum-point attribute is
introduced. For downward z gravity, this coordinate is exactly world z.
With zero gravity the convention is world z and pressure is constant.

Each contact is a constant-potential-distance plane, perpendicular to gravity.
There is no independent contact-normal attribute. A tilted static interface
between fluids of different densities under vertical gravity cannot be
represented by this hydrostatic model without additional physics. Do not
rotate gravity to accommodate a proposed contact, or mix world-z table
coordinates with potential-distance contact coordinates. Input syntax cannot
detect a numerically plausible but incorrectly converted scalar coordinate;
the author must convert every related scalar/table together.

This is pressure-hydrostatic initialization. It does not claim chemical or
thermal equilibrium for arbitrary authored composition or temperature gradients.
A single-phase declaration using ``initialPhaseName`` is rejected in the new
mode if the supplied fluid state creates multiple phases; use reviewed contact
and capillary initialization for such states.
The one-dimensional tables must have finite, strictly increasing coordinates,
finite component fractions in [0,1] summing to one at common tabulation points,
and positive absolute temperatures.

Contacts and capillarity
-----------------------

The existing phase-pressure integration is reused with the gravity component
along the new coordinate. ``equilibrationTolerance`` is an absolute pressure
tolerance in Pa. In gravity-aligned mode it also bounds the phase-pressure
mismatch at every contact; nonconvergence is an error rather than an accepted
partial initialization.

Contacts are **zero-capillary-pressure reference surfaces**. With an active
capillary-pressure law, its inversion determines the saturation transition
zone; a reference surface is not necessarily a sharp saturation boundary.
With no capillarity, contacts define sharp cell-centre phase partitions. A
cell centre on a contact belongs to the phase above it in potential distance.
The projection uses compensated sums and product residuals to reduce cancellation
of large world-coordinate terms. The upper-phase tie band is bounded by
``8*epsilon*sum(abs(u_i*x_i))`` for represented geometry plus
``2*epsilon*max(abs(zeta),abs(contact))`` for scalar comparison. Distinct
positive contact gaps at or below four times the global selected-cell bound
are rejected as too close to resolve; a tolerance must not swallow a physical
thin layer. This convention cannot restore precision absent from the mesh.
It does not change existing mesh-quality/face-area tolerances for sliver cells.
No fractional-volume treatment of crossing cells is introduced. In zero
gravity, a no-capillary sharp partition can coexist with constant pressure;
with capillarity, the active constitutive law still determines saturation.

Legacy mode preserves world-z scalar coordinates and its interpretation of
upward-z and zero gravity. Unaffected no-capillary decks retain their previous
results. This change also repairs proven capillary initialization bugs in
**both** coordinate modes: gas/oil target pressure must use the gas slot with
its native sign, and the stored primary pressure must account for the active
capillary model's reference phase. These cases intentionally produce corrected
results rather than retaining physically nonhydrostatic legacy outputs.

For active capillarity, phase curves are matched to the absolute pressure
tolerance in both modes. After inversion, primary pressure is reconstructed
from a present phase and the evaluated capillary pressure. Mobile-phase
pressure consistency is checked using strictly positive saturation and the
actual positive relative permeability, including three-phase blending. No
small-saturation or reported residual-saturation proxy substitutes for actual
mobility. The normal solver EOS/composition/capillary update is then performed
and the realized mobile-phase pressure residual is checked again. A state
that becomes mobile through a floating-point endpoint roundtrip or authored
component-density chopping fails explicitly when it exceeds the pressure
contract. The authored ``allowLocalCompDensityChopping`` setting is preserved;
exact-absence configurations may require setting it to zero. This is a strict
phase-pressure criterion, not a claim that every tiny resulting flux matters
physically. See :doc:`HydrostaticCapillaryCorrectness`.

Native set selection and supported solvers
-----------------------------------------

``setNames`` is optional and defaults to ``{all}``. In gravity-aligned mode it
may name existing authoritative GEOS element sets. Native membership is used
unchanged; a Workbench preview or node set is not an exact element-set export.
Legacy mode accepts only ``{all}``.

A selected set must have target cells globally on this solver's target regions.
Missing/empty selections and overlapping *distinct* initializers are errors.
Overlapping sets of the same initializer are allowed. Cells outside its target
union are left unchanged. Internal tables are scoped by the full subregion
path. One initializer cannot span different fluid-model names; such regions
need separate, non-overlapping initializers and physically reviewed interface
conditions. They are not automatically pressure-matched across separate
initializer definitions.

Gravity-aligned **active-capillary** initialization requires
``phaseContacts``, a spatially uniform ``TableCapillaryPressure`` model,
and a supported fixed-composition fluid. The self-consistent numerical
extension supports ``InvariantImmiscibleFluid`` and ``DeadOilFluid`` and
matches capillarity and EOS at the same primary pressure. See
:doc:`SelfConsistentCapillaryEquilibrium` for its capability version, closure,
convergence and constitutive scope. Without that separate extension, only
the pressure-independent invariant-fluid active-capillary path is certified.
Single-phase flow without capillary coupling supports its existing
compressible fluid models.

The implementation is provided by ``SinglePhaseBase`` and
``CompositionalMultiphaseBase`` and their supported constitutive dispatches.
Other flow-solver families reject this mode instead of silently ignoring it.
CPU/MPI regression coverage includes SinglePhaseFVM (isothermal and thermal)
and CompositionalMultiphaseFVM; GPU and every derived solver have not been
independently validated by this change.

Machine-readable discovery
--------------------------

``--capabilities --format=json`` exposes
``gravity-aligned-hydrostatic-initialization: true`` and
``hydrostaticInitializationVersion: 1``. Authoring clients must check both,
and the authoritative schema enum, before generating the new mode.
``native-set-hydrostatic-initialization`` describes selection of existing GEOS
sets. ``arbitrary-plane-phase-initialization`` and
``generated-set-phase-initialization`` remain false.
Clients requiring compressible active-capillary initialization must separately
check ``self-consistent-capillary-initialization`` and
``selfConsistentCapillaryInitializationVersion: 1`` and the reported fluid/model
scope. Coordinate support alone does not certify a selected constitutive model.

The global input catalog exposes the actual ``coordinateSystem`` enum/default
and ``setNames`` metadata. Its schema digest changes when this feature is
added. Structured input unit metadata remains unknown where the wrapper
registry does not declare it; descriptions here do not turn unknown catalog
units into inferred declarations.

Regression command
------------------

From the repository root, run::

  python3 scripts/testGravityAlignedEquilibrium.py \
    --geos /path/to/updated/geosx-wrapper \
    --baseline /path/to/pre-initializer/geosx-wrapper \
    --mpiexec /path/to/mpirun \
    --output /new/evidence/directory

The Python environment requires VTK. Tests use actual solver outputs rather
than reproducing the initializer as a mock: legacy exact equivalence,
analytic unequal-density two/three-phase columns and contact continuity,
rigid rotations and translated datum/tables, temperature/composition,
capillary transitions, zero gravity, singleton target coordinates,
single-phase/thermal pressure, selected sets, rejection paths, a converged
flow step, and two-rank MPI covariance.
