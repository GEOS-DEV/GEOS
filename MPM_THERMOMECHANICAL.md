# Unified MPM thermomechanical energy and local temperature/plasticity coupling

This revision replaces the two independently advanced particle energy ledgers.
`particleInternalEnergy` is the authoritative **total specific internal energy**.
The thermal stage still precedes the constitutive stage; it supplies a thermal
predictor, not the final particle temperature. The constitutive stage solves the
local stress/temperature feedback before accepting the end-of-step state.

## Supported constitutive contract and scope

The implemented contract is constant-stiffness isotropic elasticity with
corotational small elastic strains and additive plastic/thermal strain increments.
It is enabled for the `Thermal` and `LinearThermal` versions of:

* `ElasticIsotropic`;
* `PerfectlyPlastic`.

These retain the existing mechanical kernels. The adapter queries their elastic
stiffness, not the elastoplastic algorithmic tangent. `thermalYieldSoftening`
adds an optional temperature-dependent yield strength to the perfectly plastic J2 model so the
same-step temperature/yield feedback can be represented explicitly.

Other previously registered thermal aliases still exist in the catalog, but
**are rejected when thermal transport is enabled**. In particular, this revision
does not supply a thermodynamic adapter for anisotropic elasticity, polymers with
temperature-dependent moduli, damage/fracture, Graphite/Geomechanics, EOS models,
or finite-strain hyperelasticity, or the existing `VonMisesJ` logarithmic-pressure
implementation. Applying the quadratic elastic energy used here
to those models would misidentify stored energy and dissipation. They require
model-specific free-energy, stored-history, and dissipation implementations.
This patch is not a claim that all original aliases now implement complete
thermomechanics. Mechanical runs with thermal transport disabled retain their
previous paths.

The supported formulation uses reference density and reference-volume mechanical
work. It is a corotational small-strain energy approximation, not an exact
finite-strain Cauchy-stress energy formulation. Large volumetric strains require
an appropriate work-conjugate finite-strain material adapter. A converged local
solve also does not cure mesh-dependent localization in a softening material.
Frictional, cohesive, damage, and numerical-damping heat sources remain outside
this plastic-work coupling.

## Meaning of the energy state

For the supported constant-stiffness model, let `S = inverse(C)`, `sigma` be the
stress in the corotated material frame, `rho0` the reference density, `alpha` the
isotropic linear thermal expansion coefficient (zero if disabled), and `h(T)`
the existing caloric policy:

```
u(T, sigma, retained) = h(T) + u_el(sigma) + T*a(sigma) + retained
u_el(sigma)           = 0.5 * sigma:S:sigma / rho0
a(sigma)              = alpha * trace(sigma) / rho0
```

The `T*a` term is the thermoelastic contribution obtained from `u = psi + T*s`
with a quadratic elastic Helmholtz potential and thermal eigenstrain. It provides
reversible thermoelastic heating/cooling, in addition to irreversible plastic
heating. `h(T)` remains the constant/linear caloric relation supplied by the
templated thermal policy. Its reference can be negative. Total internal energy
is never clipped at zero.

Only `particleInternalEnergy` is the total energy. Additional fields are:

* `particleRetainedPlasticEnergy`: the accumulated unconverted fraction of
  plastic work; a constitutive history component needed to recover temperature;
* `particlePlasticWork`: cumulative specific plastic work, a diagnostic;
* `particleThermalReservoirEnergy`: cumulative specific energy supplied by a
  prescribed-temperature reservoir, also a diagnostic.

Do **not** add these fields to `particleInternalEnergy` when computing system
energy. They are components or cumulative transfers already represented in it.
There is no persistent `particleThermalEnergy` field.

At frozen stress/history, the transport energy derivative is `c(T) + a(sigma)`.
The transport stage uses this capacity and inverts the full state relation after
adding heat. Elastic and retained energy are not interpreted as heat. The local
mechanical energy solve subsequently includes the stress dependence of the
energy state. Nonpositive frozen-state capacity is an error.

## Mechanical work, plastic work, and temperature feedback

All tensors used for work are first placed in the same corotated frame. Voigt
strain uses engineering shear; no extra shear factors appear in stress-dot-strain
work. The existing kernel receives all three normal thermal strains, including
in plane strain.

```
sigma_bar       = 0.5 * (sigma_n + sigma_new)
deps_p          = deps - alpha*(T_new - T_n)*I - S:(sigma_new - sigma_n)
dw              = sigma_bar:deps / rho0
dwp             = sigma_bar:deps_p / rho0
retained_new    = retained_n + (1 - beta)*dwp
```

This infers the actual plastic increment from constant elastic compliance and the
returned stress; it also works for `PerfectlyPlastic`, which does not store a
plastic strain tensor. Substantially negative plastic work is rejected. Tiny
negative roundoff is set to zero. The stress-work quadrature is midpoint: a
large increment crossing first yield still requires temporal refinement.

For the default `mechanicalCoupling="Coupled"`, the solver finds the particle
endpoint temperature satisfying

```
u(T_new, sigma_new(T_new), retained_new(T_new)) = u_predictor + dw(T_new)
```

The stress update is reevaluated at each trial temperature, so a yield surface
changed by plastic heating is included in the same increment. The safeguarded
scalar energy solve brackets the root and uses secant/bisection steps. It must
satisfy an energy residual; a small temperature change alone is not convergence.
Failure to bracket/converge stops with the material and local particle index;
reduce the physical timestep. This implementation does not automatically roll
back the entire MPM physical step and retry it.

Every trial restores an independent indexed byte snapshot of **all particle-sized
constitutive arrays**, not only old stress. This includes plastic strain, yield
strength, new/old stress, and adapter history. The accepted temperature is replayed
once from that snapshot. The transaction restores only that particle and never
rewinds neighboring particles. `Wrapper::clone()` is intentionally not used as a
deep copy: in this repository it aliases underlying data.

The first law adds mechanical work once. `beta*dwp` is a partition of that work,
not a second source added on top of it. The legacy pre/post stress-work routine
is bypassed whenever thermal transport is enabled.

## Optional lagged coupling

`mechanicalCoupling="Lagged"` evaluates the return at the thermal predictor,
then recovers the final particle temperature from the energy balance while
holding the accepted stress/history fixed. Heating is present at the endpoint;
its yield-surface effect is deferred. Thermal expansion associated with that
additional temperature change is also carried to the next stress update using
the stored temperature consumed by mechanics.

This is an explicit splitting approximation for comparison and weak feedback.
A mechanical CFL bound and a thermal diffusion bound do not guarantee coupling
accuracy. Refine the physical timestep and compare stress, plastic work,
temperature, and localization. Prefer `Coupled` when thermal softening is strong.
The local coupled solve does not require a global implicit momentum solve.

## Isothermal particles and conservation diagnostics

The transport stage saves which material-following particles were thermostatted
on the final subinterval. Their constitutive update uses the prescribed endpoint
temperature. The required reservoir increment is

```
du_bath = u(target, new_stress, new_retained) - u_predictor - dw
```

For ordinary plastic heating this is negative: the bath removes heat. Sources,
transport, and mechanical reservoir exchange are all recorded. A prescribed
isothermal history does not discard plastic energy silently.

`diagnostics="1"` adds a globally reduced `MPMThermomechanics` line with mechanical
work, converted plastic heat, reservoir energy, and the maximum local specific
energy residual. The transport-stage residual remains a separate check of its
conservative transfer. Across the full step, compare total internal energy change
with external heat + mechanical work + reservoir exchange; internal contact and
conduction transfers cancel. No global kinetic/internal-energy conservation claim
is made for the existing explicit grid dynamics and its stress-work discretization.

## Inputs

Existing heat sources, particle selectors, isothermal histories, global/pair
contact resistance, DFG contact handling, and MPI thermal grid synchronization
remain available.

```xml
<Thermal enabled="1" mechanicalCoupling="Coupled"
         couplingTolerance="1e-10" maxCouplingIterations="100"
         defaultContactMode="Insulated" diagnostics="1"/>

<ThermalPerfectlyPlastic name="solid"
    defaultDensity="7800" defaultBulkModulus="1.6e11"
    defaultShearModulus="8e10" defaultYieldStress="4e8"
    thermalSpecificHeat="500" thermalConductivity="45"
    thermalReferenceTemperature="300"
    plasticHeatFraction="0.9" thermalYieldSoftening="0.002"
    solverThermalExpansion="1" defaultDrainedLinearTEC="1.2e-5"/>
```

These values demonstrate syntax, not a calibrated material. The yield law is
`yield(T) = referenceYield * exp[-thermalYieldSoftening*(T-Tr)]`.
`thermalYieldSoftening` has units `1/K`, defaults to zero, and is supported only
by the `PerfectlyPlastic` adapter. The original per-particle yield field is captured once as
the reference and follows particles through restart/migration/subdivision.
`plasticHeatFraction` is in `[0,1]` and defaults to `0.9`; it is a modeling choice,
not a universal material constant. A value of zero stores all plastic work;
one converts all of it to heat.

`couplingTolerance` is a relative specific-energy tolerance, scaled by
`max(1, abs(c(T_predictor)*max(1,T_predictor)))`. `maxCouplingIterations` bounds
bracketing and root iterations separately. Extremely large softening can remove
an admissible root or create multiple branches; reduce the timestep and examine
the material law rather than accepting a failed return.

PFW passes through the new material attributes and these solver dictionary keys:

```python
pfw["thermal"] = {
    "mechanicalCoupling": "Coupled", "couplingTolerance": 1e-10,
    "maxCouplingIterations": 100, "diagnostics": 1,
}
```

Use a consistent unit system for modulus, density, energy, heat capacity, power,
and conductivity. The two new material inputs are dimensionless and `1/K`,
respectively, so they do not acquire length/mass/time conversion factors.

## Restart and integration boundaries

New fields and the reference-yield field are restartable particle/constitutive
state and use the existing indexed particle transfer/subdivision mechanisms.
New checkpoints identify energy formulation version 2. **Old thermal checkpoints
are rejected** because their caloric energy and independent work ledger cannot
be relabeled as a total energy without a material-specific conversion. Start a
fresh simulation after applying the patch. No checkpoint converter is supplied.

The thermal stage still runs at 13.5; local thermomechanics runs through the normal
stress update at stage 17. Particle temperature and material history are committed
before solver dependency/output updates. The local implementation runs on the
host, with explicit wrapper moves before/after host work, just like the existing
reference thermal operator. It is not a device-optimized constitutive solver;
HIP memory movement and MPI integration require validation in a full GEOS build.

## Verification and code map

`src/coreComponents/constitutive/solid/ThermomechanicalLaw.hpp` implements the
stored-energy components, work partition, elastic-compliance inversion, and local
energy root helper. `Thermal.hpp` adds the material capabilities and perfectly plastic J2 yield law.
`SolidMechanicsMPMThermomechanics.cpp` implements constitutive rollback/replay,
local solving, thermostat exchange, and diagnostics. `SolidMechanicsMPMThermal.cpp`
now transports total energy and uses frozen-state inversion. The standard material
and solver registration paths expose the inputs; both CMake lists include the new
sources/tests.

The standalone thermomechanical test uses an independent radial-return material
with the production energy/root helpers. It checks elastic load/unload, shear
work factors, plastic heat fractions 0/0.5/0.9/1, energy balance, final-temperature
yield consistency, timestep refinement of lagged coupling, nonlinear heat capacity,
isothermal reservoir heat, reversible thermoelastic heating, stored-energy
inversion, negative reference energies, and solve failures.

```
g++ -std=c++17 -O2 -Wall -Wextra -Werror -Wfloat-equal -pedantic \
  -I src/coreComponents \
  src/coreComponents/physicsSolvers/solidMechanics/unitTests/testMPMThermomechanics.cpp \
  -o /tmp/testMPMThermomechanics
/tmp/testMPMThermomechanics
```

At packaging: 940 thermomechanical numerical checks, the existing 743 thermal
checks, and three Python input/formatter tests passed. The thermomechanical test
also passed AddressSanitizer/UndefinedBehaviorSanitizer with leak detection
disabled because LeakSanitizer cannot inspect processes in this container. These
do **not** exercise
the GEOS wrapper serialization, full material dispatch, checkpoint I/O, particle
migration/subdivision, MPI, or HIP path. The full GEOS executable/TPL build is not
available in this workspace; no full build or end-to-end run is claimed.

Before production use, build with the target host configuration and run a
homogeneous adiabatic plasticity test (with and without softening), an isothermal
plasticity test, elastic thermal expansion/compression tests, and a same-version
restart. Compare one/multiple ranks, host/HIP, and subdivided/non-subdivided runs.
Check timestep convergence even when every nonlinear solve converges.
