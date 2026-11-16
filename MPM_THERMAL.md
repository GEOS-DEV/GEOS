# MPM thermal transport

The unified energy/coupling revision is documented in
[MPM_THERMOMECHANICAL.md](MPM_THERMOMECHANICAL.md). Its material support and restart
restrictions supersede those of the original transport-only extension below.

This extension adds Fourier heat conduction, conservative particle/grid thermal
transfer, interface resistance, particle heating/cooling, and prescribed particle
temperatures to the explicit `SolidMechanics_MPM` solver. It is integrated once per
physical step and works independently of the mechanical PIC/FLIP/XPIC/FMPM choice
and the mechanical Pairwise/PGS/Newton contact iterations.

Apply the cumulative patch to the unmodified source in
`GEOS-feature-crook5-multifield(4).zip`, from the repository root:

```bash
git apply --check geos_mpm_unified_thermomechanics_cumulative.patch
git apply geos_mpm_unified_thermomechanics_cumulative.patch
```

The earlier foundation/operators patches are not prerequisites. Do not apply this
cumulative patch on top of them. Rebuild GEOS after applying it. `Thermal.cpp` is
included in both the full constitutive build and `GEOS_ENABLE_CONSTITUTIVE_MPM_ONLY`.

## Thermal material composition

`Thermal<MechanicalModel, ThermalLaw>` adds a caloric/transport policy without
changing the mechanical kernel type or copying its state arrays. Mechanical
dispatch deliberately uses the original mechanical base; thermal transport
queries the independent `ThermalMaterial` capability. A distinct mechanical
catalog identity preserves named mechanical update/event paths.

Registered mechanical bases are `ElasticIsotropic`, `CeramicDamage`,
`SurfaceInformedPolymer`, `StrainHardeningPolymer`, `Graphite`, `Geomechanics`,
`Hyperelastic`, `HyperelasticMMS`, `Chiumenti`, `VonMisesJ`, `PerfectlyPlastic`,
and `ElasticTransverseIsotropic`. Prefix a base with `Thermal` for constant thermal
properties or `LinearThermal` for temperature-dependent thermal properties.

Only the `ElasticIsotropic` and `PerfectlyPlastic` thermal adapters
currently satisfy the unified energy contract. Other aliases remain registered
but are rejected when thermal transport is enabled. See MPM_THERMOMECHANICAL.md.
All active particle materials must use a supported thermal adapter when transport is enabled.
The background-grid null material can remain ordinary `ElasticIsotropic`.

Example using SI units:

```xml
<ThermalElasticIsotropic
  name="steel"
  defaultDensity="7800"
  defaultBulkModulus="1.6e11"
  defaultShearModulus="8.0e10"
  thermalSpecificHeat="500"
  thermalConductivity="45"
  thermalReferenceTemperature="300"
  defaultDrainedLinearTEC="1.2e-5"
  solverThermalExpansion="1"/>
```

Capacity is specific energy/(mass K); conductivity is power/(length K).
For the constant law, `u=c*(T-Tr)`. Reference thermal energy may be negative.
For the linear law:

```text
c(T) = thermalSpecificHeat + thermalHeatCapacitySlope*(T-Tr)
k(T) = thermalConductivity + thermalConductivitySlope*(T-Tr)
u(T) = c0*(T-Tr) + 0.5*thermalHeatCapacitySlope*(T-Tr)^2
```

The code checks positive capacity, nonnegative conductivity, and nonnegative
absolute temperature at runtime. Cooling beyond the admissible caloric branch is
an error; no thermal energy is silently clipped.

`solverThermalExpansion=1` includes isotropic thermal strain and reversible
thermoelastic energy in the supported corotational small-elastic-strain updates.
It defaults to zero. Other mechanical models require their own stored-energy
adapters before they can participate in the unified thermomechanical solve.

## Solver configuration and contact laws

Put the optional block inside `SolidMechanics_MPM`:

```xml
<Thermal enabled="1"
         defaultContactMode="Resistance"
         defaultContactResistance="1.0e-4"
         gapTolerance="1.0e-6"
         safetyFactor="0.45"
         areaQuadrature="32"
         maxSubsteps="10000"
         diagnostics="1">
  <ContactResistance name="groups01"
                     contactGroups="{ 0, 1 }" resistance="2.0e-4"/>
  <ContactResistance name="steelPolymer"
                     materials="{ steel, polymer }" resistance="5.0e-4"/>
  <ContactResistance name="specificPair"
                     contactGroups="{ 0, 2 }" materials="{ steel, polymer }"
                     resistance="7.0e-4"/>
  <ContactResistance name="excludedPair"
                     contactGroups="{ 1, 2 }" mode="Insulated"/>
  <HeatSource name="heater" particleIDs="{ 12, 19 }"
              power="0.01" startTime="0" endTime="2"/>
  <HeatSource name="cooler" contactGroups="{ 3 }" specificPower="-100"/>
  <Isothermal name="bath" targetRegions="{ ParticleRegion1 }"
              temperature="300"/>
</Thermal>
```

Resistance is **area-specific**, in area K/power (m^2 K/W in SI).
Conductance is `unique interface area / resistance`. Priority is combined
material/group selectors, then group pair, then material pair, then the global
law. Side associations in combined selectors are retained when reversing a pair.
Equally specific matching rules are rejected rather than resolved by input order.
Material names refer to named constitutive instances, not material class names
or particle-region numeric indices.

Modes are `Resistance`, `Perfect`, and `Insulated`. In `Resistance` mode a zero
resistance means exact perfect nodal contact. The default is `Resistance` with
zero resistance. Open interfaces are insulated; `gapTolerance` is a fraction of
the minimum active cell size. Heat exchange is independent of relative velocity,
friction, mechanical contact activation, and whether a contact impulse exists.

Thermal slots retain both the mechanical field and the material identity. Intact
same-group DFG splits remain thermally continuous. Damaged DFG interfaces use
the existing separability classification and the chosen thermal contact law.
Separate materials sharing one mechanical field are thermally bonded unless an
explicit material-pair rule is supplied. Thermal insulation does not change
mechanical contact constraints.

The interface geometry uses opposed reconstructed normals and surface positions;
the fallback estimates a surface from the center-of-mass separation and grid
spacing, as in the current mechanical contact approach. Each eligible binary
patch receives one planar nodal area, clipped to nodal support and the physical
domain. The same area transfers equal/opposite energy on its two sides. Plane
strain uses the represented physical thickness. Linear/CPDI and cubic B-spline
grid bases are supported by the area quadrature. Increase `areaQuadrature` to
check area convergence.

## Particle sources, thermostats, and lifecycle

Conditions accept `particleIDs`, `contactGroups`, `particleColors`, and
`targetRegions`. Selector categories are intersected; list entries within one
category are alternatives. Empty selectors mean all particles. Membership is
material-following, fixed when a particle is first encountered, and copied to
split children. Initial region selectors use bare region names, without the
mesh-body prefix.

Heat-source contributions add:

```text
P = specificPower*currentMass + volumetricPower*currentVolume
    + power*currentMass/originalSelectedParticleMass
```

Positive values heat and negative values cool. `power` is per original selected
particle; subdivision distributes it by child mass, preserving the total.
Sources act directly on selected particles, so P2G does not spread a prescribed
source to unrelated neighbors.

An isothermal particle remains a conducting heat reservoir. Its temperature is
enforced before each thermal projection and after the energy update. The
optional `temperatureRate` gives `temperature + temperatureRate*(t-startTime)`.
Pure uniformly prescribed nodal slots also supply Dirichlet data to the implicit
contact graph. The required reservoir heat is recorded separately, including
sources on thermostatted particles. Conflicting temperatures selecting the same
particle are rejected. Mixed fixed/free particle supports impose particle-level
constraints; check spatial and substep convergence near unresolved reservoirs.

`particleInternalEnergy` is now the sole total specific internal energy.
Temperature recovery excludes stored elastic, thermoelastic, and retained plastic
energy using the supported material contract. Plastic heat is a partition of
stress work, not an extra addition to total energy. Isothermal particles exchange
mechanically generated heat with their reservoirs. The optional lagged mode
and the default local coupled solve are described in MPM_THERMOMECHANICAL.md.
Start a fresh run: old thermal checkpoints are rejected, not silently converted.
Keep thermal condition/material names and their order unchanged for new restarts.

## Discretization, parallel execution, and limits

The thermal stage executes after grid/contact geometry is available, before
particle motion and the constitutive update. Geometry remains frozen over the
thermal substeps. Capacity-weighted P2G, particle-quadrature Fourier assembly, and
energy-increment G2P use the existing effective particle shape maps. Retained
supports are normalized with their derivatives so constant fields and energy
balance survive tiny-capacity filtering and particle/grid-node coincidences.

Bulk conduction is explicit. A bound on the actual capacity-scaled discrete
operator controls substeps, rather than using only the mechanical CFL estimate.
Source cooling and condition start/end times also limit substeps. Contact
exchange solves all eligible nodal pairs together with backward Euler. Perfect
edges are eliminated by union, not a division by zero or a large penalty.

Only owned particles scatter. Additive quantities use ghost-to-owner reduction
and owner-to-ghost replacement. Node owners solve contact; solved energy
increments are replaced onto ghosts, never summed a second time.

This is a **host thermal implementation**, including on HIP/CUDA builds. Explicit
memory transitions preserve compatibility with device mechanical kernels. It is
not a GPU-optimized transport kernel. Dense thermal grid storage scales with
nodes times mechanical fields times materials. Rigid stages use full body-color
identities and include interior particles; their dense thermal storage scales
with body count. Compact/GPU thermal storage is a future performance improvement.

Thermal contact uses the branch's current planar surface reconstruction. It does
not resolve multiple corner normals or reconstruct finite triple-junction patches.
Explicit insulated-pair rules can exclude a known non-interface pair. Radiation,
gas-gap conduction, and long-range transfer across an open gap are outside this
solid-contact law.

`diagnostics=1` reports each transport stage's total internal-energy change, source heat, reservoir
heat, and balance residual. These are globally reduced; compare the residual
against the relevant energy scale and floating-point roundoff.

## ParticleFileWriter

Thermal material aliases use the existing mechanical material generators.
For an existing dictionary-compatible material, set its `model` to the thermal
alias and add thermal properties. The `thermal` PFW option accepts a complete
`<Thermal ...>` XML string or a settings dictionary:

```python
pfw["thermal"] = {
    "defaultContactResistance": 1.0e-4,
    "diagnostics": 1,
    "heatSources": [{"name": "heater", "particleIDs": [12, 19], "power": 0.01}],
    "isothermal": [{"name": "bath", "contactGroups": [0], "temperature": 300}],
    "contactResistance": [
        {"name": "pair", "contactGroups": [0, 1], "resistance": 2.0e-4}
    ],
}
```

Use a consistent unit system. PFW's common mm/us/mg system is not SI:
500 J/(kg K) becomes `5e-4`; 45 W/(m K) becomes `4.5e-8`;
1e-4 m^2 K/W becomes `1e8`; 0.01 W becomes `1e-8`.

## Verification

The C++ numerical test is part of the MPM unit-test build. It can also run without
GEOS third-party libraries:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Werror -Wfloat-equal -pedantic \
  -I src/coreComponents \
  src/coreComponents/physicsSolvers/solidMechanics/unitTests/testMPMThermal.cpp \
  -o /tmp/testMPMThermal
/tmp/testMPMThermal
python3 inputFiles/materialPointMethod/thermalTransport/test_inputs.py
```

Generate small full-solver verification jobs with standard Python:

```bash
python3 inputFiles/materialPointMethod/thermalTransport/generate.py /tmp/thermal-source --case source
cd /tmp/thermal-source
/path/to/geos -i thermal.xml > run.log 2>&1
python3 /path/to/GEOS/inputFiles/materialPointMethod/thermalTransport/check.py .
```

Cases are `source`, `bulk`, `groups`, `dfg`, `materials`, `isothermal`, and
`insulated`. They exercise analytical particle heating, an insulated cosine
diffusion mode, stationary contact, DFG/group equivalence, material rules,
thermostat reactions, and excluded interfaces. Repeat with finer `--cells` and
area quadrature for convergence, and compare one-rank/MPI and host/HIP runs.
The generated mechanical contact is disabled deliberately to verify that
stationary thermal exchange does not depend on mechanical contact impulses.

Validation performed when preparing the patch: standalone C++ numerical checks,
sanitizer checks, Python/XML generation checks, and clean patch application.
The full GEOS executable, its third-party libraries, and HIP runtime were not
available in the patch-generation environment; the generated full-solver jobs
are provided for that validation and have not been represented as executed.
