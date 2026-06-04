Immiscible water flash
----------------------

The immiscible water flash is an efficient three-phase equilibrium model designed for systems containing an aqueous phase alongside hydrocarbon liquid and vapour phases. Instead of solving a fully coupled three-phase equation of state problem, this model vastly simplifies the phase split by assuming water forms a strictly pure aqueous phase and is completely immiscible in the hydrocarbon phases, while hydrocarbons are entirely insoluble in the aqueous phase.

Methodology
~~~~~~~~~~~

In this formulation, the mixture composition is explicitly separated into a water component and a hydrocarbon mixture before any thermodynamic equilibrium calculations occur. The water component is identified using its assigned component name (e.g., ``H2O``). 

Given a total feed mole fraction :math:`z_i` for each component, the mole fraction of the aqueous phase, :math:`V_{aq}`, is determined to be exactly equal to the total feed mole fraction of the water component:

.. math::
    V_{aq} = z_{w}

The composition of the aqueous phase is fixed at 100% water. The remaining feed defines the total hydrocarbon mole fraction, :math:`z_{hc}`:

.. math::
    z_{hc} = 1 - z_{w}

The feed composition of the remaining hydrocarbon mixture is then normalised:

.. math::
    z_{i,hc} = \frac{z_i}{z_{hc}}

for all non-water components, while the normalised hydrocarbon feed fraction for water is set to zero.

Equilibrium calculation
~~~~~~~~~~~~~~~~~~~~~~~

During the simulation step, if the total hydrocarbon fraction :math:`z_{hc}` is negligibly small, the model bypasses the flash and assigns the system as a single-phase aqueous fluid.

If hydrocarbons are present, the model performs a rigorous two-phase negative flash exclusively on the normalised hydrocarbon mixture (:math:`z_{i,hc}`) using the chosen equations of state for the liquid and vapour phases. 

This embedded flash calculation determines the intermediate hydrocarbon vapour fraction, :math:`V_{V,hc}`, and the corresponding compositions of the hydrocarbon liquid (:math:`x_i`) and hydrocarbon vapour (:math:`y_i`) phases using the standard Rachford-Rice equation and iterative fugacity updates.

Once the two-phase hydrocarbon flash converges, the overall phase fractions for the entire three-phase system are scaled back using the total hydrocarbon mole fraction:

.. math::
    V_{V} = z_{hc} V_{V,hc}

.. math::
    V_{L} = z_{hc} (1 - V_{V,hc})

The derivatives of the phase fractions and compositions with respect to pressure, temperature, and total composition are analytically transformed using the chain rule to account for this scaling and normalisation.

Recommendations and pitfalls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The immiscible water flash is highly computationally efficient and is recommended for multi-phase reservoir simulations where water solubility in the hydrocarbon phases (and vice versa) is negligible. By avoiding a full three-phase fugacity solver, it drastically reduces simulation time and improves non-linear solver stability.

Pitfalls: Because this model strictly enforces total immiscibility, it cannot accurately model scenarios where mutual solubility is physically significant, such as high-pressure CO2 injection into saline aquifers (where CO2 dissolves heavily into brine, and water vaporises into the gas stream). For such complex mutual solubility problems, using a two-phase flash paired with the Soreide-Whitson equation of state might be required.

Model parameters
~~~~~~~~~~~~~~~~

The immiscible water flash fluid models are assigned using specific XML blocks that pair the three-phase flash model with distinct viscosity and density models for the hydrocarbon and aqueous phases.

Here is an example demonstrating how to specify the parameters for a three-phase system using the immiscible water flash paired with the Lohrenz-Bray-Clark viscosity model for hydrocarbons and a distinct constant-compressibility formulation for the pure water phase:

.. code-block:: xml

    <CompositionalThreePhaseLohrenzBrayClarkViscosity
        name="fluid"
        phaseNames="{ liquid, vapour, aqueous }"
        equationsOfState="{ PengRobinson, PengRobinson, PengRobinson }"
        componentNames="{ CH4, C10, H2O }"
        componentMolarWeight="{ 1.604e-02, 1.422e-01, 1.801e-02 }"
        componentCriticalPressure="{ 4.599e+06, 2.103e+06, 2.206e+07 }"
        componentCriticalTemperature="{ 1.905e+02, 6.177e+02, 6.471e+02 }"
        componentAcentricFactor="{ 1.142e-02, 4.880e-01, 3.443e-01 }"
        waterReferencePressure="1.01325e5"
        waterReferenceTemperature="293.15"
        waterDensity="998.2"
        waterViscosity="1.002e-3"
        waterCompressibility="4.5e-10" />

Parameters
~~~~~~~~~~~~~~~~

* ``equationsOfState``: List specifying the EoS type for each phase (liquid, vapour, aqueous). Note that the aqueous EoS entry is largely ignored by the pure-water density model but must be provided for list alignment.
