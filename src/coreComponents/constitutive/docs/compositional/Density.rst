=======
Density
=======

The accurate determination of phase density is critical for modelling fluid flow, as it directly impacts volumetric calculations, gravity drainage, and phase mobilities.

Compositional Density Model
---------------------------

The compositional density model derives the molar volume of a phase using the compressibility factor, :math:`Z`, computed directly from the cubic Equation of State (e.g., Peng-Robinson). The molar density, :math:`\rho_{molar}`, is given by:

.. math::
    \rho_{molar} = \frac{P}{Z R T}

Because standard cubic equations of state notoriously under-predict liquid phase densities, a Peneloux volume shift correction is mathematically applied. The corrected molar volume, :math:`v_{corrected}`, is adjusted by a sum of component-specific volume shift parameters, :math:`c_i`, weighted by the phase composition:

.. math::
    v_{corrected} = v_{EOS} - \sum_{i=1}^{N_c} x_i c_i

The molar density is then inverted to mass density using the molecular weight of the phase mixture. This model is recommended for all hydrocarbon phases.

* **Catalog Name:** ``CompositionalDensity``
* **Parameters:** * ``componentVolumeShift``: Component-specific volume shift parameters.

Phillips Brine Density Model
----------------------------

For aqueous phases containing dissolved salts, cubic equations of state often struggle to capture the complex ionic interactions. The Phillips Brine Density model employs a robust empirical correlation developed by Phillips et al. (1981).

This model calculates the mass density of the brine as a direct polynomial and exponential function of pressure, temperature, and salinity (molality of dissolved salts). For low salinities, the model interpolates smoothly towards pure water properties to maintain physical continuity. This is highly recommended for geothermal reservoirs or saline aquifers.

* **Catalog Name:** ``PhillipsBrineDensity``
* **Parameters:**
  * ``salinity``: Brine salinity.
  * ``saltMolarWeight``: The molar weight of the salt component.

Immiscible Water Density
------------------------

A simplified exponential density model is available for pure, immiscible water phases. It evaluates density based on isothermal compressibility and thermal expansion from a reference state:

.. math::
    \rho = \rho_{ref} \exp(c_w (P - P_{ref})) \exp(-\alpha_w (T - T_{ref}))

* **Catalog Name:** ``ImmiscibleWaterDensity``
* **Parameters:** ``waterReferencePressure``, ``waterReferenceTemperature``, ``waterDensity``, ``waterCompressibility``, ``waterExpansionCoefficient``.