=========
Viscosity
=========

Viscosity dictates the fluid's resistance to shear flow and fundamentally controls phase mobilities in porous media. The framework provides multiple approaches to evaluate viscosity.

Constant Viscosity Model
------------------------

The simplest approach applies a static, user-defined viscosity for a specified phase, entirely independent of pressure, temperature, or composition. This is recommended solely for idealised simulations, synthetic testing, or when fluid flow is dominated entirely by pressure gradients rather than mobility contrasts.

* **Catalog Name:** Selected by passing ``ConstantViscosity`` into the phase model formulation.
* **Parameters:** ``constantPhaseViscosity``.

Lohrenz-Bray-Clark (LBC) Viscosity Model
----------------------------------------

The Lohrenz, Bray, and Clark (1964) correlation is the industry-standard methodology for calculating the viscosity of compositional hydrocarbon mixtures. 

The model evaluates viscosity in two distinct steps:
1. **Dilute Gas Viscosity:** The viscosity of each pure component at low-pressure (dilute) conditions is estimated as a function of temperature using the Stiel and Thodos (1961) correlations. These pure viscosities are then blended to find the dilute mixture viscosity. The user may select between the Herning-Zipperer, Wilke, or Brokaw mixing rules.
2. **Dense Fluid Contribution:** A residual viscosity, representing the effects of elevated pressure and density, is calculated using a fourth-order polynomial dependent on the reduced density of the phase. This residual is added to the dilute gas viscosity to obtain the final phase viscosity.

This model is strongly recommended for all miscible gas, volatile oil, and general compositional hydrocarbon simulations.

* **Catalog Name:** ``LohrenzBrayClark``
* **Parameters:**
  * ``viscosityMixingRule``: Defines the dilute gas mixing rule (``HerningZipperer``, ``Wilke``, or ``Brokaw``).
  * ``componentCriticalVolume``: Critical volumes required for the reduced density calculation.

Phillips Brine Viscosity Model
------------------------------

For aqueous phases, the Phillips Brine Viscosity model uses empirical correlations from Phillips et al. (1981) to calculate viscosity as a function of temperature, pressure, and salinity. The model relies on tabulated pure water saturation viscosities and scales them using a salinity-dependent multiplier. This model is recommended for saline aquifers and CO2 sequestration.

* **Catalog Name:** ``PhillipsBrine``
* **Parameters:** Tabulated via function managers; relies on the global ``salinity`` parameter.