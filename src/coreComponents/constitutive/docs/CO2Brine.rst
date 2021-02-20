.. _CO2-EOS:

##################################################################################
CO2-brine model 
##################################################################################


Summary
=================================

The CO2-brine model currently implemented in GEOSX includes two chemical components (CO2 and water (w)) that are transported by one or two fluid phases (the gas (g) phase and the brine (l) phase).
The water component is only present in the brine phase, while the CO2 component can be present in the gas phase as well as in the brine phase.
Thus, considering the phase component fractions, ::math:`y_{c,p}` (i.e., the fraction of the mass of phase ::math:`p` represented by component ::math:`c`) the following partition matrix determines the component distribution within the two phases:

.. math::
    \begin{bmatrix}
    y_{CO2,g} & y_{CO2,l} \\
         0 & 1            \\
    \end{bmatrix}

The update of the fluid properties is done in two steps:

1) We compute the phase fractions (::math:`\nu_p`) and phase component fractions (::math:`y_{c,p}`) as a function of pressure (::math:`p`), temperature (::math:`T`), component fractions (::math:`z_c`), and a constant salinity.

2) We compute phase densities (::math:`\rho_p`) and phase viscosities (::math:`\mu_p`) as a function of pressure, temperature, the updated phase component fractions, and a constant salinity.

Once the phase fractions, phase component fractions, phase densities, phase viscosities--and their derivatives with respect to pressure, temperature, and component fractions--have been computed, the :ref:`CompositionalMultiphaseFlow` proceeds to the assembly of the accumulation and flux terms.
Note that the current implementation of the flow solver is isothermal and that the derivatives with respect to temperature are therefore not used.

The models that are used in steps 1) and 2) are reviewed in more details below.

Step 1: Computation of the phase fractions and phase component fractions
========================================================================

At initialization, GEOSX constructs a two-dimensional table storing the values of CO2 solubility in brine as a function of pressure, temperature, and a constant salinity.
The user can parameterize the construction of the table by specifying the lowest pressure, the largest pressure, and the spacing (::math:`\Delta p`).
The temperature axis can be parameterized in a similar fashion.
The user also provides the value of the constant salinity used in the computations.

Then, during the simulation, Step 1 starts with a look-up in the pre-computed table to get the CO2 solubility as a function of pressure and temperature.

Step 2: Computation of the phase densities and phase viscosities
================================================================

Gas density and viscosity
-------------------------



Brine density and viscosity 
---------------------------



Methods Description
=================================

The methods described below can be found under in classes under *coreComponent/constitutive/fluid/PVTFunctions*.
They are called in ``MultiPhaseMultiComponentFluidUpdate::Compute`` kernel wrapper and
register as *PVTFunction* and *FlashModel* entry in Catalog.

Density and viscosity methods
-----------------------------
**Methods of CO2 density calculation**
The ``SpanWagnerCO2DensityFunction`` class holds table to compute density depending
on the pressure and the temperature as input data

 .. code-block:: c

    void SpanWagnerCO2DensityFunction::Evaluation( EvalVarArgs const & pressure,
                                                   EvalVarArgs const & temperature,
                                                   arraySlice1d< EvalVarArgs const > const & GEOSX_UNUSED_PARAM( phaseComposition ),
                                                   EvalVarArgs & value,
                                                   bool useMass ) const


``input parameters:``

* pressure - pressure (Pa),  :math:`$\text{P < 800\; MPa}$`
* temperature - temperature (C) :math:`$\text{ 200\;K \leq T \leq 1100\;K}$`
* phaseComposition - dimless mole or mass fraction (-) , unused in the calculation
* useMass - flag true if using a mass formulation (rather than molar)

``return parameter:`` 

* value - density (:math:`$\text{kg/m^3}$` or :math:`$\text{mol/m^3}$`), as in Span and Wagner (1996)





**Function of CO2 viscosity calculation**
The ``FenghourCO2ViscosityFunction`` class holds table to compute viscosity depending
on the pressure and the temperature as input data.

 .. code-block:: c

    virtual void FenghourCO2ViscosityFunction:Evaluation( EvalVarArgs const & pressure,
                                                          EvalVarArgs const & temperature,
                                                          arraySlice1d< EvalVarArgs const > const & GEOSX_UNUSED_PARAM( phaseComposition ),
                                                          EvalVarArgs & value,
                                                          bool GEOSX_UNUSED_PARAM( useMass )) const:


``input parameters:``

* pressure - pressure (Pa),  :math:`$\text{P < 300\; MPa}$`
* temperature - temperature (C)  :math:`$\text{200\;K \leq T \leq 1500\;K}$`
* phaseComposition - dimless mole or mass fraction (-) , unused in the calculation
* useMass - flag true if using a mass formulation (rather than molar)


``return parameter:`` 

* value - viscosity (Pa s) as in Fenghour and Wakeman (1998)

**Method of brine density calculation**
The brine *salinity* (molal unit), :math:`$\text{ m \leq 5\; \mbox{molal}}$`, is an important parameter of the brine model.
It is set and register in the object constructor. The evaluation is then interfaced as,

 .. code-block:: c

  void BrineCO2DensityFunction:Evaluation( EvalVarArgs const & pressure,
                                           EvalVarArgs const & temperature,
                                           arraySlice1d< EvalVarArgs const > const & phaseComposition,
                                           EvalVarArgs & value,
                                           bool useMass ) const


``input parameters:``

* P - pressure (Pa),  :math:`$\text{P \leq 50\;MPa}$`
* temperature - temperature (C), :math:`$\text{10\;C \leq T \leq 350\;C}$`
* phaseComposition - dimless mole or mass fraction (-)
* useMass - flag true if using a mass formulation (rather than molar)


``return parameter:`` 

* value - density (:math:`$\text{kg/m^3}$` or :math:`$\text{mol/m^3}$`) as in Eq(4) in Phillips et al. (1981).

**Method of brine viscosity calculation**

The brine viscosity is modeled as being affine function of the temperature with coefficients depending
on salinity :math:`$\text{m}$`,

 .. math::
    \mu_{brine} = \mu_{ref}\; ( a_0(m)\; +\; b_0(m)\,T) \;\;\;\;\;\;\mbox{w/}\;\;\;\; \mu_{ref} = \mu_{w}(298 K)

It is then interfaced as,

 .. code-block:: c

    void BrineViscosityFunction::Evaluation( EvalVarArgs const & GEOSX_UNUSED_PARAM( pressure ),
                                             EvalVarArgs const & temperature,
                                             arraySlice1d< EvalVarArgs const > const & GEOSX_UNUSED_PARAM( phaseComposition ),
                                             EvalVarArgs & value,
                                             bool GEOSX_UNUSED_PARAM( useMass )) const

``input parameters:``

* pressure - pressure (Pa), unused in the calculation
* temperature - temperature (C), :math:`$\text{10\;C \leq T \leq 350\;C}$`
* phaseComposition - dimless mole or mass fraction (-) , unused in the calculation
* useMass - flag true if using a mass formulation (rather than molar), unused in the calculation


``return parameter:``

* value - viscosity (Pa s), as in Eq.(1) in Phillips et al. (1981)


Flash methos
--------------
The following method of class ``CO2SolubilityFunction`` is used as flash to compute dissolution of gaseous Co2
into brine.

 .. code-block:: c

    void CO2SolubilityFunction::Partition( EvalVarArgs const & pressure,
                                           EvalVarArgs const & temperature,
                                           arraySlice1d< EvalVarArgs const > const & compFraction,
                                           arraySlice1d< EvalVarArgs > const & phaseFraction,
                                           arraySlice2d< EvalVarArgs > const & phaseCompFraction ) const

``input parameters:``

* pressure - pressure (Pa), unused in the calculation
* temperature - temperature (C), :math:`$\text{10\;C \leq T \leq 350\;C}$`
* compFraction - dimless mole or mass fraction (-)


``return parameter:``

* phaseFraction - dimless pore space fraction (-)
* phaseCompFraction - dimless mol of mass fraction for each phase (-)


Parameters
=========================

The model is represented by ``<MultiPhaseMultiComponentFluid>`` node in the input.

The following attributes are supported:

.. include:: ../../../coreComponents/fileIO/schema/docs/MultiPhaseMultiComponentFluid.rst

Supported phase names are:

======== ===========
Value     Comment
======== ===========
gas      Gas phase
water    Water phase
======== ===========

Supported component names are:

============= ===========
Value         Comment
============= ===========
co2,CO2       co2 component
water,liquid  Water component
============= ===========

Example
=========================

.. code-block:: xml

    <Constitutive>
        <MultiPhaseMultiComponentFluid
          name="fluid1"
          phaseNames="{ gas, water }"
          componentNames="{ co2, water }"
          componentMolarWeight="{ 44e-3, 18e-3 }"
          phasePVTParaFiles="{ pvtgas.txt, pvtliquid.txt }"
          flashModelParaFile="co2flash.txt"/>
    </Constitutive>




References
=================================

- A. Fenghour and W. A. Wakeman, `The viscosity of carbon dioxide
  <https://aip.scitation.org/doi/abs/10.1063/1.556013>`__, J. Phys. Chem. Ref.
  Data, vol. 27, pp. 31-44, 1998.

- S. L. Phillips et al., `A technical databook for geothermal energy
  utilization <https://escholarship.org/content/qt5wg167jq/qt5wg167jq.pdf>`__,
  Lawrence Berkeley Laboratory report, 1981.

- R. Span and W. Wagner., `A new equation of state for carbon dioxide covering
  the fluid region from the triple-point temperature to 1100 K at pressure up
  to 800 MPa <https://aip.scitation.org/doi/abs/10.1063/1.555991>`__, J. Phys.
  Chem. Ref. Data, vol. 25, pp. 1509-1596, 1996.
