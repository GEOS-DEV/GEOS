.. _CompressibleSinglePhaseFluid:

############################################
Compressible single phase fluid model
############################################

Overview
=========================

This model represents a compressible single-phase fluid with constant compressibility
and pressure-dependent viscosity.
For thermal simulations, fluid density and viscosity are temperature-dependent, governed by a constant thermal
expansion coefficient and a constant temperature-viscosity coefficient, respectively.
These assumptions are valid for slightly compressible fluids, such as water, and some
types of oil with negligible amounts of dissolved gas.

Specifically, fluid density is computed as

.. math::
   \rho(p) = \rho_0  e^{c_\rho(p - p_0)}

for isothermal cases and,

.. math::
   \rho(p,T) = \rho_0  e^{c_\rho(p - p_0)}  e^{-\beta_\rho (T - T_0)}

for thermal cases.

where :math:`c_\rho` is compressibility, :math:`p_0` is reference pressure, :math:`\rho_0` is the
density at the reference pressure, :math:`\beta_\rho` is the fluid thermal expansion coefficient,
:math:`T_0` is reference temperature.

Similarly,

.. math::
   \mu(p) = \mu_0 e^{c_\mu(p - p_0)}

for isothermal cases and,

.. math::
   \mu(p,T) = \mu_0 e^{c_\mu(p - p_0)} e^{-\beta_\mu (T - T_0)}

for thermal cases.

where :math:`c_\mu` is viscosibility (viscosity compressibility), :math:`\mu_0` is reference viscosity,
:math:`\beta_\rho` is the fluid viscosity thermal coefficient.

Parameters
=========================

For the isothermal case, the model is represented by ``<CompressibleSinglePhaseFluid>`` node in the input.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/CompressibleSinglePhaseFluid.rst

Example
=========================

.. code-block:: xml

  <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  referencePressure="2.125e6"
                                  referenceDensity="1000"
                                  compressibility="1e-19"
                                  referenceViscosity="0.001"
                                  viscosibility="0.0"/>
  </Constitutive>

For the thermal case, the model is represented by ``<ThermalCompressibleSinglePhaseFluid>`` node in the input.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/ThermalCompressibleSinglePhaseFluid.rst

Example
=========================

.. code-block:: xml

  <Constitutive>
    <ThermalCompressibleSinglePhaseFluid  name="fluid"
                                          defaultDensity="1000"
                                          defaultViscosity="0.001"
                                          referencePressure="0.0"
                                          referenceTemperature="0"
                                          compressibility="5e-10"
                                          thermalExpansionCoeff="3e-4"
                                          viscosibility="0.0"
                                          viscosityExpansivity="6.21e-4"
                                          specificHeatCapacity="1"
                                          referenceInternalEnergy="0.99" />
  </Constitutive>
