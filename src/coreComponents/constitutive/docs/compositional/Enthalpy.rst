========
Enthalpy
========

For non-isothermal simulations, such as thermal recovery processes or geothermal operations, the accurate calculation of phase enthalpies and internal energies is required to resolve the energy balance equation.

Compositional Enthalpy Model
----------------------------

The enthalpy of a compositional fluid phase is computed using a rigorous thermodynamic departure function approach. The true enthalpy at a given pressure and temperature, :math:`H(P, T)`, is expressed as the sum of an ideal gas enthalpy at the same temperature, :math:`H_{ideal}(T)`, and an EoS-derived departure enthalpy, :math:`H_{departure}(P, T, x)`:

.. math::
    H(P, T) = H_{ideal}(T) + H_{departure}(P, T, x)

**Ideal Gas Contribution**

The ideal gas enthalpy is derived by integrating a temperature-dependent specific heat capacity polynomial, :math:`C_p(T)`, from a defined reference temperature, :math:`T_{ref}`, to the current system temperature:

.. math::
    C_p(\Delta T) = a_0 + a_1 \Delta T + a_2 \Delta T^2 + a_3 \Delta T^3 + a_4 \Delta T^4

where :math:`\Delta T = T - T_{ref}`. The integrated polynomial is added to a user-supplied reference enthalpy for each component.

**Departure Function**

The departure enthalpy accounts for real-gas behaviour and molecular interactions. It is calculated exactly from the Equation of State by evaluating the partial derivative of the EoS attractive mixture parameter, :math:`a`, with respect to temperature. 

This model is recommended whenever non-isothermal compositional physics are required. 

* **Catalog Name:** ``CompositionalEnthalpy``
* **Parameters:**
  * ``enthalpyReferenceTemperature``: The reference temperature :math:`T_{ref}`.
  * ``referenceEnthalpy``: The enthalpy of each component at the reference temperature.
  * ``componentHeatCapacityCoefficients``: The polynomial coefficients (:math:`a_0` to :math:`a_4`) for the heat capacity of each component.