Enthalpy
========

For non-isothermal simulations, such as thermal recovery processes or geothermal operations, the accurate calculation of phase enthalpies and internal energies is required to resolve the energy balance equation.

Compositional enthalpy model
----------------------------

The compositional enthalpy model for a multiphase fluid system computes the total enthalpy of the phase as the sum of an ideal (polynomial) enthalpy and a real fluid correction (or departure enthalpy), which is derived from the equation of state.

The total enthalpy of a fluid phase is defined as:

.. math::
    H_{total} = H_{ideal} + H_{dep}

Ideal gas contribution
~~~~~~~~~~~~~~~~~~~~~~

The ideal enthalpy represents the enthalpy of the mixture assuming it behaves as an ideal gas, relying on the heat capacity of the individual components. The specific heat capacity :math:`C_{p,i}` for each component :math:`i` is modeled as a 4th-order Poling polynomial (Poling et al., 2001) based on the temperature difference from a reference state, :math:`\Delta T = T - T_{ref}`:

.. math::
    C_{p,i}(T) = a_{i,0} + a_{i,1}\Delta T + a_{i,2}\Delta T^2 + a_{i,3}\Delta T^3 + a_{i,4}\Delta T^4

To find the ideal enthalpy for a component, this heat capacity is integrated with respect to temperature from the reference temperature to the current temperature, and then added to a baseline reference enthalpy :math:`H_{ref,i}`:

.. math::
    H_{ideal,i}(T) = H_{ref,i} + \Delta T \left( a_{i,0} + \frac{a_{i,1}}{2}\Delta T + \frac{a_{i,2}}{3}\Delta T^2 + \frac{a_{i,3}}{4}\Delta T^3 + \frac{a_{i,4}}{5}\Delta T^4 \right)

The total ideal enthalpy for the phase is simply the mole-fraction weighted sum of the individual component ideal enthalpies:

.. math::
    H_{ideal} = \sum_i x_i H_{ideal,i}(T)

Departure function
~~~~~~~~~~~~~~~~~~

The dimensionless enthalpy departure is given by:

.. math::
    \frac{H_{dep}}{RT} = Z - 1 + \frac{A + T\frac{\partial A}{\partial T}}{B(\delta_1 - \delta_2)} \ln \left( \frac{Z + \delta_1 B}{Z + \delta_2 B} \right)

where :math:`Z` is the compressibility factor of the mixture.

This model is recommended whenever non-isothermal compositional physics are required.

Mass formulation conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~

When the simulator is configured to use the mass formulation, the underlying inputs and thermodynamic calculations remain entirely in the molar formulation. The final phase molar enthalpy :math:`H_{total}` [J/mol] is converted to a mass-based phase enthalpy :math:`H_{mass}` [J/kg] by dividing by the phase molar weight :math:`M_{phase}`:

.. math::
    H_{mass} = \frac{H_{total}}{M_{phase}}

where the phase molar weight is the mole-fraction weighted sum of the component molar weights :math:`M_i`:

.. math::
    M_{phase} = \sum_i x_i M_i

Because all inputs must be provided in molar units, if a user's specific heat capacity coefficients or reference enthalpies are originally defined in a mass formulation (e.g., :math:`\hat{a}_{i,k}` in [J/(kg K^n)]), they must be scaled by the respective component's molar weight :math:`M_i` [kg/mol] prior to input to yield the molar coefficients :math:`a_{i,k}`:

.. math::
    a_{i,k} = \hat{a}_{i,k} M_i

Implementing the Michaelides Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The the Michaelides (1981) enthalpy model which is used for the CO2-brine model (:ref:`CO2-EOS`) can be implemented using the polynomial framework by setting specific values for the water component. 

The rational function :math:`Y(m)` below serves as a universal algebraic template to compute the required simulator inputs as a function of the salt molality :math:`m` [mol/kg]. Because the underlying thermodynamic derivations rely on fixed mixing rules and molecular weights, the calculation for every target parameter condenses into this single functional form:

.. math::
    Y(m) = \frac{p_0 + p_1 m + p_2 m^2 + p_3 m^3}{m + 55.5247}

To determine the input value for a specific simulator parameter, locate that parameter in the table below and substitute its corresponding row of coefficients (:math:`p_0` through :math:`p_3`) along with your desired molality :math:`m` into the :math:`Y(m)` equation. The evaluated result provides the direct numerical input required for the value.

The values in the table are hard-coded to a reference temperature of 293.15K and molar weights of 18.01 g/mol for water and 58.44 g/mol for salt.

=================== ============= ============= ============= =============
Parameter           :math:`p_0`   :math:`p_1`   :math:`p_2`   :math:`p_3`
=================== ============= ============= ============= =============
:math:`H_{ref,w}`   6.49471e+04   -2.97460e+04  -1.19193e+04  8.46747e+02
:math:`a_{w,0}`     4.64595e+03   4.69199e+02   2.29413e+02   -1.56721e+01
:math:`a_{w,1}`     -7.53304e+00  -5.72452e+00  -2.69521e+00  1.84321e-01
:math:`a_{w,2}`     3.73590e-02   1.20531e-02   9.02790e-03   -6.15000e-04
:math:`a_{w,3}`     0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00
:math:`a_{w,4}`     0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00
=================== ============= ============= ============= =============

.. note::
    Because the original Michaelides enthalpy model relies on a 3rd-order polynomial fit with respect to temperature (Michaelides, 1981), the higher-order simulator heat capacity coefficients (:math:`a_{w,3}` and :math:`a_{w,4}`) evaluate exactly to zero.

This formulation is bounded by specific thermodynamic assumptions: it achieves about ±3% accuracy relative to standard geothermal fluid tables, remains valid for saturation temperatures up to 320 °C (with possible divergence below 0 °C), and applies to salinity levels up to roughly 5-6 mol/kg using equivalent NaCl content. It further assumes fixed molar masses of 18.01 g/mol for water and 58.44 g/mol for NaCl, meaning any deviation without polynomial refitting will introduce scaling errors. Finally, the model is limited to liquid-phase (brine) enthalpy, and phase transitions such as flashing to vapor must be handled separately using appropriate gas mixing rules.

Parameters
~~~~~~~~~~

* ``enthalpyReferenceTemperature``: (Optional) A scalar defining the reference temperature :math:`T_{ref}` with a default of 298.15 K (25 C). Unit: [K].
* ``referenceEnthalpy``: (Optional) A 2D array mapping the reference enthalpy for each component within each phase at the reference temperature. The dimensions must match ``[number of phases, number of components]``. Defaults to zero. Unit: [J/mol].
* ``componentHeatCapacityCoefficients``: (Required) A 2D array containing the 5 polynomial coefficients (:math:`a_0` to :math:`a_4`) for the specific heat capacity of each component. The dimensions must strictly be ``[number of components, 5]``. To maintain dimensional consistency with a heat capacity in [J/(mol K)] and temperature difference in [K], the units for the specific coefficients are:

  * :math:`a_{i,0}`: [J/(mol K)]
  * :math:`a_{i,1}`: [J/(mol K^2)]
  * :math:`a_{i,2}`: [J/(mol K^3)]
  * :math:`a_{i,3}`: [J/(mol K^4)]
  * :math:`a_{i,4}`: [J/(mol K^5)]
