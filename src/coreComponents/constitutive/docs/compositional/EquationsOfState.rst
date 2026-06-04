Equations of state
------------------

An Equation of State (EoS) mathematically relates the pressure, volume, temperature, and composition of a fluid. In compositional modelling, the EoS serves a dual purpose: 
1. It calculates the fugacity coefficients, :math:`\phi_i`, required to establish thermodynamic phase equilibrium.
2. It calculates the compressibility factor, :math:`Z`, required to determine the phase density.

The fugacity coefficient of a component is derived from the exact thermodynamic relationship involving the integration of the EoS volume departure with respect to pressure. 

Cubic equation of state
~~~~~~~~~~~~~~~~~~~~~~~

The Peng-Robinson (PR) (Peng and Robinson, 1976; Robinson and Peng, 1978) and Soave-Redlich-Kwong (SRK) (Soave, 1972) equations of state are widely used for hydrocarbon systems encompassing natural gas, gas condensates, and volatile oils. The fundamental pressure-volume-temperature (PVT) relationships for these equations are defined as:

* Peng-Robinson (PR):
  
  .. math::
      P = \frac{RT}{v - b} - \frac{a(T)}{v(v + b) + b(v - b)}

* Soave-Redlich-Kwong (SRK):
  
  .. math::
      P = \frac{RT}{v - b} - \frac{a(T)}{v(v + b)}

In these PVT relationships, several intermediate parameters hold specific physical meanings:

* :math:`v` represents the molar volume of the fluid phase.
* :math:`R` is the universal gas constant, serving as the fundamental scaling factor relating macroscopic thermodynamic state variables to the molar energy.
* :math:`a(T)` is the temperature-dependent attractive parameter, representing the intermolecular attractive forces between the fluid molecules.
* :math:`b` is the repulsive co-volume parameter, representing the effective physical volume occupied by the fluid molecules themselves.

To efficiently solve these equations in a multicomponent mixture, they are rewritten into a generalized cubic equation for the compressibility factor, :math:`Z`:

.. math::
    Z^3 + E_2 Z^2 + E_1 Z + E_0 = 0

where the coefficients are defined by the dimensionless mixture parameters :math:`A` and :math:`B`:

.. math::
    E_2 = (\delta_1 + \delta_2 - 1)B - 1

.. math::
    E_1 = A + \delta_1 \delta_2 B^2 - (\delta_1 + \delta_2)B(B + 1)

.. math::
    E_0 = -[AB + \delta_1 \delta_2 B^2(B + 1)]

The constants :math:`\delta_1` and :math:`\delta_2` determine the specific EoS used:

* Peng-Robinson (PR): :math:`\delta_1 = 1 + \sqrt{2}`, :math:`\delta_2 = 1 - \sqrt{2}`
* Soave-Redlich-Kwong (SRK): :math:`\delta_1 = 0`, :math:`\delta_2 = 1`

The mixture parameters :math:`A` and :math:`B` are calculated using standard Van der Waals mixing rules over the phase mole fractions, :math:`x_i`:

.. math::
    A = \sum_i \sum_j x_i x_j (1 - k_{ij}) \sqrt{A_i A_j}

.. math::
    B = \sum_i x_i B_i

where :math:`k_{ij}` are the binary interaction coefficients (BICs). The pure component parameters :math:`A_i` and :math:`B_i` are defined as:

.. math::
    A_i = \Omega_a \frac{P_{r,i}}{T_{r,i}^2} \alpha_i

.. math::
    B_i = \Omega_b \frac{P_{r,i}}{T_{r,i}}

where :math:`P_{r,i} = P/P_{c,i}` and :math:`T_{r,i} = T/T_{c,i}` are the reduced pressure and temperature of component :math:`i`. 

The dimensionless constants :math:`\Omega_a` and :math:`\Omega_b` are derived mathematically by applying the critical point constraints (where the first and second pressure-volume derivatives vanish) to dictate the weight of the attractive and repulsive forces. The exact values implemented in the code are:

* Peng-Robinson: :math:`\Omega_a = 0.457235529`, :math:`\Omega_b = 0.077796074`
* Soave-Redlich-Kwong: :math:`\Omega_a = 0.42748`, :math:`\Omega_b = 0.08664`

The alpha function, :math:`\alpha_i`, accounts for the temperature dependence of the attractive parameter and relies on the acentric factor, :math:`\omega_i`, of each component. For both equations, it takes the form:

.. math::
    \alpha_i = \left( 1 + m_i \left( 1 - \sqrt{T_{r,i}} \right) \right)^2

The parameter :math:`m_i` is computed differently depending on the chosen equation of state:

* Peng-Robinson:
  If :math:`\omega_i < 0.49`:
  
  .. math::
      m_i = 0.37464 + 1.54226 \omega_i - 0.26992 \omega_i^2
      
  If :math:`\omega_i \ge 0.49`:
  
  .. math::
      m_i = 0.3796 + 1.485 \omega_i - 0.164423 \omega_i^2 + 0.016666 \omega_i^3

* Soave-Redlich-Kwong:
  
  .. math::
      m_i = 0.480 + 1.574 \omega_i - 0.176 \omega_i^2

Finally, the fugacity coefficient, :math:`\phi_i`, for each component is calculated using the derived roots for :math:`Z`:

.. math::
    \ln \phi_i = (Z - 1) \frac{B_i}{B} - F - G \left( 2 \sum_j x_j (1 - k_{ij}) \sqrt{A_i A_j} - A \frac{B_i}{B} \right) E

where the intermediate terms are defined as:

.. math::
    E = \ln \left( \frac{Z + \delta_1 B}{Z + \delta_2 B} \right), \quad F = \ln(Z - B), \quad G = \frac{1}{(\delta_1 - \delta_2)B}

Soreide-Whitson equation of state
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Soreide-Whitson model represents a specialized modification of the Peng-Robinson EoS tailored specifically for aqueous-hydrocarbon mixtures, such as those encountered in carbon sequestration (CO2-brine systems) or water-flooding in oil reservoirs. It adapts the standard cubic EoS by capturing the "salting-out" effect, where increased brine salinity reduces the solubility of hydrocarbons and CO2 in the aqueous phase. It is highly recommended whenever accurate modelling of gas dissolution in saline water is critical.

This model introduces two primary modifications to the standard Peng-Robinson formulation, applying special treatment exclusively to the water component and its interactions within the aqueous phase and follows closely the description given by Soreide and Whitson (1992):

First, the standard alpha function, :math:`\alpha_i`, is replaced specifically for the water component by a salinity-dependent formulation. If :math:`c_{sw}` represents the brine salinity, the water alpha function is evaluated as:

.. math::
    \alpha_w = \left[ 1 + 0.4530 \left( 1 - T_{r,w} \left( 1 - 0.0103 c_{sw}^{1.1} \right) \right) + 0.0034 \left( T_{r,w}^{-3} - 1 \right) \right]^2

where :math:`T_{r,w}` is the reduced temperature of water.

Second, unlike standard equations of state that use constant binary interaction coefficients (BICs), this model calculates dynamic aqueous phase BICs, :math:`k_{ij}^{AQ}`, exclusively for pairs where one component is water and the other is a hydrocarbon or specific non-hydrocarbon. 

For water-hydrocarbon pairs, the interaction coefficient is modeled as a polynomial function of the hydrocarbon's reduced temperature, :math:`T_{r,i}`, its acentric factor, :math:`\omega_i`, and the brine salinity:

.. math::
    k_{ij}^{AQ} = A_0 (1 + \alpha_0 c_{sw}) + A_1 T_{r,i} (1 + \alpha_1 c_{sw}) + A_2 T_{r,i}^2 (1 + \alpha_2 c_{sw})

The constants for this polynomial depend on the acentric factor. A special treatment is applied for methane, identified by an acentric factor :math:`\omega_i \le 0.02` (methane's acentric factor is approximately 0.0108). To achieve higher accuracy for methane-brine systems, a specifically tuned set of empirical constants is used instead of the generalized hydrocarbon constants. 

For standard hydrocarbons (:math:`\omega_i > 0.02`):

* :math:`A_0 = 1.1120 - 1.7369 \omega_i^{-0.1}`
* :math:`A_1 = 1.1001 + 0.8360 \omega_i`
* :math:`A_2 = -0.15742 - 1.0988 \omega_i`
* :math:`\alpha_0 = 0.017407`
* :math:`\alpha_1 = 0.033516`
* :math:`\alpha_2 = 0.011478`

For methane (:math:`\omega_i \le 0.02`):

* :math:`A_0 = 1.616705 - 1.884534 \omega_i^{-0.1}`
* :math:`A_1 = 0.8157014 + 0.8723632 \omega_i`
* :math:`A_2 = -0.0887821 - 0.8767864 \omega_i`
* :math:`\alpha_0 = 0.09988448`
* :math:`\alpha_1 = 0.1485516`
* :math:`\alpha_2 = 0.2111324`

For specific non-hydrocarbon components interacting with water, the dynamic BICs are calculated using explicitly tuned correlations based on the component type:

For carbon dioxide (CO2):

.. math::
    k_{ij}^{AQ} = A_0 a_0 + A_1 a_1 T_{r,i} + A_2 a_2

where :math:`A_0 = -0.31092`, :math:`A_1 = 0.23580`, :math:`A_2 = -21.2566`, :math:`a_0 = 1.0 + 0.15587 c_{sw}^{0.7505}`, :math:`a_1 = 1.0 + 0.17837 c_{sw}^{0.979}`, and :math:`a_2 = \exp(-6.7222 T_{r,i} - c_{sw})`.

For nitrogen (N2):

.. math::
    k_{ij}^{AQ} = A_0 a_0 + A_1 a_1 T_{r,i}

where :math:`A_0 = -1.70235`, :math:`A_1 = 0.44338`, :math:`a_0 = 1.0 + 0.025587 c_{sw}^{0.75}`, and :math:`a_1 = 1.0 + 0.08126 c_{sw}^{0.75}`.

For hydrogen sulfide (H2S):

.. math::
    k_{ij}^{AQ} = A_0 + A_1 T_{r,i}

where :math:`A_0 = -0.20441` and :math:`A_1 = 0.23426`.

For hydrogen (H2) the correlation due to Chabab et al. (2024) is used:

.. math::
    k_{ij}^{AQ} = D_0 (1 + a_0 c_{sw}) + D_1 T_{r,i} (1 + a_1 c_{sw}) + D_2 \exp(D_3 T_{r,i})

where :math:`D_0 = -2.11917`, :math:`D_1 = 0.14888`, :math:`D_2 = -13.01835`, :math:`D_3 = -0.43946`, :math:`a_0 = -0.0226322`, and :math:`a_1 = -0.0044736`.

Component naming requirements
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because the Soreide-Whitson EoS applies highly specific correlations based on the exact type of component interacting with water, it relies on the component's assigned name to identify its type. This name matching is case-insensitive.

To trigger the correct dynamic binary interaction coefficients, you must use the following exact names for the special components:

* Water: ``H2O``
* Carbon dioxide: ``CO2``
* Nitrogen: ``N2``
* Hydrogen sulfide: ``H2S``
* Hydrogen: ``H2``

If a component's name does not match any of these predefined strings (for example, if you name it "CarbonDioxide", "C1", or "Methane"), the model will automatically default to categorising it as a generic hydrocarbon (``hc``). While this is the intended and correct behavior for actual hydrocarbons (where methane is safely identified via its acentric factor), using an unrecognized name for a special non-hydrocarbon (like CO2) will cause the simulator to mistakenly apply the generic hydrocarbon-water correlation instead of its highly tuned specific correlation.

Model parameters
~~~~~~~~~~~~~~~~

Equations of state are assigned per-phase within the specific fluid model XML block (e.g., ``<CompositionalTwoPhaseFluidPhillipsBrine>``).

The ``equationsOfState`` attribute takes a list specifying the EoS type for each phase corresponding to the order defined in the ``phaseNames`` list. Standard component properties necessary for the EoS calculations must also be provided in matching component order.

Here is an example demonstrating how to specify the parameters for a two-phase system using the Soreide-Whitson EoS for the aqueous phase and Peng-Robinson for the gas phase:

.. code-block:: xml

    <CompositionalTwoPhaseFluidPhillipsBrine
      name="brine"
      phaseNames="{ liquid, gas }"
      equationsOfState="{ SoreideWhitson, PengRobinson }"
      componentNames="{ CH4, CO2, H2O }"
      componentCriticalPressure="{ 4.59920e+06, 7.37730e+06, 2.20640e+07 }"
      componentCriticalTemperature="{ 1.90564e+02, 3.04128e+02, 6.47096e+02 }"
      componentCriticalVolume="{ 9.86278e-05, 9.41185e-05, 5.59480e-05 }"
      componentAcentricFactor="{ 1.14200e-02, 2.23940e-01, 3.44300e-01 }"
      componentMolarWeight="{ 1.60425e-02, 4.40095e-02, 1.80153e-02 }"
      salinity="0.6" />

* ``equationsOfState``: List specifying the EoS type (e.g., ``PengRobinson``, ``SoaveRedlichKwong``, ``SoreideWhitson``) corresponding to each phase.
* ``componentCriticalPressure``: Critical pressure of each component [Pa].
* ``componentCriticalTemperature``: Critical temperature of each component [K].
* ``componentAcentricFactor``: Acentric factor of each component.
* ``componentVolumeShift``: Volume shift parameter for each component (optional).
* ``componentBinaryCoeff``: Flattened :math:`N \times N` matrix of binary interaction coefficients (optional).
* ``salinity``: Brine salinity [mol/kg] utilized by the Soreide-Whitson aqueous model (optional, defaults to 0.0).