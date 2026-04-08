.. _CO2-EOS:

##################################################################################
CO2-brine model 
##################################################################################


Summary
=======

The CO2-brine model implemented in GEOS includes two components (CO2 and H2O) that are transported by one or two fluid phases (the brine phase and the CO2 phase).
We refer to the brine phase with the subscript :math:`\ell` and to the CO2 phase with the subscript :math:`g` (although the CO2 phase can be in supercritical, liquid, or gas state).
Both the water component and the CO2 component can be present in the CO2 phase as well as in the brine phase depending on the solubility model selected.
Thus, considering the molar phase component fractions, :math:`y_{c,p}` (i.e., the fraction of the molar mass of phase :math:`p` represented by component :math:`c`) the following partition matrix determines the component distribution within the two phases:

.. math::
    \begin{bmatrix}
    y_{H2O,\ell} & y_{CO2,\ell} \\
    y_{H2O,g}    & y_{CO2,g}    \\
    \end{bmatrix}

The update of the fluid properties is done in two steps:

1. The phase fractions (:math:`\nu_p`) and phase component fractions (:math:`y_{c,p}`) are computed as a function of pressure (:math:`p`), temperature (:math:`T`), overall component fractions (:math:`z_c`), and a constant salinity.

2. Computation of the phase properties as functions of pressure, temperature, the updated phase component fractions, at a constant salinity.

  a. Phase densities (:math:`\rho_p`) 
  b. Phase viscosities (:math:`\mu_p`) 
  c. Phase enthalpies (:math:`H_p`) for thermal models.

Once the phase fractions, phase component fractions, phase densities, phase viscosities, phase enthalpies--and their derivatives with respect to pressure, temperature, and component fractions--have been computed, the :ref:`CompositionalMultiphaseFlow` proceeds to the assembly of the accumulation and flux terms.

.. note::
   The fluid model will always calculate the derivatives with respect to temperature irrespective of whether the flow solver is thermal or not. For isothermal simulations, the derivatives are ignored.

The models that are used in steps 1) and 2) are reviewed in more details below.

Step 1: Computation of the phase fractions and phase component fractions (flash)
================================================================================

At initialization, GEOS performs a preprocessing step to construct a two-dimensional table storing the values of CO2 solubility in brine and water solubility in the gas phase as functions of pressure, temperature, and a constant salinity. Solubility is calculated and provided as moles of the component dissolved or vaporized per mass of the phase. For CO2 this will be moles of CO2 dissolved per kg of water and for water this will be moles of H2O vapourised in the CO2 phase per kg of CO2.

There are 3 alternatives to providing solubility data for the phase partition calculation.

Duan and Sun (2003) model
-------------------------
When the solubility model of Duan and Sun (2003) is selected the solubility of CO2 in brine is calculated while the vaporization of H2O into the gas phase is ignored (zero solubility of H2O in gas). Solubility values are tabulated for the selected pressure and temperature points using the Duan and Sun (2003) model calculates CO2 solubility by combining a thermodynamic framework for the liquid phase with an accurate equation of state for the vapor phase, enabling precise predictions across wide ranges of temperature, pressure, and salinity.

Specifically, we solve the following nonlinear CO2 equation of state (equation (A1) in Duan and Sun, 2003) for each pair :math:`(p,T)` to obtain the reduced volume, :math:`V_r`.

.. math::
   \frac{p_r V_r}{T_r} &= 1 + \frac{a_1 + a_2/T^2_r + a_3/T^3_r}{V_r} 
   + \frac{a_4 + a_5/T^2_r + a_6/T^3_r}{V^2_r} + \frac{a_7 + a_8/T^2_r + a_9/T^3_r}{V^4_r} \\
   &+ \frac{a_{10} + a_{11}/T^2_r + a_{12}/T^3_r}{V^5_r} 
   + \frac{a_{13}}{T^3_r V^2_r} \big( a_{14} + \frac{a_{15}}{V^2_r} \big) \exp( - \frac{a_{15}}{V^2_r} )

where :math:`p_r = p / p_{crit}` and :math:`T_r = T / T_{crit}` are respectively the reduced pressure and the reduced temperature.
We refer the reader to Table (A1) in Duan and Sun (2003) for the definition of the coefficients :math:`a_i` involved in the previous equation. 
Using the reduced volume, :math:`V_r`, we compute the fugacity coefficient of CO2, :math:`\ln_{\phi}(p,T)`, using equation (A6) of Duan and Sun (2003).
To conclude this preprocessing step, we use the fugacity coefficient of CO2 to compute and store the solubility of CO2 in brine, :math:`s_{CO2}`, using equation (6) of Duan and Sun (2003):

.. math::
   \ln \frac{ y_{CO2} P }{ s_{CO2} } = \frac{\Phi_{CO2}}{RT} - \ln_{\phi}(p,T) + \sum_c 2 \lambda_c m + \sum_a 2 \lambda_a m + \sum_{a,c} \zeta_{a,c} m^2

where :math:`\Phi_{CO2}` is the chemical potential of the CO2 component, :math:`R` is the gas constant, and :math:`m` is the salinity.
The mole fraction of CO2 in the vapor phase, :math:`y_{CO2}`, is computed with equation (4) of Duan and Sun (2003).
Note that the first, third, fourth, and fifth terms in the equation written above are approximated using equation (7) of Duan and Sun (2003) as recommended by the authors.

Spycher and Pruess (2005) model
-------------------------------
The partitioning behavior of CO2 between aqueous and gas phases in chloride brines is computed using the model developed by Spycher and Pruess (2005), which builds upon earlier work by Spycher et al. (2003). This model accounts for both the solubility of CO2 in brine and the vaporization of H2O into the gas phase, offering a comprehensive thermodynamic treatment.

The calculation begins with the determination of equilibrium constants :math:`k_{CO2}` and :math:`k_{H2O}` at a reference pressure :math:`P_0 = 1 \, \text{bar}`. These constants are expressed as temperature-dependent correlations:

.. math::

   \log_{10} k_{CO2}(P_0, T) = a_{CO2} + b_{CO2} T + c_{CO2} T^2 + d_{CO2} T^3

.. math::

   \log_{10} k_{H2O}(P_0, T) = a_{H2O} + b_{H2O} T + c_{H2O} T^2 + d_{H2O} T^3

The coefficients :math:`a, b, c, d` are provided in Table 2 of Spycher et al. (2003). These equilibrium constants are then corrected for pressure using:

.. math::

   k_c(P, T) = k_c(P_0, T) \exp\left( \frac{(P - P_0) \overline{V_c}}{RT} \right)

where :math:`\overline{V_c}` is the average partial molar volume of the pure condensed component :math:`c`, representing either CO2 or H2O.

A fitted Redlich-Kwong equation of state is employed to compute the fugacity coefficients :math:`\Phi_{CO2}` and :math:`\Phi_{H2O}`, with parameters provided in Table 1 of Spycher et al. (2003).

By equating the fugacities derived from the reaction equilibrium and the definitions of fugacity and partial pressure, the following relationships are obtained:

.. math::

   \Phi_{H2O} y_{H2O,g} P = k_{H2O}(P, T) a_{H2O,\ell}

.. math::

   \Phi_{CO2} y_{CO2,g} P = k_{CO2}(P, T) a_{CO2,\ell}

where :math:`a` denotes the activity of the component in the liquid phase.

For systems with low solubility, activities can be approximated by mole fractions. This leads to the following expressions, as given by equations (8) and (10) of Spycher et al. (2003):

.. math::

   y_{H2O,g} = \frac{K_{H2O}(P_0, T) (1 - x_{CO2,\ell})}{\Phi_{H2O} P} 
   \exp\left( \frac{(P - P_0) \bar{V}_{H2O}}{RT} \right)

.. math::

   x_{CO2,\ell} = \frac{\Phi_{CO2}(1 - y_{H2O,g}) P}{55.508\, K_{CO2}(P_0, T)} 
   \exp\left( -\frac{(P - P_0) \bar{V}_{CO2}}{RT} \right)

These expressions are used to compute the mole fractions :math:`y_{H2O,g}` and :math:`x_{CO2,\ell}` in the gas and liquid phases, respectively. These are then converted to solubilities as

.. math::

   s_{CO2} = \frac{x_{CO2,\ell}}{1-x_{CO2,\ell}}M_{H2O}

.. math::

   s_{H2O} = \frac{y_{H2O,g}}{1-y_{H2O,g}}M_{CO2}

where :math:`M_{H2O}` is the molecular weight of water and :math:`M_{CO2}` is the molecular weight of CO2. 

Tabulated solubilities
----------------------
Alternatively, the user may provide solubility tables generated externally to GEOS. These must be supplied as two‑dimensional lookup tables (TableFunction) with pressure and temperature as the coordinate axes. A CO2 solubility table is required, while the water vaporization table is optional. If a water vaporization table is not provided, water vaporization is neglected (i.e., zero vaporization is assumed).

For all tables, pressure must be specified in pascals (Pa), temperature in kelvin (K), and solubility as moles of solute per unit mass of solvent (mol/kg).

Flash calculation
-----------------
The phase split calculation determines the equilibrium compositions of the two components in each of the two pgases of the fluid system. Given the overall mixture mole fractions :math:`z_{CO2}` and :math:`z_{H2O}`, the calculation isolates the liquid (brine) phase mole fractions, :math:`y_{CO2,\ell}` and :math:`y_{H2O,\ell}`, and the gas phase mole fractions, :math:`y_{CO2,g}` and :math:`y_{H2O,g}`. This flash formulation rigorously accounts for the mutual solubility of both components without assuming that either phase is predominantly a single component.

The given mutual solubility data, :math:`s_{CO2}` and :math:`s_{H2O}`, are initially provided in units of moles of solute per kilogram of solvent. The first step is to convert these values into molar ratios (moles of solute per mole of solvent) by multiplying them by the molar mass of the respective solvents, :math:`M_{H2O}` and :math:`M_{CO2}` (in kg/mol). This conversion yields the mole-to-mole solubilities :math:`S_{CO2}` and :math:`S_{H2O}`:

.. math::

   S_{CO2} &= s_{CO2} \times M_{H2O} \\
   S_{H2O} &= s_{H2O} \times M_{CO2}

By considering a total system size of 1 mole, the overall mole fractions :math:`z_{CO2}` and :math:`z_{H2O}` are equivalent to the total moles of each component. Let :math:`n_{CO2,\ell}` and :math:`n_{H2O,\ell}` represent the moles of CO2 and H2O in the liquid phase, and :math:`n_{CO2,g}` and :math:`n_{H2O,g}` represent the moles in the gas phase. The overall material balances are:

.. math::

   n_{CO2,\ell} + n_{CO2,g} &= z_{CO2} \\
   n_{H2O,\ell} + n_{H2O,g} &= z_{H2O}

The phase equilibrium is dictated by the molar solubility ratios:

.. math::

   n_{CO2,\ell} &= S_{CO2} \cdot n_{H2O,\ell} \\
   n_{H2O,g} &= S_{H2O} \cdot n_{CO2,g}

Substituting the material balances into the equilibrium relations produces a 2x2 linear system for the unknown moles of the solutes in their secondary phases, :math:`n_{CO2,\ell}` and :math:`n_{H2O,g}`:

.. math::

   n_{CO2,\ell} - S_{CO2} (z_{H2O} - n_{H2O,g}) &= 0 \\
   n_{H2O,g} - S_{H2O} (z_{CO2} - n_{CO2,\ell}) &= 0

Rearranging these equations reveals the standard 2x2 formulation for the flash calculation, which can be expressed in matrix form as:

.. math::

   \begin{bmatrix}
   1 & S_{CO2} \\
   S_{H2O} & 1
   \end{bmatrix}
   \begin{bmatrix}
   n_{CO2,\ell} \\
   n_{H2O,g}
   \end{bmatrix}
   =
   \begin{bmatrix}
   S_{CO2} \cdot z_{H2O} \\
   S_{H2O} \cdot z_{CO2}
   \end{bmatrix}

Solving this 2x2 system provides analytical expressions for :math:`n_{CO2,\ell}` and :math:`n_{H2O,g}`. The determinant of this system is :math:`1 - S_{CO2} S_{H2O}`. The solutions are:

.. math::

   n_{CO2,\ell} &= \frac{S_{CO2} (z_{H2O} - S_{H2O} z_{CO2})}{1 - S_{CO2} S_{H2O}} \\
   n_{H2O,g} &= S_{H2O} (z_{CO2} - n_{CO2,\ell})

With :math:`n_{CO2,\ell}` and :math:`n_{H2O,g}` calculated, the remaining moles are immediately found using the material balance equations:

.. math::

   n_{H2O,\ell} &= z_{H2O} - n_{H2O,g} \\
   n_{CO2,g} &= z_{CO2} - n_{CO2,\ell}

Finally, the phase-specific mole fractions are determined by normalizing the moles in each phase by the total moles of that specific phase:

.. math::

   y_{CO2,\ell} &= \frac{n_{CO2,\ell}}{n_{CO2,\ell} + n_{H2O,\ell}} \\
   y_{H2O,\ell} &= \frac{n_{H2O,\ell}}{n_{CO2,\ell} + n_{H2O,\ell}} \\
   y_{CO2,g} &= \frac{n_{CO2,g}}{n_{CO2,g} + n_{H2O,g}} \\
   y_{H2O,g} &= \frac{n_{H2O,g}}{n_{CO2,g} + n_{H2O,g}}

The overall phase fractions for the liquid and gas phases, :math:`\nu_\ell` and :math:`\nu_g`, can be calculated by summing the moles of the individual components within each respective phase. Because the formulation assumes a total system size of 1 mole, the phase fractions are directly equal to these sums: :math:`\nu_\ell = n_{CO2,\ell} + n_{H2O,\ell}` and :math:`\nu_g = n_{CO2,g} + n_{H2O,g}`.

.. note::
   If any calculated mole fraction (or intermediate molar quantity) evaluates to a negative number or falls below the numerical minimum threshold, it implies that the corresponding component has been entirely depleted from that phase. When this occurs, the algorithm overrides the calculation and clamps the composition for that component to zero, effectively denoting the physical absence of that phase component.

.. note::
   If the determinant :math:`1 - S_{CO2} S_{H2O}` approaches zero or falls below a minimum threshold for division, the 2x2 system is considered not soluble (ill-conditioned or non-physical). In this scenario, the flash calculation cannot proceed safely, and the algorithm will throw an input error indicating a failure to calculate solubility at the given pressure and temperature.
   
Step 2: Computation of the phase densities, viscosities and enthalpies
======================================================================

CO2 phase density, viscosity and enthalpy
-----------------------------------------

In GEOS, the computation of the CO2 phase density and viscosity is entirely based on look-up in precomputed tables.
The user defines the pressure (in Pascal) and temperature (in Kelvin) axis of the density table to be internally generated.

A correlation due to Span and Wagner (1996) is used to populate the table of CO2 phase densities and is valid for pressures less than :math:`8 \times 10^8` Pascal and temperatures less than 1073.15 Kelvin. Using the user defined parameters, GEOS internally constructs a two-dimensional table storing the values of density as a function of pressure and temperature.
As explained in the work of Span and Wagner (1996), the density is calculated by solving the following nonlinear Helmholtz energy equation for each pair :math:`(p,T)` to obtain the value of density, :math:`\rho_{g}`:

.. math::
   \frac{p}{RT\rho_{g}} = 1 + \delta \phi^r_{\delta}( \delta, \tau )

where :math:`R` is the gas constant, :math:`\delta := \rho_{g} / \rho_{crit}` is the reduced CO2 phase density, and :math:`\tau := T_{crit} / T` is the inverse of the reduced temperature.
The definition of the residual part of the energy equation, denoted by :math:`\phi^r_{\delta}`, can be found in equation (6.5), page 1544 of Span and Wagner (1996).
The coefficients involved in the computation of :math:`\phi^r_{\delta}` are listed in Table (31), page 1544 of Span and Wagner (1996).   
These calculations are done in a preprocessing step.

Similarly, a correlation due to Fenghour and Wakeham (1998) is used to populate the CO2 phase viscosity table and this correlation is valid for pressures less than :math:`3 \times 10^8` Pascal and temperatures less than 1493.15 Kelvin.  
This table is populated as explained in the work of Fenghour and Wakeham (1998) by computing the CO2 phase viscosity, :math:`\mu_g`, as follows:

.. math::
   \mu_{g} = \mu_{0}(T) + \mu_{excess}( \rho_{g}, T ) + \mu_{crit}( \rho_{g}, T )  
   
The "zero-density limit" viscosity, :math:`\mu_{0}(T)`, is computed as a function of temperature using equations (3), (4), and (5), as well as Table (1) of Fenghour and Wakeham (1998).
The excess viscosity, :math:`\mu_{excess}( \rho_{g}, T )`, is computed as a function of temperature and CO2 phase density (computed as explained above) using equation (8) and Table (3) of Fenghour and Wakeham (1998).
We currently neglect the critical viscosity, :math:`\mu_{crit}`.
These calculations are done in a preprocessing step.

For thermal simulations, the CO2 phase enthalpy is calculated using the Span and Wagner (1996) equation of state, which is based on a dimensionless Helmholtz free energy formulation. Again this is calculated in a preprocessing step and a table is precomputed for lookup.

The total dimensionless Helmholtz free energy, :math:`\phi`, is the sum of the ideal gas contribution (:math:`\phi^0`) and the residual contribution (:math:`\phi^r`):

.. math::

   \phi(\tau, \delta) = \phi^0(\tau, \delta) + \phi^r(\tau, \delta)

The independent state variables are evaluated as dimensionless, reduced parameters:

.. math::

   \delta = \frac{\rho}{\rho_c} \quad \text{and} \quad \tau = \frac{T_c}{T}

where:

* :math:`\rho` is the phase density.
* :math:`\rho_c` is the critical density of CO2.
* :math:`T` is the absolute temperature.
* :math:`T_c` is the critical temperature of CO2.

The specific enthalpy, :math:`h`, is derived from the partial derivatives of the Helmholtz free energy using the following thermodynamic relation:

.. math::

   h(T, \rho) = R T \left[ 1 + \tau \left( \frac{\partial \phi^0}{\partial \tau} + \frac{\partial \phi^r}{\partial \tau} \right) + \delta \frac{\partial \phi^r}{\partial \delta} \right]

where:

* :math:`R` is the specific gas constant for CO2.
* :math:`\frac{\partial \phi^0}{\partial \tau}` is the partial derivative of the ideal gas contribution with respect to :math:`\tau`.
* :math:`\frac{\partial \phi^r}{\partial \tau}` is the partial derivative of the residual contribution with respect to :math:`\tau`.
* :math:`\frac{\partial \phi^r}{\partial \delta}` is the partial derivative of the residual contribution with respect to :math:`\delta`.

During the simulation, the update of CO2 phase density, viscosity and enthalpy is simply done with a look-up in the precomputed tables. 

Brine density and viscosity using Phillips correlation
-------------------------------------------------------

The computation of the brine density involves a tabulated correlation presented in Phillips et al. (1981). 
The user specifies the (constant) salinity and defines the pressure and temperature axis of the brine density table in the form:

The pressure must be in Pascal and must be less than :math:`5 \times 10^7` Pascal.
The temperature must be in Kelvin and must be between 283.15 and 623.15 Kelvin.
Using these parameters, GEOS performs a preprocessing step to construct a two-dimensional table storing the brine density, :math:`\rho_{\ell,table}` for the specified salinity as a function of pressure and temperature using the expression:

.. math::
 
   \rho_{\ell,table} &= A + B x + C x^2 + D x^3 \\
   x &= c_1 \exp( a_1 m ) + c_2 \exp( a_2 T ) + c_3 \exp( a_3 P )

We refer the reader to Phillips et al. (1981), equations (4) and (5), pages 14 and 15 for the definition of the coefficients involved in the previous equation.
This concludes the preprocessing step.

Then, during the simulation, the brine density update proceeds in two steps.
First, a table look-up is performed to retrieve the value of density, :math:`\rho_{\ell,table}`.
Then, in a second step, the density is modified using the method of Garcia (2001) to account for the presence of CO2 dissolved in brine as follows:

.. math::

   \rho_{\ell} = \rho_{\ell,table} + M_{CO2} c_{CO2} - c_{CO2} \rho_{\ell,table} V_{\phi}

where :math:`M_{CO2}` is the molecular weight of CO2, :math:`c_{CO2}` is the concentration of CO2 in brine, and :math:`V_{\phi}` is the apparent molar volume of dissolved CO2.
The CO2 concentration in brine is obtained as:

.. math::

   c_{CO2} = \frac{y_{CO2,\ell} \rho_{\ell,table}}{M_{H2O}(1-y_{CO2,\ell})} 

where :math:`M_{H2O}` is the molecular weight of water. 
The apparent molar volume of dissolved CO2 is computed as a function of temperature using the expression:

.. math::

   V_{\phi} = 37.51 - 9.585 \times 10^{-2} T + 8.740 \times 10^{-4} T^2 - 5.044 \times 10^{-7} T^3

The brine viscosity is controlled by a salinity parameter provided by the user.
During the simulation, the brine viscosity is updated as a function of temperature using the analytical relationship of Phillips et al. (1981):

.. math::
   \mu_{\ell} = a T + b

where the coefficients :math:`a` and :math:`b` are defined as:

.. math::
   a &= \mu_{w}(T) \times 0.000629 (1.0 - \exp( -0.7 m ) ) \\
   b &= \mu_{w}(T) (1.0 + 0.0816 m + 0.0122 m^2 + 0.000128 m^3) 
   
where :math:`\mu_{w}` is the pure water viscosity computed as a function of temperature,
and :math:`m` is the user-defined salinity (in moles of NaCl per kg of brine).


Brine density and viscosity using Ezrokhi correlation
-------------------------------------------------------

Brine density :math:`\rho_l` is computed from pure water density :math:`\rho_w` at specified pressure and temperature corrected by Ezrokhi correlation presented in Zaytsev and Aseyev (1993):

.. math::
   log_{10}(\rho_l) &= log_{10}(\rho_w(P, T)) + A(T) x_{CO2,\ell} \\
   A(T) &= a_0 + a_1T +  a_2T^2,

where :math:`a_0, a_1, a_2` are correlation coefficients defined by the user, while :math:`x_{CO2,\ell}` is the mass fraction of the CO2 component in brine, computed from molar fractions as

.. math::
   x_{CO2,\ell} = \frac{M_{CO2}y_{CO2,\ell}}{M_{CO2}y_{CO2,\ell} + M_{H2O}y_{H2O,\ell}},

Pure water density is computed according to:

.. math::
   \rho_w = \rho_{w,sat}(T) e^{c_w * (P-P_{w,sat}(T))},

where :math:`c_w` is water compressibility provided by the user with a default value of :math:`4.5 \times 10^{-10} Pa^{-1}`, while :math:`\rho_{w,sat}(T)` and :math:`P_{w,sat}(T)` are density and pressure of saturated water at a given temperature.
Both are obtained through internally constructed tables tabulated as functions of temperature and filled with the steam table data from Engineering ToolBox (2003, 2004).

Brine viscosity :math:`\mu_{\ell}` is computed from pure water viscosity :math:`\mu_w` similarly:

.. math::
   log_{10}(\mu_l) &= log_{10}(\mu_w(P, T)) + B(T) x_{CO2,\ell} \\
   B(T) &= b_0 + b_1T +  b_2T^2,

where :math:`b_0, b_1, b_2` are correlation coefficients defined by the user.

Mass fraction of CO2 component in brine :math:`x_{CO2,\ell}` is exactly as in density calculation. The dependency of pure water viscosity from pressure is ignored, and it is approximated as saturated pure water viscosity:

.. math::
   \mu_w(P, T) = \mu_{w, sat} (T),

which is tabulated using internal table as a function of temperature based on steam table data Engineering ToolBox (2004).

Enthalpy of brine
-----------------
For thermal simulations, the brine phase enthalpy is calculated using the thermodynamic correlations developed by Michaelides (1981). To optimize computational performance during the simulation, these baseline correlations are evaluated in a preprocessing step to generate an enthalpy lookup table.

The enthalpy of a pure brine solution is fundamentally determined by the molality of the salt (typically modeled as NaCl) and the temperature of the system. The pure brine enthalpy is formulated as a mass-weighted sum of the constituent enthalpies alongside an enthalpy of mixing:

.. math::

   H_{\text{brine}} = x_1 h_1(T) + x_2 h_2(T) + m \Delta h(T, m)

where :math:`x_1` and :math:`x_2` are the mass fractions of water and salt, respectively, which are derived directly from the salt molality :math:`m`. The terms :math:`h_1(T)` and :math:`h_2(T)` represent the pure component enthalpies of water and salt, calculated as third-degree polynomials with respect to temperature :math:`T`. The enthalpy of mixing, :math:`\Delta h(T, m)`, is computed using a bivariate polynomial expansion that sums over powers of temperature and molality, scaled by an appropriate unit conversion factor.

To account for the presence of dissolved non-condensible gases such as carbon dioxide, the phase enthalpy is adjusted using a thermodynamic mixing rule. During the simulation, the precomputed pure brine enthalpy and the pure :math:`\text{CO}_2` enthalpy are interpolated from the lookup tables at the current pressure and temperature. Assuming a binary mixture of brine and :math:`\text{CO}_2`, the total phase enthalpy is calculated as a composition-weighted average:

.. math::

   H_{\text{phase}} = (1 - C) H_{\text{CO}_2} + C H_{\text{brine}}

where :math:`C` is the composition fraction of the brine component. Depending on the solver configuration, this combination is performed either on a mass basis or a molar basis. If a molar basis is used, the constituent enthalpies are internally scaled by the inverse molecular weights of their respective components prior to averaging.

Parameters
=========================

The models are represented by in the input ``<CO2BrinePhillipsFluid>``, ``<CO2BrineEzrokhiFluid>`` nodes for the isothermal model and ``<CO2BrinePhillipsThermalFluid>`` and ``<CO2BrineEzrokhiThermalFluid>`` nodes for their thermal counterparts.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/CO2BrinePhillipsFluid.rst

Supported phase names are:

======== ===========
Value     Comment
======== ===========
gas      CO2 phase
water    Water phase
======== ===========

Supported component names are:

============= ===============
Value         Component
============= ===============
co2,CO2       CO2 component
water,liquid  Water component
============= ===============

The salinity is specified in moles per kg of brine using the ``salinity`` attribute.

These models utilize pre-computed lookup tables for fluid properties. These tables are constructed over defined two-dimensional grids of pressure and temperature. The user must specify the discretization of these grids using specific XML attributes.

The grid definitions rely on the following XML attributes:

* ``pressureCoordinates``: A list of pressure coordinate boundaries or explicit values.
* ``pressureInterval``: (Optional) The step size for uniform pressure discretization.
* ``temperatureCoordinates``: A list of temperature coordinate boundaries or explicit values.
* ``temperatureInterval``: (Optional) The step size for uniform temperature discretization.

The sequence of evaluation points for either pressure or temperature is generated based on the combination of the ``*Coordinates`` array and the corresponding ``*Interval`` value.

Let :math:`C = \{c_1, c_2, \dots, c_n\}` be the sequence of values provided in the ``*Coordinates`` attribute, and :math:`\Delta` be the value provided in the ``*Interval`` attribute.

There are two distinct methods for specifying the grid points:

If an interval :math:`\Delta` is specified and is strictly positive (:math:`\Delta > 0`), the grid is generated as a linearly spaced sequence.

* Only the first value (:math:`c_1`) and the last value (:math:`c_n`) of the ``*Coordinates`` array are utilized as the lower and upper bounds. 
* Any intermediate values provided in the array are **ignored**.
* The generated sequence of points :math:`x_i` is defined as:

  .. math::

     x_i = c_1 + i \cdot \Delta \quad \text{for } i = 0, 1, \dots, k

  where :math:`k` is the maximum integer such that :math:`x_k \le c_n`.

If the ``*Interval`` attribute is omitted, set to zero, or is negative (:math:`\Delta \le 0`), the grid points are taken exactly as provided in the coordinates array.

* The generated grid points exactly map to the provided array: :math:`X = C = \{c_1, c_2, \dots, c_n\}`.

The same table discretization is used for all the properties.

When utilizing the ``CO2BrineEzrokhiFluid`` model or the ``CO2BrineEzrokhiThermalFluid`` mode, the following additional attributes are required to calculate the density and viscosity corrections using the Ezrokhi correlation:

* ``ezrokhiDensityCoefficients``: A three-element array :math:`\{A_1, A_2, A_3\}` defining the empirical coefficients used in the Ezrokhi equation of state to correct brine density.
* ``ezrokhiViscosityCoefficients``: A three-element array :math:`\{B_1, B_2, B_3\}` defining the empirical coefficients used to correct brine viscosity. These are optional and if not specified, zero values are used.

The default water compressibility can be changed using the ``waterCompressibility`` attribute.

Example
=======

.. code-block:: xml

    <Constitutive>
        <CO2BrinePhillipsFluid
            name="fluid"
            phaseNames="{ gas, water }"
            componentNames="{ co2, water }"
            componentMolarWeight="{ 44e-3, 18e-3 }"
            pressureCoordinates="{1.000e+05, 2.512e+05, 6.310e+05, 1.585e+06, 3.981e+06, 1.000e+07}"
            temperatureCoordinates="{283.15, 305.15, 327.15, 349.15, 371.15, 393.15}"
            salinity="3.0" />
    </Constitutive>


.. code-block:: xml

    <Constitutive>
        <CO2BrineEzrokhiFluid
            name="fluid"
            phaseNames="{ gas, water }"
            componentNames="{ co2, water }"
            componentMolarWeight="{ 44e-3, 18e-3 }"
            pressureCoordinates="{1.0e5, 6e7}"
            pressureInterval="1e5"
            temperatureCoordinates="{283.15, 393.5}"
            temperatureInterval="5"
            salinity="1.901285269"
            ezrokhiDensityCoefficients="{0.1033, -2.2991e-5, -2.3658e-6}"
            ezrokhiViscosityCoefficients="{0, 0, 0}" />
    </Constitutive>
    
References
==========

- Z. Duan and R. Sun, `An improved model calculating CO2 solubility in pure
  water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar.
  <https://doi.org/10.1016/S0009-2541(02)00263-2>`__, Chemical Geology,
  vol. 193.3-4, pp. 257-271, 2003.

- R. Span and W. Wagner, `A new equation of state for carbon dioxide covering
  the fluid region from the triple-point temperature to 1100 K at pressure up
  to 800 MPa <https://aip.scitation.org/doi/abs/10.1063/1.555991>`__, J. Phys.
  Chem. Ref. Data, vol. 25, pp. 1509-1596, 1996.

- A. Fenghour and W. A. Wakeham, `The viscosity of carbon dioxide
  <https://aip.scitation.org/doi/abs/10.1063/1.556013>`__, J. Phys. Chem. Ref.
  Data, vol. 27, pp. 31-44, 1998.

- S. L. Phillips et al., `A technical databook for geothermal energy
  utilization <https://escholarship.org/content/qt5wg167jq/qt5wg167jq.pdf>`__,
  Lawrence Berkeley Laboratory report, 1981.

- J. E. Garcia, Density of aqueous solutions of CO2. No. LBNL-49023.
  Lawrence Berkeley National Laboratory, Berkeley, CA, 2001.

- Zaytsev, I.D. and Aseyev, G.G. Properties of Aqueous Solutions of Electrolytes, 
  Boca Raton, Florida, USA CRC Press, 1993.

- Engineering ToolBox, `Water - Density, Specific Weight and Thermal Expansion 
  Coefficients <https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html>`__,
  2003


- Engineering ToolBox, `Water - Dynamic (Absolute) and Kinematic Viscosity 
  <https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html>`__,
  2004

- E. E. Michaelides, `Thermodynamic properties of geothermal fluids.
<https://www.osti.gov/biblio/6760030>`__,
Transactions - Geothermal Resources Council, vol. 5, pp. 361-364, 1981.
