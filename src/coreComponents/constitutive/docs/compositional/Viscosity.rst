Viscosity
=========

Viscosity dictates the fluid's resistance to shear flow and fundamentally controls phase mobilities in porous media. The framework provides multiple approaches to evaluate viscosity.

Constant viscosity model
------------------------

The simplest approach applies a static, user-defined viscosity. While the viscosity does not change with state variables such as pressure, temperature, or composition, it is possible to assign different constant viscosities to the different fluid phases. This is recommended solely for idealised simulations, synthetic testing, or when fluid flow is dominated entirely by pressure gradients rather than mobility contrasts.

Parameters
~~~~~~~~~~~~

* ``constantPhaseViscosity``: The constant viscosity value applied to each phase. Unit: [Pa.s].

Lohrenz-Bray-Clark (LBC) viscosity model
----------------------------------------

The Lohrenz, Bray, and Clark correlation (Lohrenz et al., 1964) is the industry-standard methodology for calculating the viscosity of compositional hydrocarbon mixtures. 

The model evaluates viscosity by combining a low-pressure dilute gas contribution with a dense fluid residual term.

Dilute gas viscosity
~~~~~~~~~~~~~~~~~~~~~

The low-pressure, dilute gas viscosity :math:`\mu_{i}^*` of each pure component :math:`i` is estimated as a function of temperature using the correlations of Stiel and Thodos (1961). This relies on the viscosity-reducing parameter :math:`\xi_i`, defined as:

.. math::

    \xi_i = \frac{T_{c,i}^{1/6}}{M_i^{1/2} P_{c,i}^{2/3}}

where :math:`T_{c,i}` is the critical temperature [K], :math:`P_{c,i}` is the critical pressure [atm], and :math:`M_i` is the molar weight [g/mol].

For standard non-polar components, the dilute viscosity (in centipoise) depends on the reduced temperature :math:`T_{r,i} = T / T_{c,i}`:

.. math::

    \mu_{i}^* = \frac{34.0 \times 10^{-5} T_{r,i}^{0.94}}{\xi_i} \quad \text{for } T_{r,i} \le 1.5

.. math::

    \mu_{i}^* = \frac{17.78 \times 10^{-5} (4.58 T_{r,i} - 1.67)^{0.625}}{\xi_i} \quad \text{for } T_{r,i} > 1.5

For hydrogen gas (:math:`M_i < 2.1 \times 10^{-3}` kg/mol), the specific correlation is:

.. math::

    \mu_{H_2}^* = 90.71 \times 10^{-5} (0.1375 T - 1.67)^{0.625}

These pure viscosities are then blended to find the dilute mixture viscosity, :math:`\mu^*`. The user may select between the Herning-Zipperer, Wilke, or Brokaw mixing rules to aggregate the individual components with mole fractions :math:`y_i`.

The Herning-Zipperer mixing rule (Herning and Zipperer, 1936) evaluates the mixture viscosity as:

.. math::

    \mu^* = \frac{\sum_i y_i \mu_i^* \sqrt{M_i}}{\sum_i y_i \sqrt{M_i}}

The Wilke mixing rule (Wilke, 1950) evaluates the mixture viscosity as:

.. math::

    \mu^* = \sum_i \frac{y_i \mu_i^*}{\sum_j y_j \phi_{ij}}

where the interaction parameter :math:`\phi_{ij}` is defined as:

.. math::

    \phi_{ij} = \frac{\left[1 + (\mu_i^*/\mu_j^*)^{1/2} (M_i/M_j)^{-1/4}\right]^2}{\sqrt{8(1 + M_i/M_j)}}

The Brokaw mixing rule (Brokaw, R. S., 1968, Viscosity of Gas Mixtures, National Aeronautics and Space Administration) evaluates the mixture viscosity as:

.. math::

    \mu^* = \sum_i \frac{y_i \mu_i^*}{\sum_j y_j \phi_{ij} \sqrt{\mu_i^*/\mu_j^*}}

where the interaction parameter :math:`\phi_{ij}` depends on the molar weight ratios :math:`A_{ij} = M_i / M_j` and :math:`B_{ij} = [4 M_i M_j / (M_i + M_j)^2]^{0.25}`:

.. math::

    \phi_{ij} = \frac{B_{ij}}{\sqrt{A_{ij}}} \left[ 1 + \frac{A_{ij} - A_{ij}^{0.45}}{2 + 2A_{ij} + B_{ij}(1+A_{ij}^{0.45})/(1+B_{ij})} \right]

Dense fluid contribution
~~~~~~~~~~~~~~~~~~~~~~~~~

A residual viscosity, representing the effects of elevated pressure and density, is calculated using the Lohrenz et al. (1964) correlation. This evaluates a fourth-order polynomial dependent on the reduced density of the phase.

First, the mixture critical parameters are calculated using Kay's mixing rule:

.. math::

    T_c = \sum_i y_i T_{c,i}, \quad P_c = \sum_i y_i P_{c,i}, \quad V_c = \sum_i y_i V_{c,i}, \quad M = \sum_i y_i M_i

The mixture viscosity-reducing parameter :math:`\xi` and reduced density :math:`\rho_r` are defined using the phase molar density :math:`\rho_m`:

.. math::

    \xi = \frac{T_c^{1/6}}{M^{1/2} P_c^{2/3}}, \quad \rho_r = \rho_m V_c

The dense fluid polynomial :math:`f(\rho_r)` is then evaluated:

.. math::

    f(\rho_r) = 0.1023 + 0.023364 \rho_r + 0.058533 \rho_r^2 - 0.040758 \rho_r^3 + 0.0093324 \rho_r^4

This residual is added to the dilute gas viscosity to obtain the final phase viscosity :math:`\mu`:

.. math::

    \mu = \mu^* + \frac{f(\rho_r)^4 - 10^{-4}}{\xi}

This model is strongly recommended for all miscible gas, volatile oil, and general compositional hydrocarbon simulations.

Parameters
~~~~~~~~~~~~

* ``viscosityMixingRule``: Defines the dilute gas mixing rule. Valid options are ``HerningZipperer``, ``Wilke``, or ``Brokaw``.
* ``componentCriticalVolume``: The critical volumes of each component, required for the reduced density calculation. Unit: [m^3/mol].

Phillips brine viscosity model
------------------------------

For aqueous phases, the Phillips brine viscosity model uses empirical correlations from Phillips et al. (1981) to calculate viscosity as a function of temperature and salinity. In this specific implementation, the brine viscosity is evaluated independent of pressure. This model is recommended for saline aquifers and CO2 sequestration.

The model determines the final brine viscosity :math:`\mu` by applying a salinity- and temperature-dependent multiplier to the pure water viscosity :math:`\mu_w`:

.. math::

    \mu = \mu_w(T) \left[ c_0(S) + c_1(S) T \right]

where :math:`T` is the temperature in degrees Celsius and :math:`S` is the salinity in mol/kg. The pure water viscosity :math:`\mu_w(T)` is not evaluated from a closed-form analytical equation; rather, it is interpolated from a pre-computed tabulated function of pure water saturation viscosities corresponding to the given temperature.

The coefficients :math:`c_0(S)` and :math:`c_1(S)` dictate the scaling effect of the dissolved salts and are evaluated as follows:

.. math::

    c_0(S) = 1 + a S + b S^2 + c S^3

.. math::

    c_1(S) = d (1 - \exp(k S))

The correlation uses the specific empirical constants :math:`a = 0.0816`, :math:`b = 0.0122`, :math:`c = 0.000128`, :math:`d = 0.000629`, and :math:`k = -0.7`.

Parameters
~~~~~~~~~~~~

* ``salinity``: The salinity of the brine. Unit: [mol/kg].
* ``temperatureCoordinates``: List of temperature values for the interpolation of tabulated pure water properties. Unit: [K].
