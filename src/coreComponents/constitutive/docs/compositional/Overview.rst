Overview
========

The compositional fluid models implemented in this framework provide a robust and flexible approach for modelling the thermodynamic behaviour of multicomponent, multiphase fluid mixtures. These models are intended to support complex subsurface flow simulations, including primary hydrocarbon recovery, miscible gas injection, carbon sequestration, and geothermal energy extraction.

Fluid composition and phase partitioning
----------------------------------------

In this model, the fluid is described by :math:`N_c` components, with :math:`z_c` being the total mole fraction of component :math:`c`. The fluid can partition into up to three phases: a liquid phase (denoted :math:`\ell`), a vapour phase (denoted :math:`v`), and potentially an aqueous phase (denoted :math:`a`).

Therefore, by taking into account the molar phase component fractions (which is the fraction of the molar mass of phase :math:`p` represented by component :math:`c`), the following partition matrix establishes the component distribution within a two-phase (liquid-vapour) system:

.. math::

    \begin{bmatrix}
    x_1 & x_2 & x_3 & \cdots & x_{N_c} \\
    y_1 & y_2 & y_3 & \cdots & y_{N_c}
    \end{bmatrix}

where :math:`x_c` is the mole fraction of component :math:`c` in the liquid phase and :math:`y_c` is the mole fraction of component :math:`c` in the vapour phase. 

If a third aqueous phase is also present, with component mole fractions denoted by :math:`w_c`, this partition matrix naturally extends to three rows to fully describe the component distribution across all three phases:

.. math::

    \begin{bmatrix}
    x_1 & x_2 & x_3 & \cdots & x_{N_c} \\
    y_1 & y_2 & y_3 & \cdots & y_{N_c} \\
    w_1 & w_2 & w_3 & \cdots & w_{N_c}
    \end{bmatrix}

Structure of the fluid property calculations
--------------------------------------------

The evaluation of fluid properties follows a sequential thermodynamic structure designed to ensure consistency across all phase properties. At a given pressure, temperature, and overall composition, the framework executes the following overarching evaluations:

1. State Conversion: If the system is formulated using mass variables, the overall mass fractions are converted into overall mole fractions.
2. Phase Split (Flash Calculation): The thermodynamic stability of the mixture is evaluated. If unstable, a flash calculation splits the mixture into separate phases (e.g., liquid, vapour, and potentially an aqueous phase), determining the phase fractions and the composition of each individual phase.
3. Phase Properties: Using the computed phase compositions, the model evaluates properties for each specific phase:
   
   * Density: Evaluated typically via a cubic Equation of State (EoS) combined with volume shift corrections.
   * Viscosity: Evaluated using established compositional correlations or empirical models.
   * Thermal Properties: If the simulation is non-isothermal, the phase enthalpies and internal energies are computed using ideal gas heat capacities and EoS departure functions.

4. Total Fluid Properties: The individual phase properties are homogenised to yield total fluid properties (e.g., total density). Derivatives with respect to pressure, temperature, and composition are rigorously propagated using the chain rule.
