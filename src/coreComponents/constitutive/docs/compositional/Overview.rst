.. _Overview:

Overview
=========================

This model represents a full composition description of a multiphase multicomponent fluid.
Phase behavior is modeled by an equation of state (EOS) and partitioning of components into
phases is computed based on instantaneous chemical equilibrium via a two-phase flash.
Each component (species) is characterized by molar weight and critical properties that
serve as input parameters for the EOS.
See `Petrowiki`_ for more information.

In this model the fluid is described by :math:`N_c` components with :math:`z_c` being the total
mole fraction of component :math:`c`. The fluid can partition into a liquid phase, denoted :math:`\ell`,
and a vapor phase denoted by :math:`v`. Therefore, by taking into account the molar phase component
fractions, (which is the fraction of the molar mass of phase :math:`p` represented by component 
:math:`c`), the following partition matrix establishes the component distribution within the two
phases:

.. math::
    \begin{bmatrix}
    x_{1} & x_{2} & x_{3} & \cdots & x_{N_c} \\
    y_{1} & y_{2} & y_{3} & \cdots & y_{N_c} \\
    \end{bmatrix}

where :math:`x_c` is the mole fraction of component :math:`c` in the liquid phase and :math:`y_c`
is the mole fraction of component :math:`c` in the vapor phase.

The fluid properties are updated through the following steps:

1) The phase fractions (:math:`\nu_p`) and phase component fractions (:math:`x_c` and :math:`y_c`) are
computed as a function of pressure (:math:`p`), temperature (:math:`T`) and total component fractions
(:math:`z_c`).

2) The phase densities (:math:`\rho_p`) and phase viscosities (:math:`\mu_p`) are computed as a function
of pressure, temperature and the updated phase component fractions.

After calculating the phase fractions, phase component fractions, phase densities, phase viscosities,
and their derivatives with respect to pressure, temperature, and component fractions, the
:ref:`CompositionalMultiphaseFlow` then moves on to assembling the accumulation and flux terms.
