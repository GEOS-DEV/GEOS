.. _CompositionalMultiphaseFluid:

############################################
Compositional multiphase fluid model
############################################

.. include:: compositional/Overview.rst

Structure of this Documentation
-------------------------------

This documentation is divided into the following sections:

* `Phase Split`_: Details the thermodynamic stability and flash calculations (Negative Flash and K-Value Flash), including the underlying Equations of State.
* `Density`_: Describes the models used to compute molar and mass phase densities.
* `Viscosity`_: Details the models used to evaluate the resistance to flow for each phase.
* `Enthalpy`_: Explains the calculation of thermal properties.
* `Miscellaneous`_: Outlines internal unit conversions and mass-to-molar variable handling.
* `References`_: Provides the bibliography for the scientific correlations implemented in the code.

.. include:: compositional/PhaseSplit.rst

.. include:: compositional/EquationsOfState.rst

.. include:: compositional/NegativeFlash.rst

.. include:: compositional/KValueFlash.rst

.. include:: compositional/ImmiscibleWaterFlash.rst

.. include:: compositional/Density.rst

.. include:: compositional/Viscosity.rst

.. include:: compositional/Enthalpy.rst

.. include:: compositional/Miscellaneous.rst

.. include:: compositional/References.rst
