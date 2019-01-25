############################################
Black-oil fluid model
############################################

Black-oil model represents the fluid with 2 hydrocarbon components (oil and gas)
and water component, partitioned across 3 phases (also named oil, gas and water).
Phase behavior is characterized by black-oil tables, that tabulate saturated and
undersaturated oil and gas properties as functions of pressure and solution ratios.
In **dead-oil** each component occupies only one phase, whereas in **live-oil**
model gas component can dissolve in oil phase and oil component can evaporate
into gas phase.
See `Petrowiki`_ for more information.

Both types are represented by ``<BlackOilFluid>`` node in the input. Under the hood
this is a wrapper around ``PVTPackage`` library, which is included as a submodule.
In order to use the model, GEOSX must be build with ``ENABLE_PVTPACKAGE=ON``.

The following attributes are supported:

.. include:: /../coreComponents/fileIO/schema/docs/BlackOilFluid.rst

References
=========================

.. _Petrowiki: https://petrowiki.org/Phase_behavior_in_reservoir_simulation#Black-oil_PVT_models