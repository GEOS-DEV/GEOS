.. _BlackOilFluid:

############################################
Black-oil fluid model
############################################

Overview
=========================

Black-oil model represents the fluid with 2 hydrocarbon components (oil and gas)
and water component, partitioned across 3 phases (also named oil, gas and water).
Phase behavior is characterized by black-oil tables, that tabulate saturated and
undersaturated oil and gas properties as functions of pressure and solution ratios.
In **dead-oil** each component occupies only one phase, whereas in **live-oil**
model gas component can dissolve in oil phase and oil component can evaporate
into gas phase.
See `Petrowiki`_ for more information.

Usage
=========================

Both types are represented by ``<BlackOilFluid>`` node in the input.
Under the hood this is a wrapper around ``PVTPackage`` library, which is included as a submodule.
In order to use the model, GEOSX must be built with ``-DENABLE_PVTPACKAGE=ON`` (default).

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/BlackOilFluid.rst

Supported phase names are:

* oil
* gas
* water

Input example
=========================

.. code-block:: xml

  <Constitutive>
    <BlackOilFluid name="fluid1"
                   fluidType="LiveOil"
                   phaseNames="oil gas water"
                   surfaceDensities="800.0 0.9907 1022.0"
                   componentMolarWeight="114e-3 16e-3 18e-3"
                   tableFiles="pvto.txt pvtg.txt pvtw.txt"/>
  </Constitutive>


.. _Petrowiki: https://petrowiki.org/Phase_behavior_in_reservoir_simulation#Black-oil_PVT_models