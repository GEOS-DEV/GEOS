.. _CompositionalMultiphaseFluid:

############################################
Compositional multiphase fluid model
############################################

Overview
=========================

This model represents a full composition description of a multiphase multicomponent fluid.
Phase behavior is modeled by an Equation of State (EOS) and partitioning of components into
phases is computed based on instantaneous chemical equilibrium via a two- or three-phase flash.
Each component (species) is characterized by molar weight and critical properties that
serve as input parameters for the EOS.
See `Petrowiki`_ for more information.

Parameters
=========================

The model represented by ``<CompositionalMultiphaseFluid>`` node in the input.
Under the hood this is a wrapper around ``PVTPackage`` library, which is included as a submodule.
In order to use the model, GEOS must be built with ``-DENABLE_PVTPACKAGE=ON`` (default).

The following attributes are supported:

.. include:: ../../../coreComponents/schema/docs/CompositionalMultiphaseFluid.rst

Supported phase names are:

===== ===========
Value Comment
===== ===========
oil   Oil phase
gas   Gas phase
water Water phase
===== ===========

Supported Equation of State types:

===== =======================
Value Comment
===== =======================
PR    Peng-Robinson EOS
SRK   Soave-Redlich-Kwong EOS
===== =======================

Example
=========================

.. code-block:: xml

  <Constitutive>
    <CompositionalMultiphaseFluid name="fluid1"
                                  phaseNames="{ oil, gas }"
                                  equationsOfState="{ PR, PR }"
                                  componentNames="{ N2, C10, C20, H2O }"
                                  componentCriticalPressure="{ 34e5, 25.3e5, 14.6e5, 220.5e5 }"
                                  componentCriticalTemperature="{ 126.2, 622.0, 782.0, 647.0 }"
                                  componentAcentricFactor="{ 0.04, 0.443, 0.816, 0.344 }"
                                  componentMolarWeight="{ 28e-3, 134e-3, 275e-3, 18e-3 }"
                                  componentVolumeShift="{ 0, 0, 0, 0 }"
                                  componentBinaryCoeff="{ { 0, 0, 0, 0 },
                                                        { 0, 0, 0, 0 },
                                                        { 0, 0, 0, 0 },
                                                        { 0, 0, 0, 0 } }"/>
  </Constitutive>


.. _Petrowiki: https://petrowiki.spe.org/Phase_behavior_in_reservoir_simulation#Equation-of-state_models
