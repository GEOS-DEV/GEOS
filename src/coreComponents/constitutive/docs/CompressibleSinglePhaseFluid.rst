############################################
Compressible single phase fluid model
############################################

Overview
=========================

This model represents a compressible single-phase fluid with constant compressibility
and pressure-dependent viscosity.
These assumptions are valid for slightly compressible fluids, such as water, and some
types of oil with negligible amounts of dissolved gas.

Usage
=========================

The model is represented by ``<CompressibleSinglePhaseFluid>`` node in the input.

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/CompressibleSinglePhaseFluid.rst

Example
=========================

.. code-block:: xml

  <Constitutive>
    <CompressibleSinglePhaseFluid name="water"
                                  referencePressure="2.125e6"
                                  referenceDensity="1000"
                                  compressibility="1e-19"
                                  referenceViscosity="0.001"
                                  viscosibility="0.0"/>
  </Constitutive>
