.. _TemperatureDependentThermalExpansionCoefficient:


############################################################
Temperature-dependent Thermal Expansion Coefficient Model
############################################################


Overview
======================

This model gives the drained linear thermal expansion coefficient :math:`\alpha` of the solid as a
function of temperature, either as a linear law:

.. math::
   \alpha(T) = \alpha_{0} + \frac{ d\alpha }{ dT } (T - T_{0})

where
      :math:`\alpha_{0}` is the coefficient at the reference temperature :math:`T_{0}`;
      :math:`\frac{ d\alpha }{ dT }` is its derivative with respect to temperature (zero for a constant coefficient);

or as a one-dimensional ``TableFunction`` of temperature, with linear interpolation between the
table points and a constant value outside the table.

:math:`T_{0}` is a property of the material (the temperature at which :math:`\alpha = \alpha_{0}`),
not the initial temperature of the simulation.


Thermal strain
======================

With a temperature-dependent coefficient, the thermal strain is the integral of :math:`\alpha`:

.. math::
   \boldsymbol{\varepsilon}^{th} = \left( \int_{T_{init}}^{T} \alpha(s) \, ds \right) \mathbf{I}

Over a time step from :math:`T_{n}` to :math:`T`, the thermal strain increment
:math:`\left( \int_{T_{n}}^{T} \alpha(s) \, ds \right) \mathbf{I}` is subtracted from the strain increment
before the solid model is updated, so that the solid model (elastic, plastic or damage) only sees the
mechanical strain. The integral is computed exactly for both the linear and the tabulated laws, so the
result does not depend on the time step. Its derivative with respect to :math:`T` is :math:`\alpha(T)`,
which gives the temperature derivative of the stress used in the Jacobian.

When the coefficient does not depend on temperature, these expressions reduce to the usual
:math:`\alpha \, (T - T_{n})`.


Porosity
======================

In Biot-Coussy theory, the thermal term of the porosity change uses the same solid coefficient as the
stress. When the model is attached, the ``<BiotPorosity>`` update therefore uses it as well: in the
fully implicit and in the sequential (fixed-stress) schemes, each term :math:`\alpha_{\phi} \, \Delta T`
of the porosity update is replaced by the exact integral of :math:`\alpha` over the same temperature
interval, and :math:`\alpha(T)` is used in the temperature derivatives. The porosity coefficient
``defaultPorosityTEC`` must then not be given. Without the model, the porosity keeps using
``defaultPorosityTEC``, as before.


Parameters
======================

The model is defined in a ``<SolidThermalExpansion>`` block of the input XML file and must be
assigned a unique name via the ``name`` attribute. This name is used to attach the model to a porous
solid (for example ``<PorousElasticIsotropic>``) via its ``solidThermalExpansionModelName`` attribute.
When this attribute is omitted, the porous solid uses the constant coefficient
``defaultDrainedLinearTEC`` of its solid model, as before. The coefficient cannot be given in both places,
and ``defaultPorosityTEC`` of the porosity model cannot be given with the model either (see Porosity).

Either ``referenceDrainedLinearTEC`` (linear law) or ``drainedLinearTECTableName`` (tabulated law) must be given, but not both.
``referenceTemperature`` is required when ``dDrainedLinearTEC_dTemperature`` is not zero.
The table must have a single input (temperature) and use linear interpolation.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/SolidThermalExpansion.rst


Limitations
======================

- The model is supported by the porous solids (``Porous*`` models, including the damage ones).
  Other coupled solids reject it.
- The thermal poromechanics solver with embedded fractures only accepts a coefficient that does not
  depend on temperature.


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <PorousElasticIsotropic
        name="rock"
        solidModelName="rockSolid"
        porosityModelName="rockPorosity"
        permeabilityModelName="rockPerm"
        solidThermalExpansionModelName="rockTEC"/>

      <!-- linear law -->
      <SolidThermalExpansion
        name="rockTEC"
        referenceDrainedLinearTEC="2.0e-6"
        dDrainedLinearTEC_dTemperature="1.64e-9"
        referenceTemperature="0.0"/>
      ...
   </Constitutive>

or, with a table:

.. code-block:: xml

   <Constitutive>
      ...
      <SolidThermalExpansion
        name="rockTEC"
        drainedLinearTECTableName="tecTable"/>
      ...
   </Constitutive>

   <Functions>
      <TableFunction
        name="tecTable"
        inputVarNames="{ temperature }"
        coordinates="{ 0.0, 20.0, 90.0 }"
        values="{ 2.0e-6, 5.28e-6, 1.676e-5 }"/>
   </Functions>

A validation example is given in :ref:`AdvancedExampleNonLinearThermalOedometric`.
