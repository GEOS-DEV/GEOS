.. _PowerLawPermeability:


####################################################
Power Law Permeability Model
####################################################


Overview
======================

In the power law model (see `ref`_), the permeability is related to the porosity through a
phenomenological relation

.. math::
   k = k_{min} + (k_0 - k_{min}) f(\phi, \phi_0, \phi_c, n),

with

.. math::
   f = \left( \frac{\phi - \phi_c}{\phi_0 - \phi_c} \right)^n,

where :math:`\phi` is the porosity of the porous medium, :math:`\phi_0` the reference porosity
at which the permeability is the reference permeability :math:`k_0`, :math:`\phi_c` the critical
porosity below which the pore space no longer percolates, and :math:`n` the exponent of the power
law. The residual permeability :math:`k_{min}`, zero by default, is the floor the permeability
decays to once the pore space stops percolating.

This model is primarily intended for reactive transport simulations, in which the porosity evolves
with the dissolution and the precipitation of the minerals.


Parameters
======================

The Power Law Permeability Model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via the
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /docs/sphinx/datastructure/PowerLawPermeability.rst



Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <PowerLawPermeability name="matrixPerm"
                            referencePermeabilityComponents="{1.0e-15, 1.0e-15, 1.0e-15}"
                            referencePorosity="0.2"
                            criticalPorosity="0.0"
                            exponent="3.0"
                            minPermeability="1.0e-20"/>
      ...
   </Constitutive>

.. _ref: https://documentation.pflotran.org/theory_guide/mode_reactive_transport.html#changes-in-material-properties
