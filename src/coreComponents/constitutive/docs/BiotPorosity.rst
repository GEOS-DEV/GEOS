.. _BiotPorosity:


####################################################
Biot Porosity Model
####################################################


Overview
======================

Accodirng to the poroelasticity theory, the porostiy (pore volume) change with respect to pressure is expressed as the following Biot Porosity Model:

.. math::
   \phi = phi_{ref} + \alpha \nabla ( \epsilon_v - \epsilon_{ref} ) + (p - p_{ref}) / N;

where :math:`\alpha` is the Biot coefficient, :math:`\epsilon` is the strain, :math:`p` is the fluid pressure
and :math:`N = \frac{K_s}{b - \phi_{ref}}`


Parameters
======================

The Biot Porosity Model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /coreComponents/schema/docs/BiotPorosity.rst

Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <BiotPorosity name="rockPorosity"
                    grainBulkModulus="1.0e27"
                    defaultReferencePorosity="0.3"/>
      ...
   </Constitutive>
