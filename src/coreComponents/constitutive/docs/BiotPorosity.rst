.. _BiotPorosity:


####################################################
Biot Porosity Model
####################################################


Overview
======================

According to poroelasticity theory, the porosity (pore volume) change with respect to pressure is expressed via the following Biot Porosity Model:

.. math::
   \frac{d{\phi}}{dp} = {\alpha}{\frac{d{\epsilon}}{dp} + \frac{\alpha - \phi}{K_s}

where :math:`\alpha` is the Biot coefficient, :math:`\epsilon` is the volumetric strain, :math:`K_s` is the bulk modulus of the grains and :math:`p` is the pressure.


Parameters
======================

The Biot Porosity Model can be called in the
``<Constitutive>`` block of the input XML file.
This porosity model must be assigned a unique name via the
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
