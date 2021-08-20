.. _ParallelPlatesPermeability:


####################################################
Parallel Plates Permeability Model
####################################################


Overview
======================

The parallel plates permeability model defines the relationship between the effective fracture aperture and its corresponding permeability. In the model, the two fracture walls are assumed to
be smooth and parallel to each other and separated by a uniform aperture.

.. math::
   k =  \frac{a^3}{12} 
  
where :math:`a` denotes the effective fracture aperture.

Please note that :math:`k` dimensionally is not a permeability.



Parameters
======================

The Parallel Plates Permeability Model can be called in the
``<Constitutive>`` block of the input XML file.
This permeability model must be assigned a unique name via
``name`` attribute.
This name is used to assign the model to regions of the physical
domain via a ``materialList`` attribute in the ``<ElementRegions>``
block.

The following attributes are supported:

.. include:: /coreComponents/schema/docs/ParallelPlatesPermeability.rst



Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <ParallelPlatesPermeability name="fracPerm"/>
      ...
   </Constitutive>


