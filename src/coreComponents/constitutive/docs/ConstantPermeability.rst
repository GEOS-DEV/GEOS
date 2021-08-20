.. _KozenyCarmanPermeability:


####################################################
Constant Permeability Model
####################################################


Overview
======================

This model is used to define a diagonal permeability tensor that does not depend on any primary variable.

Parameters
======================


The following attributes are supported:

.. include:: /coreComponents/schema/docs/CarmanKozenyPermeability.rst


Example
=======================

.. code-block:: xml

   <Constitutive>
      ...
      <ConstantPermeability name="matrixPerm"
                            permeabilityComponents="{}""  />
      ...
   </Constitutive>
