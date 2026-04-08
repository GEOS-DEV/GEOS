.. _PermeabilitySpecification:

####################################################
Permeability specification
####################################################

Overview
======================

**PermeabilitySpecification** is an optional, higher-level XML tag that you can place in the **FieldSpecifications** block. 
It describes a 3-axis permeability on one or more element regions.
After the input deck is read, GEOS expands each **PermeabilitySpecification** into several **FieldSpecification** objects (one per region and per axis) that the rest of the code already understands.

This is convenient when you want the same **setNames**, **fieldName**, and other attributes for all three permeability components, instead of repeating three nearly identical **FieldSpecification**.

Examples
===============

The following illustrates a **PermeabilitySpecification** for a single region.

.. code-block:: xml

   <FieldSpecifications>
     ...
     <PermeabilitySpecification
       name="perm"
       setNames="{ all }"
       regionNames="{ reservoir/block1 }"
       fieldName="rockPerm_permeability"
       functionName="permFunc"
       scales="{ 9.869233e-16, 9.869233e-16, 9.869233e-16 }"/>
     ...
   </FieldSpecifications>     

The following illustrates a **PermeabilitySpecification** for multiple regions.

.. code-block:: xml

   <FieldSpecifications>
     ...
     <PermeabilitySpecification
       name="perm"
       setNames="{ all }"
       regionNames="{ reservoir/block1, reservoir/block2 }"
       fieldName="rockPerm_permeability"
       functionName="permFunc"
       scales="{ 9.869233e-16, 9.869233e-16, 9.869233e-16 }"/>
     ...
   </FieldSpecifications>     
