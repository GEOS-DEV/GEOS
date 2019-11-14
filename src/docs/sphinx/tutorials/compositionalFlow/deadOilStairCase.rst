Dead-oil staircase model
=================================
We consider a homogeneous reservoir with the staircase geometry shown in the
following figure.

Three fluid phases are present, water (w), oil (o) and gas (g) and a dead-oil
fluid model is considered  as described in :ref:`Black oil fluid model <BlackOilFluid>`.

Input file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The input file for this tutorial can be found at
```
GEOSX/examples/compositionalFlow/deadoilStaircase.xml
```
Let us now have a look at each section of the xml input file.

Solvers tag
---------------------------------------
This 3D test case involves gravity so we specify the components of the gravity vector.
Additionally, since we want to run a compositional flow simualtion we add a `CompositionalMultiphaseFlow` solver
to the `Solvers` section. A detail descriptions of all fields specific to a
compositional solver can be found ... Here, we need to make sure to specify the
region to which the solver has to be applied and the specific fluid and solid names.
Solver parameters are defined in a separate subsection.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_SOLVER_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_SOLVER_BLOCK -->

Geometry and Mesh tag
---------------------------------------
The geometry tag can be used to define the location of source as sink terms. They
are identified by defining two geometrical objects, boxes, with specific extensions.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_GEOM_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_GEOM_BLOCK -->

In this test case we construct an internal mesh formed by 16 cell blocks each one containing
5 x 5 x 3 cells.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
    :language: xml
    :start-after: <!-- START_SPHINX_INCLUDE_MESH_BLOCK -->
    :end-before: <!-- END_SPHINX_INCLUDE_MESH_BLOCK -->

Element regions tag
----------------------------------------
The cell blocks defined in the previous tag are here assigned to a specific region.
Here, we consider two regions, one identifying the high permeability channel in
which the flow occurs and one identifying the flow barriers.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_ELEMREG_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_ELEMREG_BLOCK -->

Events tag
---------------------------------------

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_EVENTS_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_EVENTS_BLOCK -->

Numerical methods tag
----------------------------------------

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
   :language: xml
   :start-after: <!-- START_SPHINX_INCLUDE_NUMMET_BLOCK -->
   :end-before: <!-- END_SPHINX_INCLUDE_NUMMET_BLOCK -->

Constitutive tag
----------------------------------------

.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
    :language: xml
    :start-after: <!-- START_SPHINX_INCLUDE_CONST_BLOCK -->
    :end-before: <!-- END_SPHINX_INCLUDE_CONST_BLOCK -->

Field Specification tag
----------------------------------------

Permeability and porosity
``````````````````````````````````````````
.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
    :language: xml
    :start-after: <!-- START_SPHINX_INCLUDE_PERM_BLOCK -->
    :end-before: <!-- END_SPHINX_INCLUDE_PERM_BLOCK -->

Initial conditions
``````````````````````````````````````````
.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
    :language: xml
    :start-after: <!-- START_SPHINX_INCLUDE_INIT_BLOCK -->
    :end-before: <!-- END_SPHINX_INCLUDE_INIT_BLOCK -->

Boundary conditions
``````````````````````````````````````````
.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
      :language: xml
      :start-after: <!-- START_SPHINX_INCLUDE_BC_BLOCK -->
      :end-before: <!-- END_SPHINX_INCLUDE_BC_BLOCK -->

Output tags
----------------------------------------
.. literalinclude:: ../../../../coreComponents/physicsSolvers/FiniteVolume/integratedTests/compositionalMultiphaseFlow/deadoil_3ph_staircase_3d.xml
    :language: xml
    :start-after: <!-- START_SPHINX_INCLUDE_OUTPUT_BLOCK -->
    :end-before: <!-- END_SPHINX_INCLUDE_OUTPUT_BLOCK -->
Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
