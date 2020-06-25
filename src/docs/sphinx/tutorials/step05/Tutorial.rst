.. _TutorialDeadOilEgg:

########################################################
Tutorial 5: Multiphase flow in the Egg model  
########################################################

**Context**

In this tutorial, we illustrate the concepts presented in :ref:`TutorialDeadOilBottomLayersSPE10`
using the `Egg model <https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/gdj3.21>`_.
We show how to set up a water injection problem relying on a Dead-Oil thermodynamic model.
The wells are placed following the description of the original test case. 

**Objectives**

In this tutorial, we re-use many GEOSX features already presented in
:ref:`TutorialDeadOilBottomLayersSPE10`, with a focus on

- how to import an external mesh with embedded geological properties (permeability) in the GMSH format (``.msh``),
- how to tune the parameters of the linear and nonlinear solvers.

**Input file**

This tutorial is based on the XML file located at

.. code-block:: console

  src/coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml

The mesh file corresponding to the Egg model is stored in the GEOSXDATA repository.
Therefore, you must first download the GEOSXDATA repository in the same folder
as the GEOSX repository to run this test case.
  
------------------------------------
GEOSX input file
------------------------------------

The XML file considered here follows the typical structure of the GEOSX input files:

 #. :ref:`Solver <Solver_tag_dead_oil_egg_model>`
 #. :ref:`Mesh <Mesh_tag_dead_oil_egg_model>`
 #. :ref:`Geometry <Geometry_tag_dead_oil_egg_model>`
 #. :ref:`Events <Events_tag_dead_oil_egg_model>`
 #. :ref:`NumericalMethods <NumericalMethods_tag_dead_oil_egg_model>`
 #. :ref:`ElementRegions <ElementRegions_tag_dead_oil_egg_model>`
 #. :ref:`Constitutive <Constitutive_tag_dead_oil_egg_model>`
 #. :ref:`FieldSpecifications <FieldSpecifications_tag_dead_oil_egg_model>`
 #. :ref:`Outputs <Outputs_tag_dead_oil_egg_model>`

.. _Solver_tag_dead_oil_egg_model:

Solvers: tuning the solution strategy
-------------------------------------

In this tutorial, we use the approach described in :ref:`TutorialDeadOilBottomLayersSPE10`
to couple reservoir flow with wells.
That is, we define a coupling solver of type **CompositionalMultiphaseReservoir**
named ``coupledFlowAndWells``. 
This coupling solver drives the simulation and is in charge of binding the following
single-physics solvers:

 - the single-physics reservoir flow solver, a solver of type **CompositionalMultiphaseFlow** named ``compositionalMultiphaseFlow`` and documented at :ref:`CompositionalMultiphaseFlow`,
 - the single-physics well solver, a solver of type **CompositionalMultiphaseWell** named ``compositionalMultiphaseWell`` and documented at :ref:`CompositionalMultiphaseWell`).

The solver information is specified in the **NonlinearSolverParameters** and
**LinearSolverParameters** XML blocks.
It should be specified within the **CompositionalMultiphaseReservoir** XML block since
the coupling solver drives the solution strategy.
Note that any solver information specified in the single-phase physics XML blocks will
not be taken into account.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_SOLVERS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_SOLVERS_END -->


.. _Mesh_tag_dead_oil_egg_model:

Specifying a reservoir mesh and defining the geometry of the wells
------------------------------------------------------------------

The **Mesh** block consists of two parts:

- The reservoir mesh is described in the **PAMELAMeshGenerator** block,
- The wells are described in the **InternalWell** blocks, with one well per block.

The reservoir mesh is imported from a ``.msh`` file that contains the mesh geometry
and also includes the permeability values in the x, y, and z directions.
These quantities must be specified using the metric unit system, i.e., in meters
for the well geometry and square meters for the permeability field.
We note that the mesh file only contains the active cells, so there is no keyword
needed in the XML file  to define them.

As in :ref:`TutorialDeadOilBottomLayersSPE10`, the geometry of the wells is defined
internally using the description provided in the **InternalWell** XML blocks.
We remind the user that this block must point to the reservoir mesh
(using the attribute ``meshName``), the corresponding well region (using
the attribute ``wellRegionName``), and the corresponding well control
(using the attribute ``wellControlName``).
The well placement implemented here follows the pattern of the original test case.
Additional details about the specification of the well geometry can be
found in :ref:`TutorialDeadOilBottomLayersSPE10`.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_MESH -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_MESH_END -->


.. _Geometry_tag_dead_oil_egg_model:

Geometry tag
-----------------

The **Geometry** XML block was used in the single-phase tutorials to specify boundary conditions.
Since we use wells and assume no-flow boundary conditions in this tutorial, the **Geometry**
block is not needed.

.. _Events_tag_dead_oil_egg_model:

Specifying events
------------------------

In the **Events** XML block, we specify two types of **PeriodicEvents**.

The periodic event named ``solverApplications`` notifies GEOSX that the
coupled solver ``coupledFlowAndWells`` has to be applied to its target
regions (here, reservoir and wells) at every time step.
The time stepping strategy has been fully defined in the **CompositionalMultiphaseReservoir**
coupling block using the ``initialDt`` attribute and the **NonlinearSolverParameters**
nested block.

As in :ref:`TutorialDeadOilBottomLayersSPE10`, we also define an output event
instructing GEOSX to write out ``.vtk`` files at the time frequency specified
by the attribute ``timeFrequency``.
The ``target`` attribute must point to the **VTK** sub-block of the **Outputs**
block (defined at the end of the XML file) by name (here, ``vtkOutput``).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_EVENTS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_EVENTS_END -->


.. _NumericalMethods_tag_dead_oil_egg_model:

Defining Numerical Methods
----------------------------------

In the ``NumericalMethods`` XML block, we instruct GEOSX to use a TPFA finite-volume
numerical scheme.
This part is identical to the corresponding section of :ref:`TutorialDeadOilBottomLayersSPE10`.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_NUMERICAL_METHODS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_NUMERICAL_METHODS_END -->


.. _ElementRegions_tag_dead_oil_egg_model:

Defining reservoir and well regions
-----------------------------------

In this section of the input file, we follow the procedure already described in
:ref:`TutorialDeadOilBottomLayersSPE10` for the definition of the reservoir and well
regions.

We associate a **CellElementRegion** named ``reservoir`` to the reservoir mesh.
Since we have imported a mesh with one region consisting of hexahedral cells, we
must set the attribute ``cellBlocks`` to ``0_HEX``.
If you use a name that is not ``0_HEX`` for this attribute, GEOSX will throw an error
at the beginning of the simulation.

We also associate a **WellElementRegion** to each well. As the **CellElementRegion**,
it contains a ``materialList`` that must point (by name) to the constitutive models
defined in the **Constitutive** XML block.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_ELEMENT_REGIONS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_ELEMENT_REGIONS_END -->


.. _Constitutive_tag_dead_oil_egg_model:

Defining material properties with constitutive laws
---------------------------------------------------------------------

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_CONSTITUTIVE_END -->


.. _FieldSpecifications_tag_dead_oil_egg_model:

Defining properties with the FieldSpecifications
---------------------------------------------------------------------

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_FIELD_SPECS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_FIELD_SPECS_END -->


.. _Outputs_tag_dead_oil_egg_model:

Specifying the output formats
----------------------------------

In this section, we request an output of the results in VTK format and an output of the restart file.
Note that the names defined here must match the names used in the **Events** XML block to define the output frequency.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/Egg/dead_oil_egg.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_OUTPUT -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_OUTPUT_END -->


All elements are now in place to run GEOSX.

------------------------------------
Runnning GEOSX
------------------------------------

In progress

------------------------------------
Visualization of results
------------------------------------

A file compatible with Paraview is produced in this tutorial.
It is found in the output folder, and usually has the extension `.pvd`.
We can load this file into Paraview directly and visualize results:


------------------------------------
To go further
------------------------------------

**Feedback on this tutorial**

This concludes the tutorial on setting up a Dead-Oil simulation in the Egg model.
For any feedback on this tutorial, please submit
a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.

**Next tutorial**

In the next tutorial :ref:`TutorialCO2FieldCaseUnstructuredGrid`, we learn how to run a
more complex test case based on an unstructured mesh.

**For more details**

  - A complete description of the reservoir flow solver is found here: :ref:`CompositionalMultiphaseFlow`.
  - The well solver is description at :ref:`CompositionalMultiphaseWell`. 
  - The available constitutive models are listed at :ref:`Constitutive`.

