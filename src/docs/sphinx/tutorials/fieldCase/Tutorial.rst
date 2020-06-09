.. _TutorialFieldCase:

#########################################
Simple field cases
#########################################

**Context**

In this tutorial, we set up a simple field case for single phase flow simulation (see :ref:`SinglePhaseFlow`).


**Objectives**

At the end of this tutorial you will know:

  - how to import external mesh information and properties,
  - how to run a specific solver (here, flow) in a specific region only,
  - the basic syntax of a solver block for wells,
  - how to control output and visualize results.


**Input file**

The xml input file for this test case is located at:

.. code-block:: console

  src/coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml


*****************
Domain definition
*****************

We consider the following mesh as a numerical support to the simulations in this tutorial:

.. image:: mesh.png
   :width: 400px

This mesh contains three continuous regions:

  - a Top region (overburden, elementary tag = 1)
  - a Middle region (reservoir layer, elementary tag = 2)
  - a Bottom region (underburden, elementary tag = 3)


The mesh is defined using the GMSH file format (see :ref:`Meshes` for more information on
the supported mesh file format). Each tetrahedron is associated to a unique tag.

.. note::

  The GMSH file format starts numbering tags from 1. In GEOSX, the numbering
  of ``CellElementRegion`` starts from 0. As a consequence, when regions
  are defined in GEOSX, we substract 1 from the GMSH tag to refer
  to the same region. The next sections of this tutorial illustrates this.



***************************
Importing the mesh in GEOSX
***************************

Here, we use the ``PAMELAMeshGenerator`` to load the mesh (see :ref:`ImportingExternalMesh`).
The syntax to import external meshes is simple : in the XML file,
the mesh ``file`` is included with its path relative to the location of the input XML file
(or absolute)
and a user-specified ``name`` for the mesh object.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_MESH -->
  :end-before: <!-- SPHINX_FIELD_CASE_MESH_END -->


************************
Running flow simulations
************************



We now run a single-phase flow simulations on this field.
We assume that the overburden and the underburden are impermeable,
and flow only happen in the reservoir.

There are two methods to achieve this regional solve.

The first solution is to define a unique ``CellElementRegion`` corresponding to the reservoir.

.. code-block:: xml

  <ElementRegion>
    <CellElementRegion name="ReservoirLayer"
                       cellBlocks="{1_TETRA}"
                       materialList="{water, rock}">
  </ElementRegion>

The second solution is to define all the ``CellElementRegions`` as they are in the GMSH file,
but defining the solvers only on the reservoir layer. In this case, the ``ElementRegion`` tag
is :

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_REGION -->
  :end-before: <!-- SPHINX_FIELD_CASE_REGION_END -->

We choose the latest in order to keep the over- and underburden for visualization.

.. note::
  The material list here was set for a single phase flow problem. This list is subject
  to change if the problem is not a single phase flow problem.


Single Phase Flow, no well
==========================

Now, we demonstrate how to run a basic flow simulation in the reservoir layer.
We do not consider any coupling with wells. Injection and production will be specified by
imposing a high pressure in the cells close to the injection area and a low pressure
in the cells close to the production area.

Material definition
-------------------

We simulate a single phase flow in the reservoir layer with two types of materials
: fluid (water) and solid (rock).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_FIELD_CASE_CONSTITUTIVE_END -->

The constitutive parameters such as the density, the viscosity, the compressibility etc. can
be modified here using the International System of Units.

.. note::
  To consider an incompressible fluid, the user has to set the compressibility to 0.

Numerical methods settings
--------------------------

Once the materials are defined, we define the numerical method that will be used in the solver.
We propose to use the classical two-point flux approximation scheme.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_NUMERICAL -->
  :end-before: <!-- SPHINX_FIELD_CASE_NUMERICAL_END -->

The ``TwoPointFluxApproximation`` node has to contain the primary field to be solved in
``fieldName``. For a flow problem it is the pressure.

Solver settings
---------------

The ``Solver`` XML tag is then set.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_SOLVER -->
  :end-before: <!-- SPHINX_FIELD_CASE_SOLVER_END -->

This node is crucial as it gather all the information previously defined. We use the classical
``SinglePhaseFVM`` (FVM for Finite Volume Method), with the two-points flux approximation
previously defined in the ``NumericalMethod`` tag. The ``targetRegions`` refers only
to the Reservoir, because we just want to solve the flow in this region. The ``fluidNames``
and ``solidNames`` refers to the previously defined materials in the ``Constitutive`` tag.

The ``NonlinearSolverParameters`` and ``SystemSolverParameters`` are then used to set the
numerical solver parameters such as the Newton tolerance and the maximum number of
iterations.

Defining geometry boxes
-----------------------

The geometry boxes are useful in order to flag some cells or nodes within it. We will use it
this GEOSX feature to mimic the injection and production.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_GEOMETRY -->
  :end-before: <!-- SPHINX_FIELD_CASE_GEOMETRY_END -->

In order to define a box, the user had to define ``xMax`` and ``xMin`` that are two opposite nodes of the box.

.. note::
  In the previous XML snippet, the box ``all`` includes the whole model. This box has to be defined for every problem
  involving the import of an external mesh

Field specification
-------------------

The next step is to specify fields, including:

  - The initial value (here, the pressure has to be initialized)
  - The static properties (here, we have to define the permeability tensor and the porosity)
  - The boundary conditions (here, the injection and production pressure have to be set)

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_GEOMETRY -->
  :end-before: <!-- SPHINX_FIELD_CASE_GEOMETRY_END -->

You may note :

 - All static parameters and initial value fields have to present the ``initialCondition`` field set to ``1``.
 - The ``objectPath`` refers to the ``ÈlementRegion`` in which the field value has
 - The ``setName`` field points to the box previously defined to apply the fields. 
 - Tensors values has to be set component by component. GEOSX deals with a diagonal matrix for
   the permeability, so the three values of the permeability tensor are set individually
   using ``component`` field.
 - ``name`` and ``fieldName`` have a different meaning. ``name`` is used to give a name to the XML block. It has
   to be unique. ``fieldName`` is the name of the field register in GEOSX. This value has to be set carefully
   in accordance with the documentation of each solver.

Defining output
---------------

The ``Output`` XML tag is used to trigger the writing of visualization files. Here we propose
to write files natively readable by Paraview.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_OUTPUT -->
  :end-before: <!-- SPHINX_FIELD_CASE_OUTPUT_END -->

.. note::
  The ``name`` keyword defines the name of the output file.

Triggering events
-----------------

The final steps is to organize events so the simulation and the outputs can be triggered.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_EVENTS -->
  :end-before: <!-- SPHINX_FIELD_CASE_EVENTS_END -->

The ``Events`` tag is associated with the ``maxTime`` keyword defining the maximum time. If this time is
ever reached or exceeded, the simulation will end.

Two ``PeriodicEvent`` are then defined. The first one is associated with the solver. The  ``forceDt`` keyword
means that there will always be time-steps of 20 seconds. The second is associated with the output. The ``timeFrequency``
keyword means that it will be executed every 100 seconds. The ``targetExactTimestep`` is set to 1, meaning that
the Event Manager will impose this event will be triggered exactly every 100 seconds.

Lauching the simulation
-----------------------

The simulation can be launched with:

.. code-block:: console

  geosx -i FieldCaseTutorial1.xml

