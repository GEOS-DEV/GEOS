.. _TutorialFieldCase:

#########################################
Running a simple field cases from scratch
#########################################

In this tutorial, we present how to run a field case by incorporating iteratively different
physics, from a simple flow problem in a continuous media to a coupled poromechanics one
in a fractured media.

We will detailed the way the input XML files need to be written in order to run those
simulations.

In GEOSX, all the simulations parameters have to appear in a XML file, between the XML tags

.. code-block:: xml
  <?xml version="1.0" ?>
  <Problem>
   ...
  </Problem>


*****************
Domain definition
*****************

We consider the following mesh that will be the support for the simulations to come.

.. image:: mesh.png
   :width: 400px

This mesh is composed of three different continuous regions.

  - Top region (overburden, elementary tag = 1)
  - Middle region (reservoir layer, elementary tag = 2)
  - Bottom region (underburden, elementary tag = 3)

The mesh is defined in the GMSH file format (see :ref:`Meshes` for more information on
the supported mesh file format). Each tetrahedra has an elementary tag associated.

.. note::
  GMSH file format is numbering its tags starting from 1. In GEOSX, the numbering
  of the ``CellElementRegion`` starts from 0. As a consequence, when the regions
  will be defined in GEOSX, we have to substract 1 from the GMSH tag to refer
  to the same region. The next sections of this tutorial will illustrate this.

Importing the mesh in GEOSX
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The mesh has to be defined using the Mesh tag (see :ref:`Meshes`). For this problem, we
will use the ``PAMELAMeshGenerator`` to load the mesh (see :ref:`ImportingExternalMesh`).
The syntax is very simple : only the file path (that is relative to the location of the
XML file) and the name of the mesh (that can modified to whatever the user want) need to
be set.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_MESH -->
  :end-before: <!-- SPHINX_FIELD_CASE_MESH_END -->

************************
Running flow simulations
************************

In this part, we will present how to run simple flow simulations to demonstrate the
GEOSX capabilities to handle reservoir simulations.
To ease the demonstration, we make the assumption that the overburden and the underburden are impermeable, the flow only happen in the reservoir. There is two ways to achieve that.

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

We choose the latest in order to keep the burdens for the visualization.

.. note::
  It is possible to define three ``CellElementRegion`` by separating the overburden from
  the underburden

.. note::
  The material list here was set for a single phase flow problem. This list is subject
  to change if the problem is not a single phase flow problem.

Single Phase Flow
^^^^^^^^^^^^^^^^^

We first simulate a single phase flow in the reservoir layer. We define the two materials
involved in the simulation : the only fluid (water) and the solid (rock).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_FIELD_CASE_CONSTITUTIVE_END -->

The constitutive parameters such as the density, the viscosity, the compressibility etc. can
be modified here using the International System of Units. 

.. note::
  To consider an incompressible fluid, the user has to set the compressibility to 0.

Once the materials are defined, we define the numerical method that will be used in the solver.
We propose to use the classical two-point flux approximation scheme.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseTutorial1.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_NUMERICAL -->
  :end-before: <!-- SPHINX_FIELD_CASE_NUMERICAL_END -->

The ``TwoPointFluxApproximation`` node has to contain the primary field to be solved in
``fieldName``. For a flow problem it is the pressure.

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
