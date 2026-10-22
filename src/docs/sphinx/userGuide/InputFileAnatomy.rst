.. _InputFileAnatomy:

###############################################################################
Anatomy of an Input File
###############################################################################

A GEOS simulation is described entirely by one or more XML input files.
There is no graphical pre-processor and no scripting language to learn: the XML
file declares *what* is simulated (mesh, materials, physics, boundary conditions)
and *when* things happen (time stepping, outputs).
This page explains how such a file is organized and how its parts are connected.
For a hands-on walkthrough of a complete file, see :ref:`TutorialSinglePhaseFlowWithInternalMesh`.


Top-level blocks
================

Every input file has a single root element ``<Problem>``.
Its children are the following blocks. Each block is optional in the schema, but a
typical simulation uses the first nine.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Block
     - Purpose
   * - ``<Mesh>``
     - The computational mesh: generated internally (``InternalMesh``, ``InternalWellbore``)
       or imported from a file (``VTKMesh``). Wells are declared here too. See :ref:`Meshes`.
   * - ``<Geometry>``
     - Simple geometric objects (``Box``, ``Cylinder``, ``Disc``, ``Rectangle``, ``ThickPlane``,
       ``CustomPolarObject``) used to define named *sets* of mesh objects, or fracture planes.
   * - ``<ElementRegions>``
     - Groups mesh cell blocks into regions (``CellElementRegion``) and assigns materials to them.
       Fracture (``SurfaceElementRegion``) and well (``WellElementRegion``) regions are also declared here.
   * - ``<Constitutive>``
     - Material models: fluids, solids, porosity, permeability, relative permeability, etc.
       See :ref:`Constitutive`.
   * - ``<NumericalMethods>``
     - Spatial discretizations: finite-volume flux approximations and finite-element spaces.
       See :ref:`NumericalMethodsManager`.
   * - ``<Solvers>``
     - Physics solvers and their nonlinear/linear solver settings. See :ref:`Solvers`.
   * - ``<FieldSpecifications>``
     - Initial conditions, boundary conditions and source terms. See :ref:`FieldSpecification`.
   * - ``<Events>``
     - The simulation timeline: when solvers, outputs and tasks run, and when the run ends.
       See :ref:`EventManager`.
   * - ``<Outputs>``
     - Output formats (``VTK``, ``Silo``, ``TimeHistory``, ``Restart``, ...). See :ref:`Outputs`.
   * - ``<Functions>``
     - Tables and analytical functions used to make inputs vary in time or space.
       See :ref:`FunctionManager`.
   * - ``<Tasks>``
     - Actions that are not physics solvers: data collection, statistics, initialization,
       constitutive model drivers. See :ref:`TasksManager`.
   * - ``<Included>``
     - Other XML files to merge into this one (see below).
   * - ``<Parameters>``
     - Named values for the XML preprocessor (see below).
   * - ``<ExternalDataSource>``
     - External data sources (e.g. VTK hierarchical data) referenced by the mesh.
   * - ``<ParticleRegions>``
     - Particle regions, used by the Material Point Method solver.
   * - ``<Benchmarks>``
     - Machine-specific settings used by the performance benchmark suite.


How blocks refer to each other
==============================

Every object in the file has a ``name`` attribute. Other objects refer to it by
that name. This is how the blocks above are connected: no block is implicitly
linked to another. The example below is a complete, runnable deck, adapted from
``inputFiles/singlePhaseFlow/3D_10x10x10_compressible_base.xml`` and its ``_smoke.xml`` companion.
Comments mark each cross-reference.

.. code-block:: xml

  <Problem>
    <Solvers>
      <!-- discretization -> NumericalMethods; targetRegions -> ElementRegions -->
      <SinglePhaseFVM name="SinglePhaseFlow"
                      discretization="singlePhaseTPFA"
                      targetRegions="{ mainRegion }">
        <NonlinearSolverParameters newtonTol="1.0e-6" newtonMaxIter="8"/>
        <LinearSolverParameters solverType="gmres" preconditionerType="amg" krylovTol="1.0e-10"/>
      </SinglePhaseFVM>
    </Solvers>

    <Mesh>
      <!-- cellBlockNames defines the cell block "cellBlock" -->
      <InternalMesh name="mesh" elementTypes="{ C3D8 }"
                    xCoords="{ 0, 10 }" yCoords="{ 0, 10 }" zCoords="{ 0, 10 }"
                    nx="{ 10 }" ny="{ 10 }" nz="{ 10 }"
                    cellBlockNames="{ cellBlock }"/>
    </Mesh>

    <Geometry>
      <!-- each geometric object defines a set with the same name -->
      <Box name="source" xMin="{ -0.01, -0.01, -0.01 }" xMax="{ 1.01, 1.01, 1.01 }"/>
      <Box name="sink"   xMin="{ 8.99, 8.99, 8.99 }"    xMax="{ 10.01, 10.01, 10.01 }"/>
    </Geometry>

    <NumericalMethods>
      <FiniteVolume>
        <TwoPointFluxApproximation name="singlePhaseTPFA"/>
      </FiniteVolume>
    </NumericalMethods>

    <ElementRegions>
      <!-- cellBlocks -> Mesh ("*" = all blocks); materialList -> Constitutive -->
      <CellElementRegion name="mainRegion" cellBlocks="{ * }" materialList="{ water, rock }"/>
    </ElementRegions>

    <Constitutive>
      <CompressibleSinglePhaseFluid name="water" defaultDensity="1000" defaultViscosity="0.001"
                                    referencePressure="0.0" compressibility="5e-10" viscosibility="0.0"/>
      <!-- a "rock" is assembled from three sub-models, referenced by name -->
      <CompressibleSolidConstantPermeability name="rock"
                                             solidModelName="nullSolid"
                                             porosityModelName="rockPorosity"
                                             permeabilityModelName="rockPerm"/>
      <NullModel name="nullSolid"/>
      <PressurePorosity name="rockPorosity" defaultReferencePorosity="0.05"
                        referencePressure="0.0" compressibility="1.0e-9"/>
      <ConstantPermeability name="rockPerm" permeabilityComponents="{ 1.0e-12, 1.0e-12, 1.0e-15 }"/>
    </Constitutive>

    <FieldSpecifications>
      <!-- objectPath -> ElementRegions/<region>/<cell block>; setNames -> Geometry or built-in sets -->
      <FieldSpecification name="initialPressure" initialCondition="1" setNames="{ all }"
                          objectPath="ElementRegions/mainRegion/cellBlock"
                          fieldName="pressure" scale="5e6"/>
      <FieldSpecification name="sourceTerm" setNames="{ source }"
                          objectPath="ElementRegions/mainRegion/cellBlock"
                          fieldName="pressure" scale="1e7"/>
      <FieldSpecification name="sinkTerm" setNames="{ sink }"
                          objectPath="ElementRegions/mainRegion/cellBlock"
                          fieldName="pressure" scale="0.0"/>
    </FieldSpecifications>

    <Events maxTime="5000.0">
      <!-- target -> /<Block>/<name> -->
      <PeriodicEvent name="solverApplications" forceDt="20.0" target="/Solvers/SinglePhaseFlow"/>
      <PeriodicEvent name="outputs" timeFrequency="100.0" target="/Outputs/vtkOutput"/>
    </Events>

    <Outputs>
      <VTK name="vtkOutput"/>
    </Outputs>
  </Problem>

The references used most often are summarized below.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Attribute
     - Found in
     - Refers to the ``name`` of
   * - ``cellBlocks``
     - ``CellElementRegion``
     - cell blocks created by the mesh (``cellBlockNames`` of ``InternalMesh``,
       or region attribute values of a ``VTKMesh``). ``{ * }`` selects all of them.
   * - ``materialList``
     - ``CellElementRegion``
     - models in ``<Constitutive>``. Solvers find the fluid, solid, etc. they need
       in the ``materialList`` of their target regions.
   * - ``solidModelName``, ``porosityModelName``, ``permeabilityModelName``
     - composite solid models (e.g. ``CompressibleSolidConstantPermeability``, ``PorousElasticIsotropic``)
     - other models in ``<Constitutive>``.
   * - ``targetRegions``
     - every solver
     - regions in ``<ElementRegions>``.
   * - ``discretization``
     - flow and mechanics solvers
     - an object in ``<NumericalMethods>``.
   * - ``flowSolverName``, ``solidSolverName``, ``wellSolverName``
     - coupled solvers (poromechanics, reservoir + wells, ...)
     - other solvers in ``<Solvers>``.
   * - ``objectPath``
     - field specifications
     - a path in the mesh hierarchy: ``ElementRegions/<region>/<cell block>`` for cell fields,
       ``nodeManager`` for nodal fields, ``faceManager`` for face fields.
   * - ``setNames``
     - field specifications
     - sets of mesh objects: geometric objects from ``<Geometry>``, sets created by the mesh
       generator, or the built-in set ``all``. See :ref:`FieldSpecification`.
   * - ``functionName``
     - field specifications
     - a function in ``<Functions>``.
   * - ``target``
     - events
     - an absolute path to a solver, output or task: ``/Solvers/<name>``, ``/Outputs/<name>``, ``/Tasks/<name>``.

A misspelled name is the most common error in a new input file, and GEOS stops
with an error when a reference cannot be resolved.


Splitting a deck: ``<Included>``
================================

Large decks are usually split into a *base* file (physics, materials, regions) and
smaller *driver* files (mesh resolution, run length) that include it:

.. code-block:: xml

  <Problem>
    <Included>
      <File name="./3D_10x10x10_compressible_base.xml"/>
    </Included>
    <!-- blocks below are merged with the included file -->
  </Problem>

Relative paths are resolved from the location of the including file.
Several input files can also be given on the command line (``geos -i a.xml -i b.xml``);
they are merged in the same way. Merge rules are described in :ref:`InputFiles`.


Parameters and units: ``<Parameters>``
======================================

``<Parameters>`` defines named values that can be reused anywhere in the file as ``$name$``.
These, along with units (e.g. ``"1 [min]"``) and inline arithmetic (e.g. ``"`$t0$+$dt$`"``), are
handled by the Python XML preprocessor, not by GEOS itself. A file using them must be
preprocessed before running; see :ref:`AdvancedXMLFeatures`.


Looking up any element or attribute
===================================

The complete list of XML elements and attributes, with types, defaults and
descriptions, is generated from the code and always up to date:

* :doc:`/docs/sphinx/datastructure/CompleteXMLSchema` lists every element. For example,
  the solver used above is documented at :ref:`XML_SinglePhaseFVM`.
* Your build of GEOS can write its own schema with ``geos -s schema.xsd``. Point your
  XML editor at this file to get validation and auto-completion (see :ref:`InputFiles`).
* ``geos -i input.xml -v`` reads and checks an input file without running the simulation.
