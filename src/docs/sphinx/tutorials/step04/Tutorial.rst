.. _TutorialDeadOilBottomLayersSPE10:

####################################################################
Tutorial 4: Multiphase flow with wells
####################################################################

**Context**

In this tutorial, we set up a multiphase, multicomponent test case (see :ref:`CompositionalMultiphaseFlow`).
The permeability field corresponds to the three bottom layers (layers 83, 84, and 85) of the SPE10 test case.
The thermodynamic behavior of the fluid mixture is specified using a Dead-Oil model.
Injection and production are performed using multi-segmented wells.

**Objectives**

At the end of this tutorial you will know:

  - how to import an external mesh with embedded geological properties (porosity and permeability) in the Eclipse format (``.grdecl``),
  - how to set up a multiphase, multicomponent simulation,
  - how to couple reservoir flow with wells.

**Input file**

This tutorial is based on the XML file located at

.. code-block:: console

  src/coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml

The mesh file used in this tutorial is not stored in the main GEOSX repository.
To run the test case specified in the XML file, you must first download the GEOSXDATA repository.  
The XML file that we are going to describe assumes that the GEOSXDATA repository has been
cloned in the same folder as the GEOSX repository.

.. note::
        GEOSXDATA is a separate repository in which we store large mesh files in order to keep the main GEOSX repository lightweight.


------------------------------------
GEOSX input file
------------------------------------

The XML file considered here follows the typical structure of the GEOSX input files:

 #. :ref:`Solver <Solver_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`Mesh <Mesh_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`Geometry <Geometry_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`Events <Events_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`NumericalMethods <NumericalMethods_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`ElementRegions <ElementRegions_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`Constitutive <Constitutive_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`FieldSpecifications <FieldSpecifications_tag_dead_oil_bottom_layers_spe10>`
 #. :ref:`Outputs <Outputs_tag_dead_oil_bottom_layers_spe10>`

.. _Solver_tag_dead_oil_bottom_layers_spe10:

Solvers: coupling reservoir flow with wells
-------------------------------------------

In GEOSX, the simulation of reservoir flow with wells is set up by combining three solvers
listed and parameterized in the **Solvers** XML block of the input file.
We introduce separately a flow solver and a well solver acting on different regions of the
domain---respectively, the reservoir region and the well regions.
To drive the simulation and bind these single-physics solvers, we also specify a *coupling solver*
between the reservoir flow solver and the well solver.
This coupling of single-physics solvers is the generic approach used in GEOSX to
define multiphysics problems.
It is illustrated in :ref:`TutorialPoroelasticity` for a poroelastic test case. 

The three solvers employed in this tutorial are:

 - the single-physics reservoir flow solver, a solver of type **CompositionalMultiphaseFlow** named ``compositionalMultiphaseFlow`` (more information on this solver at :ref:`CompositionalMultiphaseFlow`),
 - the single-physics well solver, a solver of type **CompositionalMultiphaseWell** named ``compositionalMultiphaseWell`` (more information on this solver at :ref:`CompositionalMultiphaseWell`),
 - the coupling solver that binds the two single-physics solvers above, an object of type **CompositionalMultiphaseReservoir** named ``coupledFlowAndWells``.

Let us have a closer look at the **Solvers** XML block displayed below.
Each solver has a name that can be chosen by the user and is not imposed by GEOSX.
These names are used here to point the coupling solver to the single-physics solvers
using the attributes ``flowSolverName`` and ``wellSolverName``. The name of the coupling
solver is also used in the **Events** XML block to trigger the application of the solver.
The coupling solver defines all the target regions on which the single-physics solvers
are applied, namely the reservoir region (named ``reservoir`` here) and one region for each well.

The simulation is fully coupled and driven by the coupled solver. Therefore, the time stepping
information (here, ``initialDt``, but there may be other parameters used to fine-tune the time
stepping strategy), the nonlinear solver parameters, and the linear solver parameters must be
specified at the level of the coupling solver. There is no need to specify these parameters at
the level of the single-physics solvers. 
Note that it is worth repeating the ``logLevel=1`` parameter at the level of the well solver
to make sure that a notification is issued when the well control is switched (from rate control
to BHP control, for instance).

Note that the same fluid and relative permeability
models (set with the ``fluidNames`` and ``relPermNames`` attributes) must be used in the two
single-physics solvers.
GEOSX will throw an error and terminate the simulation if it is not the case.

Take note of the specification of well constraints and controls in the single-physics
well solver with ``control``, ``targetBHP``, ``targetRate``
and ``injectionStream`` (for the composition of the multiphase injection fluid).

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_SOLVERS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_SOLVERS_END -->

.. _Mesh_tag_dead_oil_bottom_layers_spe10:

Specifying a reservoir mesh and defining the geometry of the wells
------------------------------------------------------------------

In the presence of wells, the **Mesh** block of the XML input file includes two parts:

 - a sub-block **PAMELAMeshGenerator** defining the reservoir mesh (see :ref:`TutorialSinglePhaseFlowExternalMesh` for more on this),
 - a collection of sub-blocks **InternalWell** defining the geometry of the wells.

In this tutorial, the reservoir mesh is imported from an Eclipse ``.grdecl`` mesh file that
describes the mesh. It also contains the value of the three components of the permeability
(in the x, y, and z directions) and the value of the porosity for each cell.
The import is requested in the **PAMELAMeshGenerator** XML sub-block. The mesh description
must be done in meters, and the permeability field must be specified in square meters (not in Darcy or milliDarcy).
More information about the mesh importer can be found in :ref:`Meshes`.

Each well is defined internally (i.e., not imported from a file) in a separate **InternalWell**
XML sub-block. An **InternalWell** sub-block must point to the reservoir mesh that the well perforates
using the attribute ``meshName``, to the region corresponding to this well using the attribute
``wellRegionName``, and to the control of this well using the attribute ``wellControl``.

In this tutorial, the five wells have the same structure, with one vertical segment
discretized into four well cells.
We define three perforations along the well (one perforation for each layer of the reservoir mesh).
The location of the perforations is found internally using the linear distance along the wellbore
from the top of the well, specified by the attribute ``distanceFromHead``.
It is the responsibility of the user to make sure that there is a perforation in the bottom cell
of the well mesh otherwise an error will be thrown and the simulation will terminate. The well
geometry must be specified in meters.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_MESH -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_MESH_END -->

.. _Geometry_tag_dead_oil_bottom_layers_spe10:

Geometry tag
-----------------

The **Geometry** XML block was used in the previous tutorials to specify boundary conditions.
Since we use wells and assume no-flow boundary conditions in this tutorial, the **Geometry**
block is not needed.

.. _Events_tag_dead_oil_bottom_layers_spe10:

Specifying events
------------------------

In the **Events** XML block of this tutorial, we specify three types of **PeriodicEvents**
serving different purposes: solver application, result output, and restart file generation.

The periodic event named ``solverApplications`` triggers the application of the solvers
on their target regions. 
For a coupled simulation, this event must point to the coupling solver by name.
The name of the coupling solver is ``coupledFlowAndWells`` and was defined in the **Solvers** block.
The time step is initialized using the ``initialDt`` attribute of the coupling solver.
Then, if the solver converges in more than a certain number of nonlinear iterations (by default, 40% of the
maximum number of nonlinear iterations), the time step will be increased until it reaches the maximum
time step size specified with ``maxEventDt``. 
If the time step fails, the time step will be cut. The parameters defining the time stepping strategy
can be finely tuned by the user in the coupling solver block.
Note that all times are in seconds.

The output event forces GEOSX to write out the results at the frequency specified by the attribute
``timeFrequency``.
Here, we choose to output the results using the VTK format (see :ref:`TutorialSinglePhaseFlowExternalMesh`
for a tutorial that uses the Silo output file format).
Using ``targetExactTimestep=1`` in this XML block forces GEOSX to adapt the time stepping to
ensure that an output is generated exactly at the time frequency requested by the user.
In the ``target`` attribute, we must use the name defined in the **VTK** XML tag
inside the **Output** XML section, as documented at the end of this tutorial (here, ``vtkOutput``).

The restart event instructs GEOSX to write out one or multiple restart files at set times in the simulation.
Restart files are used to restart a simulation from a specific point in time.
These files contain all the necessary information to restoring all internal
variables to their exact state at the chosen time, and continue the simulation from here on.
Here, the ``target`` attribute must contain the name defined in the **Restart** XML sub-block 
of the **Output** XML block (here, ``restartOutput``).


More information about events can be found at :ref:`EventManager`.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_EVENTS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_EVENTS_END -->





.. _NumericalMethods_tag_dead_oil_bottom_layers_spe10:

Defining Numerical Methods
----------------------------------

In the **NumericalMethods** XML block, we select a two-point flux approximation (TPFA) finite-volume scheme to
discretize the governing equations on the reservoir mesh.
TPFA is currently the only numerical scheme that can be used with a flow solver of type
**CompositionalMultiphaseFlow**.
There is no numerical scheme to specify for the well mesh.


.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_NUMERICAL_METHODS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_NUMERICAL_METHODS_END -->


.. _ElementRegions_tag_dead_oil_bottom_layers_spe10:

Defining reservoir and well regions
-----------------------------------

In the presence of wells, it is required to specify at least two types of element regions in
the **ElementRegions** XML block: **CellElementRegion** and **WellElementRegion**.

Here, we define a **CellElementRegion** named ``reservoir`` corresponding to the
reservoir mesh.
The attribute ``cellBlocks`` must be set to ``DEFAULT_HEX`` to point this element region
to the hexahedral mesh corresponding to the bottom layers of SPE10 that was
imported with the **PAMELAMeshGenerator**.
Note that ``DEFAULT_HEX`` is a name internally defined by the mesh importer to denote the only
hexahedral region of the reservoir mesh.
We refer to :ref:`TutorialSinglePhaseFlowExternalMesh` for a discussion on
hexahedral meshes in GEOSX.

.. note::
        If you use a name that is not ``DEFAULT_HEX`` for this attribute, GEOSX will throw an error at the beginning of the simulation.


The **CellElementRegion** must also point to the constitutive models that are used to update
the dynamic rock and fluid properties in the cells of the reservoir mesh.
The names ``fluid``, ``rock``, and ``relperm`` used for this in the ``materialList``
correspond to the attribute ``name`` of the **Constitutive** block.

We also define five **WellElementRegions** corresponding to the five wells.
These regions point to the well meshes defined in the **Mesh** XML block
and to the constitutive models introduced in the **Constitutive** block.
As before, this is done using the names chosen by the user when these blocks
are defined.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_ELEMENT_REGIONS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_ELEMENT_REGIONS_END -->



.. _Constitutive_tag_dead_oil_bottom_layers_spe10:

Defining material properties with constitutive laws
---------------------------------------------------------------------

For a simulation performed with the **CompositionalMultiphaseFlow** 
physics solver, at least three types of constitutive models must be
specified in the **Constitutive** XML block:

- a fluid model describing the thermodynamics behavior of the fluid mixture,
- a relative permeability model,
- a rock compressibility model.

All these models use SI units exclusively.
A capillary pressure model can also be specified in this block but is omitted here for simplicity.
  
Here, we introduce a fluid model  describing a simplified mixture thermodynamic behavior.
Specifically, we use a Dead-Oil model by placing the XML tag **BlackOilFluid** with the attribute
``fluidType=DeadOil``. Other fluid models can be used with the **CompositionalMultiphaseFlow**
solver, as explained in :ref:`FluidModels`.

With the tag **BrooksCoreyRelativePermeability**, we define a relative permeability model.
A list of available relative permeability models can be found at
:ref:`RelativePermeabilityModels`.

The properties are chosen to match those of the original SPE10 test case.
Note that the current fluid model implemented in GEOSX only supports unit
viscosities (i.e., 0.001 Pa.s). Therefore, we set the end point of the oil
relative permeability to 0.1 to preserve the mobility ratio of the SPE10 test case.
     
.. note::
        The names and order of the phases listed for the attribute ``phaseNames`` must be identical in the fluid model (here, **BlackOilFluid**) and the relative permeability model (here, **BrooksCoreyRelativePermeability**). Otherwise, GEOSX will throw an error and terminate. 

We also introduce a model to define the rock compressibility. This step is similar
to what is described in the previous tutorials (see for instance
:ref:`TutorialSinglePhaseFlowWithInternalMesh`).

We remind the reader that the attribute ``name`` of the constitutive models defined here
must be used in the **ElementRegions** and **Solvers** XML blocks to point the element
regions and the physics solvers to their respective constitutive models.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_CONSTITUTIVE_END -->




.. _FieldSpecifications_tag_dead_oil_bottom_layers_spe10:

Defining properties with the FieldSpecifications
---------------------------------------------------------------------

In the **FieldSpecifications** section, we define the initial conditions as well
as the geological properties not read from the mesh file, such as the porosity field in this tutorial.
All this is done using SI units.
Here, we focus on the specification of the initial conditions for
a simulation performed with the **CompositionalMultiphaseFlow** solver.
We refer to :ref:`TutorialSinglePhaseFlowWithInternalMesh` for a more general
discussion on the **FieldSpecification** XML blocks.

For a simulation performed with the **CompositionalMultiphaseFlow** solver,
we have to set the initial pressure as well as the initial global component
fractions (in this case, the oil, gas, and water component fractions).
The ``component`` attribute of the **FieldSpecification** XML block must use the
order in which the ``phaseNames`` have been defined in the **BlackOilFluid**
XML block. In other words, ``component=0`` is used to initialize the oil
global component fraction, ``component=1`` is used to initialize the gas global
component fraction, and ``component=2`` is used to initialize the water global
component fraction, because we previously set ``phaseNames="{oil, gas, water}"``
in the **BlackOilFluid** XML block. Since the SPE10 test case only involves
two-phase flow, we set the initial component fraction of gas to zero.

There is no initialization to perform in the wells since the well
properties are initialized internally using the reservoir initial conditions.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_FIELD_SPECS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_FIELD_SPECS_END -->




.. _Outputs_tag_dead_oil_bottom_layers_spe10:

Specifying the output formats
----------------------------------

In this section, we request an output of the results in VTK format and an output of the restart file.
Note that the names defined here must match the names used in the **Events** XML block to define the output frequency.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/fluidFlow/benchmarks/SPE10/dead_oil_spe10_layers_83_84_85.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_OUTPUT -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_BOTTOM_SPE10_OUTPUT_END -->


All elements are now in place to run GEOSX.


------------------------------------
Running GEOSX
------------------------------------

The first few lines appearing to the console are indicating that the XML elements are read and registered correctly:

.. code-block:: console

  Adding Solver of type CompositionalMultiphaseReservoir, named coupledFlowAndWells
  Adding Solver of type CompositionalMultiphaseFlow, named compositionalMultiphaseFlow
  Adding Solver of type CompositionalMultiphaseWell, named compositionalMultiphaseWell
  Adding Mesh: PAMELAMeshGenerator, mesh
  Adding Mesh: InternalWell, wellProducer1
  Adding Mesh: InternalWell, wellProducer2
  Adding Mesh: InternalWell, wellProducer3
  Adding Mesh: InternalWell, wellProducer4
  Adding Mesh: InternalWell, wellInjector1
  Adding Event: PeriodicEvent, solverApplications
  Adding Event: PeriodicEvent, vtk
  Adding Event: PeriodicEvent, restarts
  Adding Output: VTK, vtkOutput
  Adding Output: Restart, restartOutput
  Adding Object CellElementRegion named reservoir from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion1 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion2 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion3 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion4 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion5 from ObjectManager::Catalog.
		
This is followed by the creation of the 39600 hexahedral cells of the imported mesh:  

.. code-block:: console

  0 >>> **********************************************************************
  0 >>>                          PAMELA Library Import tool                   
  0 >>> **********************************************************************
  0 >>> ECLIPSE GRDECL FORMAT IDENTIFIED
  0 >>> *** Importing Eclipse mesh format file bottom_layers_of_SPE10
  0 >>> ---- Parsing bottom_layers_of_SPE10.grdecl
  0 >>>      o SPECGRID or DIMENS Found
  0 >>>      o COORD Found
  0 >>>      o ZCORN Found
  0 >>>      o PORO Found
  0 >>>      o PERMX Found
  0 >>>      o PERMY Found
  0 >>>      o PERMZ Found
  0 >>> *** Converting into internal mesh format
  0 >>> 39600  total GRDECL hexas
  0 >>> 39600  initially set as active hexas
  0 >>> 0  active but flat hexas (->deactivated)
  0 >>> 0  active but ill-shaped hexas (->deactivated)
  0 >>> 0  duplicated hexas
  0 >>> 39600  will actually be used in the GEOSX mesh
  0 >>> *** Filling mesh with imported properties
  0 >>> *** Creating Polygons from Polyhedra...
  0 >>> 132840 polygons have been created
  0 >>> *** Done
  0 >>> *** Perform partitioning...
  0 >>> TRIVIAL partitioning...
  0 >>> Ghost elements...
  0 >>> Clean mesh...
  0 >>> *** Done...
  0 >>> Clean Adjacency...
  0 >>> *** Done...
			
When ``Running simulation`` is shown, we are done with the case set-up and
the code steps into the execution of the simulation itself:

.. code-block:: console
		
  Time: 0s, dt:1000s, Cycle: 0
    Attempt:  0, NewtonIter:  0
    ( R ) = ( 6.88e+04 ) ; 
    Attempt:  0, NewtonIter:  1
    ( R ) = ( 2.73e+03 ) ; 
    Last LinSolve(iter,res) = (   1, 2.22e-16 ) ; 
    Attempt:  0, NewtonIter:  2
    ( R ) = ( 3.60e+00 ) ; 
    Last LinSolve(iter,res) = (   1, 2.22e-16 ) ; 
    Attempt:  0, NewtonIter:  3
    ( R ) = ( 1.30e-02 ) ; 
    Last LinSolve(iter,res) = (   1, 2.22e-16 ) ; 
    Attempt:  0, NewtonIter:  4
    ( R ) = ( 3.65e-07 ) ; 
    Last LinSolve(iter,res) = (   1, 2.22e-16 ) ; 
  coupledFlowAndWells: Newton solver converged in less than 8 iterations, time-step required will be doubled.

------------------------------------
Visualization of results
------------------------------------

A file compatible with Paraview is produced in this tutorial.
It is found in the output folder, and usually has the extension `.pvd`.
More details about this file format can be found here
`here <https://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format>`_.
We can load this file into Paraview directly and visualize results:

|pic1| |pic2|

.. |pic1| image:: pressure.gif
   :width: 45%

.. |pic2| image:: saturation.gif
   :width: 45%

------------------------------------
To go further
------------------------------------

**Feedback on this tutorial**

This concludes the tutorial on setting up a Dead-Oil simulation in a channelized permeability field.
For any feedback on this tutorial, please submit
a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.

**Next tutorial**

In :ref:`TutorialDeadOilEgg`, we learn how to run a test case based on the Egg model.

**For more details**

  - A complete description of the reservoir flow solver is found here: :ref:`CompositionalMultiphaseFlow`.
  - The well solver is description at :ref:`CompositionalMultiphaseWell`. 
  - The available constitutive models are listed at :ref:`Constitutive`.
