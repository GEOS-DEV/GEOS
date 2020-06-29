.. _TutorialCO2FieldCaseUnstructuredGrid:

##################################################
Tutorial 6: CO2 injection into unstructured grid 
##################################################

**Context**
In this tutorial, we go on with our previous field case (see :ref:`TutorialFieldCase`) adding a CO2 injection well in the highest point of the reservoir. 

**Objectives**

At the end of this tutorial you will know:

 - how to set up a Co2 injection scenario
 - how to add well coupling into the domain 
 - how to run a case using mpi-parallelism

**Input file**

The xml file for this test case is located at :

.. code-block:: console

  src/coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml

It contains the follwing tags:

 #. :ref:`Solver <Solver_tag_co2_field_case>`
 #. :ref:`Mesh <Mesh_tag_co2_field_case>`
 #. :ref:`Geometry <Geometry_tag_co2_field_case>`
 #. :ref:`Events <Events_tag_co2_field_case>`
 #. :ref:`NumericalMethods <NumericalMethods_tag_co2_field_case>`
 #. :ref:`ElementRegions <ElementRegions_tag_co2_field_case>`
 #. :ref:`Constitutive <Constitutive_tag_co2_field_case>`
 #. :ref:`FieldSpecifications <FieldSpecifications_tag_co2_field_case>`
 #. :ref:`Outputs <Outputs_tag_co2_field_case>`

Domain definition
--------------------

We consider the field case mesh as a numerical support to the simulations with a single injection point:

.. image:: mesh_andWell-2.png
   :width: 600px

This mesh contains three continuous regions:

  - a Top region (overburden, elementary tag = 1)
  - a Middle region (reservoir layer, elementary tag = 2)
  - a Bottom region (underburden, elementary tag = 3)

A single injection wellbore will be at the center of the reservoir. The picture shows an example of the 8-core METIS partitioning used to lauch the simulation.


------------------------------------
GEOSX input file
------------------------------------

.. _Solver_tag_co2_field_case:

Defining a solver
-----------------
Let us inspect the **Solver** XML tags.
It consists in 3 blocks *CompositionalMultiphaseFlow*, *CompositionalMultiphaseWell* and *CompositionalMultiphaseReservoir*, which are respectively handling solution from multiphase flow in the reservoir, multiphase flow in the wells and coupling between those two parts.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_SOLVER -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_SOLVER_END -->

In the *CompositionalMultiphaseFlow* (:ref:`CompositionalMultiphaseFlow`), a classical multiphase compositional solver is detailed, including a TPFA discretization, reference to fluid data through *fluidNames*, to solid data through *soldiNames* and to relative permeability model through *relPermNames* attributes.

The *CompositionalMultiphaseWell* (:ref:`CompositionalMultiphaseWell`)  consists in wellbore specifications (see :ref:`TutorialDeadOilBottomLayersSPE10` for detailed tutorial on wells integration). As its reservoir counterpart, it includes references to fluid and relative permeabilities models, but also defines *WellControls* sub-tag, that can specified injector and producer `control` spliting between BHP-controlled or rate-controlled. Alongside with that attribute are the *targetBHP* and *targetRate*, that specify the maximal admissible pressure and rate for the well. The injector specific attribute, *injectionStream*, descibes the composition of the injected mixture.

The coupling section *CompositionalMultiphaseReservoir* describes the binding between those two previous element (see :ref:`TutorialPoroelasticity` for detailed tutorial on coupling physics in GEOSX). In addition to bound to the previously described blocks through *flowSolverName* and *wellSolverName* sub-tags, it contains the *initialDt* starting time-step size value and defined the *NonlinearParameters* and *LinearSolverParameters* that are used to control newton-loop and linear solver behaviors.(see :ref:`LinearSolvers` for a detailed description of linear solvers attributes) 

.. _Mesh_tag_co2_field_case:

Specifying a computational mesh
---------------------------------
The **Mesh** tag is used as in previous tutorials to import field mesh either internally (:ref:`TutorialSinglePhaseFlowWithInternalMesh`) or externally (:ref:`TutorialSinglePhaseFlowExternalMesh`). In the current tutorial, it will also be of paramount importance as it will be where the *InternalWell* multi-segmented wells will be defined. Apart from the `name` identifier attribute and their `wellRegionName` (:ref:`ElementRegions <ElementRegions_tag_co2_field_case>`) and `wellControlsName` (:ref:`Solver <Solver_tag_co2_field_case>`) binding attributes, `polylineNodeCoords` and `polylineSegmentConn` attributes will be used to define path of the wellbore and connection between those nodes. The `numElementPerSegment` is discretizing the wellbore's segments while the `radius` attribute stands for the wellbore radius. (:ref:`TutorialDeadOilBottomLayersSPE10` for details on wellbore use). Once the wellbore are defined and discretized, the place for different *Perforation* can be defined using curvilinear distance from the head of the wellbore `distanceFromHead`.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_MESH -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_MESH_END -->

.. note::
        It is the responsibility of the user to make sure that there is a perforation in the bottom cell of the well mesh otherwise an error will be thrown and the simulation will terminate.


.. _Geometry_tag_co2_field_case:

Geometry tag
----------------

  Here we use the usual "infinity" box to put all the domain in the *all* set.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_GEOMETRY -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_GEOMETRY_END -->


.. _Events_tag_co2_field_case:

Specifying events
------------------------
        
   The solver is applied as a recurent event, whose target is refered as **Solver/coupledFlowAnWells** name-tag. The outputs are periodically written every 11 days and 24 hours, constraining schedule to match exactly this date. The output path to data is specified as a *target* of this *PeriodicEvent*.


.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_EVENTS -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_EVENTS_END -->

An other periodic event is also defined under the name `restarts`. It consists of saved checkpoints every 116 days, whose physical output folder name will be defined under the **Output** tag.

.. _NumericalMethods_tag_co2_field_case:

Defining Numerical Methods
----------------------------------

The ``TwoPointFluxApproximation`` is chosen as our fluid equation discretization. The node should specify
        -the primary field to solve for as ``fieldName``. For a flow problem, this field is the pressure. 
        - the ``boundaryFieldName`` is used to specify boundary object for imposing boundary conditions.
        -the ``coefficientName`` is used during TPFA transmissibilities construction.

      Here we specified ``targetRegions`` as we only solve flow for reservoir.

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_NUMERICAL -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_NUMERICAL_END -->

.. _ElementRegions_tag_co2_field_case:

Defining regions in the mesh
-----------------------------------
As in :ref:`TutorialFieldCase`, the **ElementRegions** tag allows us to split over-and-under burden and reservoir into two different entities. What's new here is the addition of *WellElementRegion* that are not bound to a cellBlock and that contains a list of materials present.


.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_REGION -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_REGION_END -->


.. _Constitutive_tag_co2_field_case:

Defining material properties with constitutive laws
-------------------------------------------------------
Under the **Constitutive** tag, three items can be found:
        - *MultiPhaseMultiComponentFluid* which will allow us to define phases names, component molar weights and caracteristic behaviors such as viscosity and density dependecies with respect to pressure and temperature. 
        - *PoreVolumeCompressibleSolid* which contains all the data needed to model rock compressibility behavior.
        - *BrooksCoreyRelativePermeability* which set the relative permeability for each phase along with its end-point value, residual volume fraction and Corey exponent.
          
.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_CONSTITUTIVE_END -->

One can notice that PVT data, required by *MultiPhaseMultiComponentFluid*, are brought in to be able to model behavior of Co2 in liquid and gas phase with repect to pressure and temperature variations. These *pvtgas.txt* and *pvtliquid.txt* are composed as follow

.. code:: 

        DensityFun SpanWagnerCO2Density 1e6 1.5e7 5e4 94 96 1
        ViscosityFun FenghourCO2Viscosity 1e6 1.5e7 5e4 94 96 

.. code::
  
        DensityFun BrineCO2Density 1e6 1.5e7 5e4 94 96 1 0
        ViscosityFun BrineViscosity 0

The first keyword is an identifier for either density or viscosity model generated in GEOSX at runtime. It is followed by an identifier for the type of the model (see ref:`PVTModels`) before the lower, upper and step increment values for pressure and temperature range. The trailing 0 for BrineCO2Density entry is the salinity of the brine (see :ref:`CO2-EOS`)

.. note::
  The *0* value for *BrineViscosity* indicates that liquid Co2 viscosity is constant with respect to pressure and temperature. 

.. _FieldSpecifications_tag_co2_field_case:

Defining properties with the FieldSpecifications
---------------------------------------------------------------------
As in previous tutorials, **FieldSpecifications** tag is the place where to declare all the scoped fields such as directional permeability, reference porosity, initial pressure and compositions. These can be homogeneously and permanantly fixed or use *TableFunctions* via `functionName` attributes to activate  modifiers and to set them heterogeneously. (see :ref:`TutorialFieldCase` for details)

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_FIELD -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_FIELD_END -->

.. _Outputs_tag_co2_field_case:

Specifying the output formats
----------------------------------

The **Outputs** XML tag is used to trigger the writing of visualization and restart files.
Here, we write files in a format natively readable by Paraview under the tag *VTK*. 

.. literalinclude:: ../../../../coreComponents/physicsSolvers/multiphysics/integratedTests/FieldCaseCo2InjTutorial.xml
  :language: xml
  :start-after: <!-- SPHINX_FIELD_CASE_Co2_OUTPUT -->
  :end-before: <!-- SPHINX_FIELD_CASE_Co2_OUTPUT_END -->

A *Restart* tag can also be specified. In conjunction with a *PeriodicEvent*, it allows resuming computation from a checkpoint in time. 

.. 
        ------------------------------------------------
        Using Functions to specify dependent properties
        ------------------------------------------------
        The *TableFunctions* are used to map heterogeneous properties (here permeabiity and porosity) onto the mesh via the reading of regular grid coordinate in *xlin.geos*, *ylin.geos* and *zlin.geos* and a ND-table at *voxelFile* adress containing heterogenous data. (see :ref:`TutorialFieldCase` for more details on those) 


------------------------------------
Runnning GEOSX
------------------------------------

The simulation can be launched with on 8-cores using MPI-parallelism: 

.. code-block:: console

  mpirun -np 8 geosx -i FieldCaseCo2InjTutorial.xml -x 2 -y 2 -z 2

Then, the multicore loading will rank-prefix the usual GEOSX mesh pre-treatment output as

.. code-block:: console

        1 >>> **********************************************************************
        1 >>>                          PAMELA Library Import tool                   
        1 >>> **********************************************************************
        1 >>> GMSH FORMAT IDENTIFIED
        1 >>> *** Importing Gmsh mesh format...
        3 >>> **********************************************************************
        6 >>> **********************************************************************
        6 >>>                          PAMELA Library Import tool                   
        6 >>> **********************************************************************
        6 >>> GMSH FORMAT IDENTIFIED
        3 >>>                          PAMELA Library Import tool                   
        3 >>> **********************************************************************
        3 >>> GMSH FORMAT IDENTIFIED
        6 >>> *** Importing Gmsh mesh format...
        3 >>> *** Importing Gmsh mesh format...
        7 >>> **********************************************************************
        7 >>>                          PAMELA Library Import tool                   
        7 >>> **********************************************************************
        7 >>> GMSH FORMAT IDENTIFIED
        7 >>> *** Importing Gmsh mesh format...
        0 >>> **********************************************************************
        0 >>>                          PAMELA Library Import tool                   
        0 >>> **********************************************************************
        0 >>> GMSH FORMAT IDENTIFIED
        0 >>> *** Importing Gmsh mesh format...
        5 >>> **********************************************************************
        5 >>>                          PAMELA Library Import tool                   
        5 >>> **********************************************************************
        5 >>> GMSH FORMAT IDENTIFIED
        5 >>> *** Importing Gmsh mesh format...
        4 >>> **********************************************************************
        4 >>>                          PAMELA Library Import tool                   
        4 >>> **********************************************************************
        4 >>> GMSH FORMAT IDENTIFIED
        4 >>> *** Importing Gmsh mesh format...
        2 >>> **********************************************************************
        2 >>>                          PAMELA Library Import tool                   
        2 >>> **********************************************************************
        2 >>> GMSH FORMAT IDENTIFIED
        2 >>> *** Importing Gmsh mesh format...
        7 >>> Reading nodes...
        4 >>> Reading nodes...
        1 >>> Reading nodes...
        5 >>> Reading nodes...
        0 >>> Reading nodes...
        6 >>> Reading nodes...
        3 >>> Reading nodes...
        2 >>> Reading nodes...
        6 >>> Done0
        6 >>> Reading elements...
        5 >>> Done0
        5 >>> Reading elements...
        4 >>> Done0
        4 >>> Reading elements...
        3 >>> Done0
        3 >>> Reading elements...
        1 >>> Done0
        1 >>> Reading elements...
        7 >>> Done0
        7 >>> Reading elements...
        0 >>> Done0
        0 >>> Reading elements...
        2 >>> Done0

A restart from a checkpoint file `FieldCaseCo2InjTutorial_restart_000000015.root` is always available thanks to the following command line :

.. code-block:: console

  mpirun -np 8 geosx -i FieldCaseCo2InjTutorial.xml -x 2 -y 2 -z 2 -r FieldCaseCo2InjTutorial_restart_000000015  

The output then shows the loading of HDF5 restart files by each core. 

.. code-block:: console

        Loading restart file /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015
        Rank 0: rankFilePattern = /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_%07d.hdf5
        Rank 0: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000000.hdf5
        Rank 4: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000004.hdf5
        Rank 5: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000005.hdf5
        Rank 7: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000007.hdf5
        Rank 1: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000001.hdf5
        Rank 3: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000003.hdf5
        Rank 2: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000002.hdf5
        Rank 6: Reading in restart file at /work/simu/geosx/test/FieldCaseCo2InjTutorial_restart_000000015/rank_0000006.hdf5


------------------------------------
Visualization of results
------------------------------------

Post-treating under Paraview, we can isolate the *Reservoir* block and focus on planes othogonal to the injection wellbore.  

.. image:: fcCo2_saturation-1.png
   :width: 600px

Closing up to the wellbore, we can see the reservoir filling up with Co2,

.. image:: fcCo2_sat.gif
   :width: 600px

Inspecting the pressure value along the same plane slices, we can check pressure rise at the injector

.. image:: fcCo2_pres.gif
   :width: 600px

------------------------------------
To go further
------------------------------------

**Feedback on this tutorial**

This concludes the Co2 injection field case tutorial.
For any feedback on this tutorial, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.

**Next tutorial**

In the next tutorial :ref:`TutorialElasticity`, we switch to mechanics and learn how to run a bending case to get familiar with mechanical problems in GEOSX.

**For more details**

  - A complete description of the reservoir flow solver is found here: :ref:`CompositionalMultiphaseFlow`.
  - The well solver is described at :ref:`CompositionalMultiphaseWell`. 
  - The available fluid constitutive models are listed at :ref:`FluidModels`.

