
.. _TutorialDeadOilEgg:

########################################################
Multiphase Flow with Wells
########################################################

**Context**

In this example, we build on the concepts presented in :ref:`TutorialDeadOilBottomLayersSPE10`
to show how to set up a multiphase water injection problem with wells in
the three-dimensional `Egg model <https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/gdj3.21>`_.
The twelve wells (four producers and eight injectors) are placed according to the description
of the original test case.

**Objectives**

In this example, we re-use many GEOSX features already presented in
:ref:`TutorialDeadOilBottomLayersSPE10`, but we now focus on:

- how to import an external mesh with embedded geological properties (permeability) in the GMSH format (``.msh``),
- how to set up the wells.

**Input file**

This example is based on the XML file located at

.. code-block:: console

  ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml

The mesh file corresponding to the Egg model is stored in the GEOSXDATA repository.
Therefore, you must first download the GEOSXDATA repository in the same folder
as the GEOSX repository to run this test case.

.. note::
        GEOSXDATA is a separate repository in which we store large mesh files in order to keep the main GEOSX repository lightweight.
   
The XML file considered here follows the typical structure of the GEOSX input files:

 #. :ref:`Solver <Solver_tag_dead_oil_egg_model>`
 #. :ref:`Mesh <Mesh_tag_dead_oil_egg_model>`
 #. :ref:`Events <Events_tag_dead_oil_egg_model>`
 #. :ref:`NumericalMethods <NumericalMethods_tag_dead_oil_egg_model>`
 #. :ref:`ElementRegions <ElementRegions_tag_dead_oil_egg_model>`
 #. :ref:`Constitutive <Constitutive_tag_dead_oil_egg_model>`
 #. :ref:`FieldSpecifications <FieldSpecifications_tag_dead_oil_egg_model>`
 #. :ref:`Outputs <Outputs_tag_dead_oil_egg_model>`
 #. :ref:`Tasks <Tasks_tag_dead_oil_egg_model>`    

.. _Solver_tag_dead_oil_egg_model:

-----------------------------------------
Coupling the flow solver with wells
-----------------------------------------

In GEOSX, the simulation of reservoir flow with wells is set up by combining three solvers
listed and parameterized in the **Solvers** XML block of the input file.
We introduce separately a flow solver and a well solver acting on different regions of the
domain---respectively, the reservoir region and the well regions.
To drive the simulation and bind these single-physics solvers, we also specify a *coupling solver*
between the reservoir flow solver and the well solver.
This coupling of single-physics solvers is the generic approach used in GEOSX to
define multiphysics problems.
It is illustrated in :ref:`TutorialPoroelasticity` for a poroelastic test case. 

The three solvers employed in this example are:

 - the single-physics reservoir flow solver, a solver of type **CompositionalMultiphaseFVM** named ``compositionalMultiphaseFlow`` (more information on this solver at :ref:`CompositionalMultiphaseFlow`),
 - the single-physics well solver, a solver of type **CompositionalMultiphaseWell** named ``compositionalMultiphaseWell`` (more information on this solver at :ref:`CompositionalMultiphaseWell`),
 - the coupling solver that binds the two single-physics solvers above, an object of type **CompositionalMultiphaseReservoir** named ``coupledFlowAndWells``.

The **Solvers** XML block is shown below.
The coupling solver points to the two single-physics solvers using the attributes
``flowSolverName`` and ``wellSolverName``.
These names can be chosen by the user and are not imposed by GEOSX.
The flow solver is applied to the reservoir and the well solver is applied to the wells,
as specified by their respective ``targetRegions`` attributes.

The simulation is fully coupled and driven by the coupled solver. Therefore, the time stepping
information (here, ``initialDt``, but there may be other parameters used to fine-tune the time
stepping strategy), the nonlinear solver parameters, and the linear solver parameters must be
specified at the level of the coupling solver.
There is no need to specify these parameters at the level of the single-physics solvers.
Any solver information specified in the single-physics XML blocks will not be taken into account.

.. note::
        It is worth repeating the ``logLevel="1"`` parameter at the level of the well solver to make sure that a notification is issued when the well control is switched (from rate control to BHP control, for instance).

Here, we instruct GEOSX to perform at most ``newtonMaxIter = "10"`` Newton iterations. 
GEOSX will adjust the time step size as follows:

- if the Newton solver converges in ``dtIncIterLimit x newtonMaxIter = 5`` iterations or fewer, GEOSX will double the time step size for the next time step,
- if the Newton solver converges in ``dtCutIterLimit x newtonMaxIter = 8`` iterations or more, GEOSX will reduce the time step size for the next time step by a factor ``timestepCutFactor = 0.1``,
- if the Newton solver fails to converge in ``newtonMaxIter = 10``, GEOSX will cut the time step size by a factor ``timestepCutFactor = 0.1`` and restart from the previous converged time step.

The maximum number of time step cuts is specified by the attribute ``maxTimeStepCuts``.
Note that a backtracking line search can be activated by setting the attribute ``lineSearchAction`` to ``Attempt`` or ``Require``.
If ``lineSearchAction = "Attempt"``, we accept the nonlinear iteration even if the line search does not reduce the residual norm.
If ``lineSearchAction = "Require"``, we cut the time step if the line search does not reduce the residual norm. 

.. note::
   To use the linear solver options of this example, you need to ensure that GEOSX is configured to use the Hypre linear solver package.


.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_SOLVERS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_SOLVERS_END -->

	       

.. _Mesh_tag_dead_oil_egg_model:

-------------------------------------------------
Mesh definition and well geometry
-------------------------------------------------

In the presence of wells, the **Mesh** block of the XML input file includes two parts:

 - a sub-block **PAMELAMesh** defining the reservoir mesh (see :ref:`TutorialSinglePhaseFlowExternalMesh` for more on this),
 - a collection of sub-blocks **InternalWell** defining the geometry of the wells.

The reservoir mesh is imported from a ``.msh`` file that contains the mesh geometry
and also includes the permeability values in the x, y, and z directions.
These quantities must be specified using the metric unit system, i.e., in meters
for the well geometry and square meters for the permeability field.
We note that the mesh file only contains the active cells, so there is no keyword
needed in the XML file  to define them.

Each well is defined internally (i.e., not imported from a file) in a separate **InternalWell**
XML sub-block. An **InternalWell** sub-block must point to the reservoir mesh that the well perforates
using the attribute ``meshName``, to the region corresponding to this well using the attribute
``wellRegionName``, and to the control of this well using the attribute ``wellControl``.
Each block **InternalWell** must point to the reservoir mesh
(using the attribute ``meshName``), the corresponding well region (using
the attribute ``wellRegionName``), and the corresponding well control
(using the attribute ``wellControlName``).

Each well is defined using a vertical polyline going through the seven layers of the
mesh, with a perforation in each layer.
The well placement implemented here follows the pattern of the original test case.
The well geometry must be specified in meters.

The location of the perforations is found internally using the linear distance along the wellbore
from the top of the well, specified by the attribute ``distanceFromHead``.
It is the responsibility of the user to make sure that there is a perforation in the bottom cell
of the well mesh otherwise an error will be thrown and the simulation will terminate.
For each perforation, the well transmissibility factors employed to compute the perforation rates are calculated
internally using the Peaceman formulation.


.. image:: egg_model.png
   :width: 400px
   :align: center          

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_benchmark.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_MESH -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_MESH_END -->


.. _Events_tag_dead_oil_egg_model:

------------------------
Events
------------------------


In the **Events** XML block, we specify four types of **PeriodicEvents**.

The periodic event named ``solverApplications`` notifies GEOSX that the
coupled solver ``coupledFlowAndWells`` has to be applied to its target
regions (here, reservoir and wells) at every time step.
The time stepping strategy has been fully defined in the **CompositionalMultiphaseReservoir**
coupling block using the ``initialDt`` attribute and the **NonlinearSolverParameters**
nested block.

We also define an output event instructing GEOSX to write out ``.vtk`` files at the time frequency specified
by the attribute ``timeFrequency``.
Here, we choose to output the results using the VTK format (see :ref:`TutorialSinglePhaseFlowExternalMesh`
for a example that uses the Silo output file format).
The ``target`` attribute must point to the **VTK** sub-block of the **Outputs**
block (defined at the end of the XML file) by name (here, ``vtkOutput``).

We define the events involved in the collection and output of the well production rates following the procedure defined in :ref:`TasksManager`.
The time history collection events trigger the collection of the well rates at the desired frequency, while the time history output events trigger the output of the HDF5 files containing the time series.
These events point by name to the corresponding blocks of the **Tasks** and **Outputs** XML blocks, respectively. Here, these names are ``wellRateCollection1`` and ``timeHistoryOutput1``.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_EVENTS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_EVENTS_END -->


.. _NumericalMethods_tag_dead_oil_egg_model:

----------------------------------
Numerical methods
----------------------------------

In the ``NumericalMethods`` XML block, we instruct GEOSX to use a TPFA finite-volume
numerical scheme.
This part is similar to the corresponding section of :ref:`TutorialDeadOilBottomLayersSPE10`, and has been adapted to match the specifications of the Egg model.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_NUMERICAL_METHODS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_NUMERICAL_METHODS_END -->


.. _ElementRegions_tag_dead_oil_egg_model:

-----------------------------------
Reservoir and well regions
-----------------------------------

In this section of the input file, we follow the procedure already described in
:ref:`TutorialDeadOilBottomLayersSPE10` for the definition of the reservoir region with multiphase constitutive models.

We associate a **CellElementRegion** named ``reservoir`` to the reservoir mesh.
Since we have imported a mesh with one region consisting of hexahedral cells, we
must set the attribute ``cellBlocks`` to ``DEFAULT_HEX``.

.. note::
        If you use a name that is not ``DEFAULT_HEX`` for this attribute, GEOSX will throw an error at the beginning of the simulation.

We also associate a **WellElementRegion** to each well. As the **CellElementRegion**,
it contains a ``materialList`` that must point (by name) to the constitutive models
defined in the **Constitutive** XML block.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_ELEMENT_REGIONS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_ELEMENT_REGIONS_END -->


.. _Constitutive_tag_dead_oil_egg_model:

-------------------
Constitutive models
-------------------

The **CompositionalMultiphaseFVM** physics solver relies on at least four types of constitutive
models listed in the **Constitutive** XML block:

- a fluid model describing the thermodynamics behavior of the fluid mixture,
- a relative permeability model,
- a rock permeability model,
- a rock porosity model.

All the parameters must be provided using the SI unit system.

This part is identical to that of :ref:`TutorialDeadOilBottomLayersSPE10`. 

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_CONSTITUTIVE_END -->


.. _FieldSpecifications_tag_dead_oil_egg_model:

-----------------------
Initial conditions
-----------------------

We are ready to specify the reservoir initial conditions of the problem in the **FieldSpecifications**
XML block.
The well variables do not have to be initialized here since they will be defined internally.

The formulation of the **CompositionalMultiphaseFVM** physics solver (documented
at :ref:`CompositionalMultiphaseFlow`) requires the definition of the initial pressure field
and initial global component fractions.
We define here a uniform pressure field that does not satisfy the hydrostatic equilibrium,
but a hydrostatic initialization of the pressure field is possible using :ref:`FunctionManager`:.
For the initialization of the global component fractions, we remind the user that their ``component``
attribute (here, 0 or 1) is used to point to a specific entry of the ``phaseNames`` attribute
in the **DeadOilFluid** block. 

Note that we also define the uniform porosity field here since it is not included in the mesh file
imported by the **PAMELAMesh**.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_FIELD_SPECS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_FIELD_SPECS_END -->


.. _Outputs_tag_dead_oil_egg_model:

-------
Outputs
-------

In this section, we request an output of the results in VTK format and an output of the rates for each producing well.
Note that the name defined here must match the name used in the **Events** XML block to define the output frequency.

.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_OUTPUTS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_OUTPUTS_END -->

.. _Tasks_tag_dead_oil_egg_model:

-------
Tasks
-------

In the **Events** block, we have defined four events requesting that a task periodically collects the rate for each producing well.
This task is defined here, in the **PackCollection** XML sub-block of the **Tasks** block.
The task contains the path to the object on which the field to collect is registered (here, a ``WellElementSubRegion``) and the name of the field (here, ``wellElementMixtureConnectionRate``).
The details of the history collection mechanism can be found in :ref:`TasksManager`. 


.. literalinclude:: ../../../../../inputFiles/compositionalMultiphaseWell/benchmarks/Egg/deadOilEgg_base_iterative.xml
  :language: xml
  :start-after: <!-- SPHINX_TUT_DEAD_OIL_EGG_TASKS -->
  :end-before: <!-- SPHINX_TUT_DEAD_OIL_EGG_TASKS_END -->
	       

All elements are now in place to run GEOSX.

------------------------------------
Running GEOSX
------------------------------------

The first few lines appearing to the console are indicating that the XML elements are read and registered correctly:

.. code-block:: console

  Adding Solver of type CompositionalMultiphaseReservoir, named coupledFlowAndWells
  Adding Solver of type CompositionalMultiphaseFVM, named compositionalMultiphaseFlow
  Adding Solver of type CompositionalMultiphaseWell, named compositionalMultiphaseWell
  Adding Mesh: PAMELAMesh, mesh
  Adding Mesh: InternalWell, wellProducer1
  Adding Mesh: InternalWell, wellProducer2
  Adding Mesh: InternalWell, wellProducer3
  Adding Mesh: InternalWell, wellProducer4
  Adding Mesh: InternalWell, wellInjector1
  Adding Mesh: InternalWell, wellInjector2
  Adding Mesh: InternalWell, wellInjector3
  Adding Mesh: InternalWell, wellInjector4
  Adding Mesh: InternalWell, wellInjector5
  Adding Mesh: InternalWell, wellInjector6
  Adding Mesh: InternalWell, wellInjector7
  Adding Mesh: InternalWell, wellInjector8
  Adding Event: PeriodicEvent, solverApplications
  Adding Event: PeriodicEvent, vtk
  Adding Output: VTK, vtkOutput
  Adding Object CellElementRegion named reservoir from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion1 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion2 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion3 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion4 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion5 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion6 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion7 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion8 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion9 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion10 from ObjectManager::Catalog. 
  Adding Object WellElementRegion named wellRegion11 from ObjectManager::Catalog.
  Adding Object WellElementRegion named wellRegion12 from ObjectManager::Catalog.
                
This is followed by the creation of the 18553 hexahedral cells of the imported mesh:  

.. code-block:: console

  0 >>> **********************************************************************
  0 >>>                          PAMELA Library Import tool                   
  0 >>> **********************************************************************
  0 >>> GMSH FORMAT IDENTIFIED
  0 >>> *** Importing Gmsh mesh format...
  0 >>> Reading nodes...
  0 >>> Done0
  0 >>> Reading elements...
  0 >>> Reading element data...
  0 >>> Number of nodes = 22227
  0 >>> Number of triangles = 0
  0 >>> Number of quadrilaterals = 0
  0 >>> Number of tetrahedra = 0
  0 >>> Number of hexahedra = 18553
  0 >>> Number of pyramids = 0
  0 >>> Number of prisms = 0
  0 >>> *** Done
  0 >>> *** Creating Polygons from Polyhedra...
  0 >>> 59205 polygons have been created
  0 >>> *** Done
  0 >>> *** Perform partitioning...
  0 >>> TRIVIAL partitioning...
  0 >>> Ghost elements...
  0 >>> Clean mesh...
  0 >>> *** Done...
  0 >>> Clean Adjacency...
  0 >>> *** Done...
                
At this point, we are done with the case set-up and
the code steps into the execution of the simulation itself:

.. code-block:: console

  Time: 0s, dt:10000s, Cycle: 0

    Attempt:  0, NewtonIter:  0
    ( Rfluid ) = (9.39e+01) ;     ( R ) = ( 1.06e+02 ) ; 
    Attempt:  0, NewtonIter:  1
    ( Rfluid ) = (2.14e+00) ;     ( R ) = ( 2.20e+00 ) ; 
    Last LinSolve(iter,res) = (   2, 3.62e-03 ) ; 
    Attempt:  0, NewtonIter:  2
    ( Rfluid ) = (3.23e-01) ;     ( R ) = ( 3.37e-01 ) ; 
    Last LinSolve(iter,res) = (   4, 1.82e-03 ) ; 
    Attempt:  0, NewtonIter:  3
    ( Rfluid ) = (1.07e-02) ;     ( R ) = ( 1.16e-02 ) ; 
    Last LinSolve(iter,res) = (   2, 6.13e-03 ) ; 
    Attempt:  0, NewtonIter:  4
    ( Rfluid ) = (7.46e-05) ;     ( R ) = ( 7.50e-05 ) ; 
    Last LinSolve(iter,res) = (   3, 5.09e-03 ) ;
    
  coupledFlowAndWells: Newton solver converged in less than 15 iterations, time-step required will be doubled.

------------------------------------
Visualization
------------------------------------

A file compatible with Paraview is produced in this example.
It is found in the output folder, and usually has the extension `.pvd`.
More details about this file format can be found 
`here <https://www.paraview.org/Wiki/ParaView/Data_formats#PVD_File_Format>`_.
We can load this file into Paraview directly and visualize results:

|pic1| |pic2|

.. |pic1| image:: pressure.gif
   :width: 45%

.. |pic2| image:: saturation.gif
   :width: 45%

We have instructed GEOSX to output the time series of rates for each producer.
The data contained in the corresponding hdf5 files can be extracted and plotted
as shown below.

.. plot::

   import matplotlib
   import matplotlib.pyplot as plt
   import numpy as np
   import h5py

   def main():

       numWells = 4
       iplt = -1
       cmap = plt.get_cmap("tab10")
    
       # Loop over the four producers 
       for iw in range(1,numWells+1):

           # File path
           hdf5FilePath = 'wellRateHistory'+str(iw)+'.hdf5'    

           # Read HDF5
           hf = h5py.File(hdf5FilePath, 'r')
           time = hf.get('Time')
           time = np.array(time)
           massRate = hf.get('wellElementMixtureConnectionRate')
           massRate = np.array(massRate)

           # Some comments about the computation of the volumetric rate here:
           # A proper oil rate constraint for individual wells is currently being implemented
           # In the meantime, the volume rate is (wrongly) computed by dividing
           # the total mass rate by the surface oil density 
        
           # Conversions
           inCubicMeters = 1/848.9
           inDays = 1.0 / (24 * 3600)
        
           # Plot HDF5 content (here, the rate at the well)
           iplt += 1
           plt.plot(time[:,0]*inDays, abs(massRate[:,0])*inCubicMeters/inDays, '-o', color=cmap(iplt), label='Producer #'+str(iw))

       plt.xlim( -1, 175)
       plt.ylim( 0, 3800)        
       plt.grid()
       plt.xlabel('time [days]')
       plt.ylabel('total rate [cubic meters per day]')
       plt.legend(bbox_to_anchor=(0.025, 0.975), loc='upper left', borderaxespad=0.)
       plt.show()

   if __name__ == "__main__":
       main()

	   
------------------------------------
To go further
------------------------------------

**Feedback on this example**

This concludes the example on setting up a Dead-Oil simulation in the Egg model.
For any feedback on this example, please submit
a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.

**For more details**

  - A complete description of the reservoir flow solver is found here: :ref:`CompositionalMultiphaseFlow`.
  - The well solver is description at :ref:`CompositionalMultiphaseWell`. 
  - The available constitutive models are listed at :ref:`Constitutive`.

