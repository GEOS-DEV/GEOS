.. _triaxialDriver:


###########################################################
Triaxial Driver
###########################################################


**Context**

The main goal of this example is to learn how to use the ``TriaxialDriver`` to examine solid mechanics models by simulating rock mechanics experiments, such as triaxial tests. 
Necessarily, prior to field scale simulations with testing model, predictions (e.g., stress-strain relationships) obtained from running ``TriaxialDriver`` should be compared with experimental measurements to calibrate material properties.


**Objectives**

At the end of this example, you will know:

 - how to set up a triaxial test scenario with the ``TriaxialDriver``,
 - how to define a material model,
 - how to specify loading conditions with ``TableFunction``,
 - how to save and postprocess simulation results.


**Input file**

The XML file for this test case is located at:

.. code-block:: console

  inputFiles/solidMechanics/triaxialDriver_ExtendedDruckerPrager.xml


This example also uses a set of table files located at:

.. code-block:: console

  inputFiles/solidMechanics/triaixalDriverTables/


A python script for post-processing the simulation results is also provided:

.. code-block:: console

  src/docs/sphinx/basicExamples/triaxialDriver/triaxialDriverFigure.py


------------------------------------------------------------------
Description of the case
------------------------------------------------------------------

Instead of launching a full finite element simulation to mimic experimental loading conditions, GEOSX provides a ``TriaxialDriver`` module to investigate constitutive behaviors and simulate laboratory tests. The triaxial driver makes it easy to interpret the mechanical response and calibrate the constitutive models with available experimental data. 

In this example, the Extended Drucker-Prager model (see :ref:`DruckerPragerExtended`) is applied to solve for elastoplastic deformation of rock samples, when subjected to controlled loading conditions. The strain hardening Extended Drucker-Prager model with an associated plastic flow rule in GEOSX is tested in this example. To replicate a typical triaxial test, table function is hereby employed to specify loading conditions in both axial and radial direction. The resulted strains and stresses in both directions are numerically calculated and can be saved into an output file in a simple ASCII format.


For this example, we focus on the ``Task``, the ``Constitutive``, and the ``Functions`` tags.

------------------------------------------------------------------
Task
------------------------------------------------------------------

In GEOSX, the ``TriaxialDriver`` is defined with a particular XML structure. 
The ``TriaxialDriver`` is added to a standard XML input deck as a solo task to the ``Tasks`` queue and added as a ``SoloEvent`` to the event queue.

For this example, we simulate the elastoplastic deformation of a confined specimen caused by external load. A homogeneous domain with one solid material is assumed. The material is named as ``ExtendedDruckerPrager``, whose mechanical properties are specified in the ``Constitutive`` section.

Different testing modes are available in the ``TriaxialDriver`` to mimic different loading conditions in various core sample tests:
 
+--------------------+-------------------------+--------------------------+---------------------------+
| **mode**           | **axial loading**       | **radial loading**       | **initial stress**        |
+--------------------+-------------------------+--------------------------+---------------------------+
| ``strainControl``  | axial strain controlled | radial strain controlled | isotropic stress using    |
|                    | with ``axialControl``   | with ``radialControl``   | ``initialStress``         |
+--------------------+-------------------------+--------------------------+---------------------------+
| ``stressControl``  | axial stress controlled | radial stress controlled | isotropic stress using    |
|                    | with ``axialControl``   | with ``radialControl``   | ``initialStress``         |
+--------------------+-------------------------+--------------------------+---------------------------+
| ``mixedControl``   | axial strain controlled | radial stress controlled | isotropic stress using    |
|                    | with ``axialControl``   | with ``radialControl``   | ``initialStress``         |
+--------------------+-------------------------+--------------------------+---------------------------+

A triaxial test is usually conducted on a confined rock sample to determine material properties. 
As shown, a conventional triaxial test is described using the ``mode="mixedControl"`` testing mode in the ``TriaxialDriver``.

In a triaxial test, the testing sample is under confining pressure (radial stresses) and subjected to increased axial load. Therefore, a stress control ``radialControl="stressFunction"`` is defined in the radial direction to impose confining pressure; whereas a strain control ``axialControl="strainFunction"`` is applied in the axial direction to represent axial compression.
 
The initial stress is specified by ``initialStress="-10.e6"``. To ensure the static equilibrium at the first timestep, its value should be consistent with the initial set of applied stresses defined in axial or radial loading functions.
This stress has a negative value due to the negative sign convention for compressive stress in GEOSX.


``steps="200" `` defines the number of load steps and ``output="simulationResults.txt"`` specifies an output file to which the simulation results will be written. 


.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaxialDriver_ExtendedDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_TASK -->
    :end-before: <!-- SPHINX_TASK_END -->


Besides triaxial test, volumetric and Oedometer tests can be handled by using the ``strainControl`` mode, if properly defining loading controls in axial and radial direction. A volumetric test can be modelled by setting the axial and radial control functions to the same strain function, whereas an oedometer test can be running by constraining the radial strain to zero.


------------------------------
Constitutive laws
------------------------------

Any solid material model implemented in GEOSX can be called by the ``TriaxialDriver``.

For this problem, Extended Drucker Prager model ``ExtendedDruckerPrager`` is used to describe the mechanical behavior of an isotropic material, when subjected to external loading.
As for the material parameters, ``defaultInitialFrictionAngle``, ``defaultResidualFrictionAngle`` and ``defaultCohesion`` denote the initial friction angle, the residual friction angle, and cohesion, respectively, as defined by the Mohr-Coulomb failure envelope.
As the residual friction angle ``defaultResidualFrictionAngle`` is larger than the initial one ``defaultInitialFrictionAngle``, a strain hardening model is adopted, whose hardening rate is given as ``defaultHardening="0.0001"``. 
If the residual friction angle is set to be less than the initial one, strain weakening will take place. 
Setting ``defaultDilationRatio="1.0"`` corresponds to an associated flow rule.


.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaxialDriver_ExtendedDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_MATERIAL -->
    :end-before: <!-- SPHINX_MATERIAL_END -->


All constitutive parameters such as density, viscosity, bulk modulus, and shear modulus are specified in the International System of Units.


------------------------------------------------------------------
Functions
------------------------------------------------------------------

The ``TriaxialDriver`` uses a simple form of time-stepping to advance through the loading steps, allowing for simulating both rate-dependent and rate-independent models.  

In this case, users should define two different time history functions (``strainFunction`` and ``stressFunction``) to describe loading conditions in axial and radial direction respectively. More specifically, the table functions in this example include the temporal variations of radial stress and axial strain, which rely upon the external files in the table directory.


.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaxialDriver_ExtendedDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_FUNCTION -->
    :end-before: <!-- SPHINX_FUNCTION_END -->


The ``strainFunction`` TableFunction is an example of a 1D interpolated function, which describes the strain as a function of time ``inputVarNames="{ time }"``. This table is defined using a single coordinateFile:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaixalDriverTables/time.geos
  :language: none

And a single voxelFile:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaixalDriverTables/axialStrain.geos
  :language: none


Similarly, the correlation between the confining stress and time is given through the ``stressFunction``, which is defined using the same coordinateFile:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaixalDriverTables/time.geos
  :language: none

And a different voxelFile:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaixalDriverTables/radialStress.geos
  :language: none


For this specific test, the confining stress is kept as a constant, which equals to the ``initialStress``. 
Instead of monotonic change of axial load, two loading/unloading cycles are defined in the ``strainFunction``. 
In this way, both plastic loading and elastic unloading can be detected and observed. 

Note that ``stressFunction`` and ``strainFunction`` have negative values for a compressive test.


------------------------------------------------------------------
Mesh
------------------------------------------------------------------

Even discretization is not required for the ``TriaxialDriver``, a dummy mesh should be defined to pass all the necessary checks when initializing GEOSX and running the module. A dummy mesh should be created in the ``Mesh`` section and assigned to the ``cellBlocks`` in the ``ElementRegions`` section. 


.. literalinclude:: ../../../../../inputFiles/solidMechanics/triaxialDriver_ExtendedDruckerPrager.xml
    :language: xml
    :start-after: <!-- SPHINX_MESH -->
    :end-before: <!-- SPHINX_MESH_END -->


Once calibrated, the testing constitutive models can be easily extended for full field scale simulation by adding solver, discretization, and boundary condition blocks to the xml file. And this example can also be running with full GEOSX model and generate identical results as provided by the ``TriaxialDriver``.


---------------------------------
Running  TriaxialDriver
---------------------------------

The ``TriaxialDriver`` is launched as any other GEOSX simulation by using the following command:

.. code-block:: sh

   path/to/geosx -i triaxialDriver_ExtendedDruckerPrager.xml


The running log appears to the console to indicate if the case can be successfully executed or not:
 
.. code-block:: sh

   Max threads: 32
   MKL max threads: 16
   GEOSX version 0.2.0 (HEAD, sha1: bb16d72)
   Adding Event: SoloEvent, triaxialDriver
      TableFunction: strainFunction
      TableFunction: stressFunction
   Adding Mesh: InternalMesh, mesh1
   Adding Object CellElementRegion named dummy from ObjectManager::Catalog.
   Total number of nodes:8
   Total number of elems:1
   Rank 0: Total number of nodes:8
	   dummy/cellBlock01 does not have a discretization associated with it.
   Time: 0s, dt:1s, Cycle: 0
   Cleaning up events
   Umpire            HOST sum across ranks:   23.2 KB
   Umpire            HOST         rank max:   23.2 KB
   total time                         0.435s
   initialization time                0.053s
   run time                           0.004s


---------------------------------
Inspecting results
---------------------------------

The simulation results are saved in a text file, which is named as ``simulationResults.txt``.
This output file is format with a brief header to mark the meaning of each column. Each row corresponds to one timestep of the driver, starting from initial conditions in the first row. 


.. code:: sh

  # column 1 = time
  # column 2 = axial_strain
  # column 3 = radial_strain_1
  # column 4 = radial_strain_2
  # column 5 = axial_stress
  # column 6 = radial_stress_1
  # column 7 = radial_stress_2
  # column 8 = newton_iter
  # column 9 = residual_norm
  0.0000e+00 0.0000e+00 0.0000e+00 0.0000e+00 -1.0000e+07 -1.0000e+07 -1.0000e+07 0.0000e+00 0.0000e+00 
  2.5000e-02 -1.0000e-04 2.5000e-05 2.5000e-05 -1.1500e+07 -1.0000e+07 -1.0000e+07 1.0000e+00 0.0000e+00 
  5.0000e-02 -2.0000e-04 5.0000e-05 5.0000e-05 -1.3000e+07 -1.0000e+07 -1.0000e+07 1.0000e+00 0.0000e+00 
  7.5000e-02 -3.0000e-04 7.5000e-05 7.5000e-05 -1.4500e+07 -1.0000e+07 -1.0000e+07 1.0000e+00 0.0000e+00 
  1.0000e-01 -4.0000e-04 1.0000e-04 1.0000e-04 -1.6000e+07 -1.0000e+07 -1.0000e+07 1.0000e+00 0.0000e+00 
  ...

Please note that the file contains two columns for radial strain (``radial_strain_1`` and ``radial_strain_2``) and two columns for radial stress (``radial_stress_1`` and ``radial_stress_2``). For isotropic materials, the stresses and strains along the two radial axes would be the same. However, the stresses and strains in the radial directions could potentially differ for the cases with anisotropic materials and true-triaxial loading conditions.

This output file can be processed and visualized using any tools. As an example, with the provided python script, the simulated stress-strain curve, p-q diagram and relationship between volumetric strain and axial strain can be plotted, which could be then validated with experimental observations:


.. plot:: docs/sphinx/basicExamples/triaxialDriver/triaxialDriverFigure.py


------------------------------------------------------------------
To go further
------------------------------------------------------------------


**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOSX/GEOSX/issues>`_.
