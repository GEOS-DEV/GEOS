.. _particleFileWriter:

############################################
particleFileWriter ("pfw"; MPM Only)
############################################


Overview
========

The ``particleFileWriter.py`` script automates the setup of GEOS-MPM simulations by reading a user-defined input file, constructing geometry objects,
generating particle files, creating the GEOS XML input, and optionally writing and submitting machine-specific batch scripts.
It provides a flexible workflow for defining materials, particle fields, solver parameters, and restart behavior,
while supporting both fresh runs and continuation runs from prior restart states.
Together, these components allow users to efficiently prepare, launch, and manage GEOS-MPM jobs across a variety of simulation configurations.
.. _script_setup_and_initialization:

- :ref:`Script Setup and Initialization <script_setup_and_initialization>`

- :ref:`Reading User Input and Dependencies <reading_user_input_and_dependencies>`

- :ref:`Parameter Definitions and Processing <parameter_definitions_and_processing>`

  - :ref:`Parameter Dictionary <parameter_dictionary>`
  - :ref:`Processing pfw Inputs <process_pfw_inputs>`

- :ref:`Geometry and Domain Setup <geometry_and_domain_setup>`

  - :ref:`Object Construction <object_construction>`
  - :ref:`Material and Subregion Setup <material_and_subregion_setup>`

- :ref:`Validation <validation>`

  - :ref:`Error Checking <error_checking>`

- :ref:`Particle Generation <particle_generation>`

  - :ref:`Create Particle File <create_particle_file>`

- :ref:`Run Configuration <run_configuration>`

  - :ref:`Batch Run Configuration <batch_run_configuration>`
  - :ref:`Continuation of Previous Job <continuation_of_previous_job>`
  - :ref:`Create Input File <create_input_file>`
  - :ref:`Create Batch Script <create_batch_script>`
  - :ref:`Submit Job <submit_job>`
  - :ref:`Automatic Restart Checking <automatic_restart_checking>`

.. _script_setup_and_initialization:

Script Setup and Initialization
-------------------------------

This section of ``particleFileWriter.py`` sets up the execution environment, imports required libraries, parses user input, and configures machine-specific behavior for running GEOS-MPM simulations.

The script begins by importing standard scientific Python libraries (e.g., ``numpy``, ``matplotlib``), utilities for system interaction (e.g., ``os``, ``subprocess``, ``platform``), and custom geometry tools from ``pfw_geometryObjects``. These dependencies support geometry generation, particle placement, visualization, and job submission workflows.

Command-line input is handled using ``argparse``, requiring the user to provide a Particle File Writer (PFW) input file. This file defines the simulation setup and is the primary driver for generating GEOS-MPM input and particle data.

Machine-specific configuration is then performed by detecting the current system hostname and matching it against a predefined list of supported machines (e.g., ``dane``, ``tuolumne``). Based on this, the script sets:

- the number of cores per node
- scheduler type (Flux vs. SLURM)
- MPI execution parameters (rank and total processes)

This allows the script to adapt automatically to different HPC environments.

Finally, a helper function ``compute_similarity`` is defined, which measures how similar two strings are using character-level differences. This is typically used for matching or validating user-provided names (e.g., geometry or material identifiers).

.. raw:: html

   <h4><u>Notes</u></h4>

- The script assumes execution in an HPC environment and may behave differently depending on the detected machine.
- MPI is initialized only when not using the Flux scheduler.
- Warnings related to NumPy string comparisons are suppressed to avoid excessive console output.


:ref:`Back to top <particleFileWriter>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">

.. _reading_user_input_and_dependencies:

Reading User Input and Dependencies
-----------------------------------

This section initializes the Particle File Writer workflow by locating user-specific settings, reading the requested PFW input file, and preparing any additional dependencies needed for the run.
It first determines the current username and attempts to import a corresponding ``userDefs_<username>.py`` file, which provides machine- and user-specific paths such as ``geosPath`` and ``pfwPath``. These settings allow the script to find the GEOS installation and the local Particle File Writer utilities.

The script then processes the user-supplied input file name and scans that file for any lines tagged with ``#[pfw_dependency]``. These tags identify additional Python files that must be copied into the current working directory before the input file is imported. Each dependency is copied safely using a temporary filename and rename step, followed by a short visibility check to guard against stale filesystem metadata on Lustre systems. After invalidating Python’s import cache, the script imports the input file as a module and begins building run-specific metadata such as a timestamp and the current working directory.

.. raw:: html

   <h4><u>Notes</u></h4>

- The ``userDefs_<username>.py`` file must exist and contain valid definitions for ``geosPath`` and ``pfwPath``.
- Dependency files listed with ``#[pfw_dependency]`` are copied automatically before the input file is imported.
- The temporary copy-and-rename pattern helps ensure that copied files are complete and visible before Python tries to import them.
- This section prepares the environment for the rest of the script by making sure all required input definitions and helper modules are available.

:ref:`Back to top <particleFileWriter>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _parameter_definitions_and_processing:

Parameter Definitions and Processing
------------------------------------

This section defines the configurable parameters used by the Particle File Writer and explains how user-provided values in the ``pfw`` dictionary are interpreted, checked, and converted into GEOS solver settings.

=====

.. _parameter_dictionary:


Parameter Dictionary
^^^^^^^^^^^^^^^^^^^^

This section defines the full set of configurable parameters used by the Particle File Writer (PFW).
Each entry in the ``parameters`` dictionary contains:

- a default value
- a flag indicating whether the parameter should be written to the GEOS MPM XML solver input

If a parameter is not explicitly provided in the user input file, its default value is used.



Core Execution and Job Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``runDebug``
     - Run in debug mode (short test runs)
     - False
   * - ``generateParticleFile``
     - Generate particle file before running
     - True
   * - ``runContinuation``
     - Continue from previous simulation
     - False
   * - ``restartJobDir``
     - Directory for restart files
     - "."
   * - ``restartCycleNum``
     - Restart cycle index
     - 0

HPC / Job Submission Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``mNodes``, ``mCores``
     - Number of nodes and cores
     - 1
   * - ``mPartition``
     - Scheduler partition
     - pbatch
   * - ``mWallTime``
     - Job wall time
     - 00:30:00
   * - ``mSubmitJobs``
     - Automatically submit jobs
     - False
   * - ``autoRestart``
     - Automatically restart failed jobs
     - False

Domain and Discretization
~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``xmin, xmax, ymin, ymax, zmin, zmax``
     - Domain bounds
     - must be user defined
   * - ``nI, nJ, nK``
     - Grid resolution
     - 5
   * - ``ppc``
     - Particles per cell
     - 2
   * - ``periodic``
     - Periodic boundary conditions
     - [False, False, False]

Solver and Time Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``timeIntegrationOption``
     - Time integration scheme
     - None
   * - ``updateMethod``, ``updateOrder``
     - MPM update scheme
     - None
   * - ``cflFactor``
     - CFL stability factor
     - None
   * - ``initialDt``
     - Initial timestep
     - None
   * - ``endTime``
     - Total simulation time
     - 1.0

Output and History
~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``reactionHistory``
     - Enable reaction output
     - None
   * - ``reactionWriteInterval``
     - Reaction output frequency
     - None
   * - ``boxAverageHistory``
     - Enable box-average output
     - None
   * - ``plotInterval``
     - Plot output interval
     - None
   * - ``restartInterval``
     - Restart output interval
     - None

Physics and Material Behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``materials``
     - List of material models
     - None
   * - ``materialPropertyString``
     - Material definition block
     - None
   * - ``bodyForce``
     - External body force
     - None
   * - ``frictionCoefficient``
     - Contact friction
     - None
   * - ``damageFieldPartitioning``
     - Damage field control
     - None

Contact and Cohesive Behavior
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :widths: 30 40 30
   :header-rows: 1

   * - Parameter
     - Description
     - Default
   * - ``enableCohesiveFailure``
     - Enable cohesive fracture
     - 0
   * - ``maxCohesiveNormalStress``
     - Cohesive strength (normal)
     - 0.01
   * - ``maxCohesiveShearStress``
     - Cohesive strength (shear)
     - 0.01
   * - ``characteristicNormalDisplacement``
     - Cohesive length scale
     - 0.01

=====

Particle File Fields
~~~~~~~~~~~~~~~~~~~~

The following defines which fields are written to the particle file and their ordering:

.. code-block:: python

   particleFieldOrder = [
       "Velocity",
       "MaterialType",
       "ContactGroup",
       "SurfaceFlag",
       "Damage",
       "Porosity",
       "Temperature",
       "StrengthScale",
       "CZTag", # cohesive zone = CZ
       "RVector",
       "MaterialDirection",
       "SurfaceNormal",
       "SurfacePosition",
       "SurfaceTraction",
       "ShrinkageFlag"
   ]

.. raw:: html

   <h4><u>Notes</u></h4>

- The order of ``particleFieldOrder`` must match the order used when writing particle data.
- The ``parameters`` dictionary acts as a centralized configuration system for both simulation setup and XML generation.
- Parameters marked for XML inclusion are automatically written to the GEOS solver input file when specified.

=====

.. _process_pfw_inputs:

Processing ``pfw`` Inputs
^^^^^^^^^^^^^^^^^^^^^^^^^

This section reads the ``pfw`` dictionary from the user input file, applies default values where needed, and prepares the solver parameter string that will later be written into the GEOS MPM XML input. It also promotes recognized parameters to global variables so they can be accessed throughout the Particle File Writer script.

For each parameter listed in the predefined ``parameters`` dictionary, the script checks whether the user supplied a value in ``job.pfw``.
If so, that value is used; otherwise, the default value is retained. Parameters marked for XML inclusion are converted into XML-style attribute strings and collected in ``parameterStrings``. Any dependency files identified earlier are added temporarily to the ``pfw`` dictionary for bookkeeping, then removed before XML generation so they are not passed through as solver parameters.

The script also checks for user-provided parameters that are not in the known parameter list. If an unknown parameter is detected, it searches for a similar known name using ``compute_similarity`` and prints a suggestion when possible. These unknown parameters are still added to the XML parameter string, which helps preserve flexibility while warning the user about possible typos.

Finally, the script verifies that required particle fields are present when explicit surface contact or cohesive-zone (CZ) options are enabled.
If needed, ``SurfaceNormal`` and ``SurfacePosition`` are added automatically to ``particleFileFields``. The section then finalizes the XML solver parameter string, trims ``particleFieldOrder`` to include only fields requested by the user, and initializes the interior-domain discretization variables used later in the script.

.. raw:: html

   <h4><u>Notes</u></h4>

- The ``pfw`` dictionary is the main user-facing parameter interface for the script.
- Recognized parameters are assigned as global variables so later sections of the script can use them directly.
- Unknown parameters are not discarded; they are passed through to the XML with a warning.
- Contact and cohesive-zone (CZ) settings may require additional particle fields, which are added automatically if missing.
- ``ppcx``, ``ppcy``, and ``ppcz`` default to ``ppc`` when directional particle-per-cell values are not explicitly provided.

:ref:`Back to top <particleFileWriter>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">

.. _geometry_and_domain_setup:


Geometry and Domain Setup
-------------------------

This section defines how simulation geometry is constructed and how material and subregion information is prepared prior to particle generation.

=====

.. _object_construction:


Object Construction
^^^^^^^^^^^^^^^^^^^

The script determines how geometric objects are generated for the simulation. Users can either directly provide a list of objects through the ``pfw`` dictionary or define a ``make_objects`` function in the input file.

If a ``make_objects`` function is present, the script inspects its signature:

- If no arguments are required, all MPI ranks construct the full set of objects.
- If arguments are present, each rank constructs only the objects within its assigned spatial slice (based on the x-domain).

This rank-aware construction improves performance for simulations with many objects (e.g., granular systems or Voronoi tessellations) by reducing redundant geometry generation.

If neither ``objects`` nor ``make_objects`` is defined, the script issues an error.

.. raw:: html

   <h4><u>Notes</u></h4>

- Rank-based object construction is recommended for large simulations.
- Users can explicitly pass ``objects=[]`` if no geometry is required.
- Geometry objects must support the required interface (e.g., ``getSubregions``).

:ref:`Back to top <particleFileWriter>`

=====

.. _material_and_subregion_setup:

Material and Subregion Setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section initializes material definitions and associates geometry subregions with material and particle-type identifiers.

The script first checks that a ``materials`` array is provided in the input file. If not, execution is aborted.
The materials list is then formatted into a string representation suitable for inclusion in the GEOS XML input.

Next, the script gathers all subregions from the defined geometry objects. Each subregion corresponds to a specific material index and particle type. These are collected and grouped by material to determine how many unique particle types exist per material.

The total number of subregions is computed and reported, providing a measure of geometric and material complexity.

Finally, the script initializes the ``particleRefinement`` array, which controls per-material particle resolution.
If not specified by the user, it defaults to a value of 1 for all materials.

.. raw:: html

   <h4><u>Notes</u></h4>

- The ``materials`` array is required and defines the material models used in the simulation.
- Each geometry object contributes subregions that map to material indices and particle types.
- The number of subregions affects particle generation and material assignment.
- ``particleRefinement`` can be used to locally refine particle resolution for specific materials.


:ref:`Back to top <particleFileWriter>`


.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">

  


.. _validation:

Validation
----------

This section performs basic validation checks on the simulation configuration before proceeding with particle generation and job execution.

=====

.. _error_checking:


Error Checking
^^^^^^^^^^^^^^

This section performs basic validation checks on the simulation configuration before proceeding with particle generation and job execution. These checks are only executed on the root MPI rank (``rank == 0``) to avoid redundant output across processes.

Currently, the script verifies that plane strain simulations are configured correctly. Specifically:

- If ``planeStrain == 1``, the number of grid cells in the z-direction (``nK``) must be equal to 3.
- If ``planeStrain == 1``, the domain partitioning in the z-direction (``zpar``) must be equal to 1.

If either of these conditions is violated, an error message is printed and execution is halted. The script terminates using:

- ``comm.Abort()`` when running under MPI
- ``sys.exit()`` when running under the Flux scheduler

.. raw:: html

   <h4><u>Notes</u></h4>

- These checks help ensure that plane strain simulations are properly constrained in the out-of-plane direction.
- The requirement ``nK = 3`` reflects the use of a minimal three-cell thickness for plane strain MPM formulations.
- The requirement ``zpar = 1`` ensures that no parallel decomposition is applied in the constrained direction.
- Additional validation checks could be added here to enforce consistency in other simulation configurations.



:ref:`Back to top <particleFileWriter>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">


.. _particle_generation:


Particle Generation
-------------------

This section describes the core workflow used to generate the GEOS-MPM particle file.

=====

.. _create_particle_file:


Create Particle File
^^^^^^^^^^^^^^^^^^^^

This section generates the GEOS-MPM particle file by discretizing the simulation domain, testing each candidate particle location against the user-defined geometry objects, 
assigning particle properties, and writing the resulting particle data to file. For parallel runs, each MPI rank generates particles for its assigned portion of the domain, 
and the per-rank files are merged into a single particle file at the end.

=====

Domain and Particle Discretization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The script first defines the particle file name and computes the interior grid resolution, taking into account periodic boundaries and plane strain settings. It then calculates:

- interior cell counts in each direction
- grid cell spacing
- total domain bounds including ghost cells
- number of particles in each direction
- particle dimensions

These quantities determine the candidate particle locations that will be tested against the geometry.

For plane strain simulations, only one particle layer is used in the z-direction. Otherwise, the particle count in each direction is based on the particles-per-cell settings.

.. raw:: html

   <h4><u>Notes</u></h4>

- ``nI``, ``nJ``, and ``nK`` represent the interior discretization.
- ``NI``, ``NJ``, and ``NK`` include ghost cells where needed for the GEOS XML input.
- In plane strain, the particle thickness is taken as the full out-of-plane cell thickness.

=====

Particle File Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``generateParticleFile`` is enabled, the script opens a rank-local particle file and prepares for particle generation. It also computes a characteristic ``surfaceDepth`` used to determine whether a particle lies near an object boundary. This value depends on the grid spacing and the particles-per-cell refinement.

The script then computes representative partition-center indices for each processor partition in x, y, and z. These are printed for diagnostic purposes and can help users verify that the domain decomposition is behaving as expected.

.. raw:: html

   <h4><u>Notes</u></h4>

- Each MPI rank writes to a temporary file named ``mpmParticleFile_<job>_<rank>``.
- ``surfaceDepth`` is used by geometry objects to classify particles as interior, surface, cohesive, or fully damaged.

=====

Particle Generation Loop
~~~~~~~~~~~~~~~~~~~~~~~~

The main particle-generation loop iterates through candidate particle locations across the assigned portion of the domain. For each particle center:

#. The physical particle coordinates are computed.
#. The script determines which geometry objects should be checked.
#. The particle is tested against each object using ``object.isInterior(pt, surfaceDepth)``.
#. If the particle lies inside an object, the script retrieves the particle's material and state information.
#. If requested, the particle may be further refined into smaller sub-particles according to ``particleRefinement``.

When ``sortObjects`` is enabled, the script first reduces the object list to those whose x-bounds overlap the current x-slice. This can significantly improve performance for simulations with many objects.

.. raw:: html

   <h4><u>Notes</u></h4>

- Only the first matching geometry object is used for a given particle location.
- Setting ``mat = -1`` provides a way to define voids or defects, since particles with negative material ID are not written.
- Progress information is printed during particle creation, including an estimated remaining time.

=====

Assigned Particle Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For each accepted particle, the script assigns a set of fields based on the object definition and the contents of ``particleFileFields``. These may include:

- particle type
- velocity
- material type
- contact group
- surface flag
- damage
- porosity
- temperature
- strength scaling
- cohesive-zone tag ("CZ flag")
- characteristic particle size tensor (``RVector``)
- material direction
- surface normal
- surface position
- surface traction
- shrinkage flag

Whenever possible, the script first checks for object methods such as ``getVelocity`` or ``getDamage``. If these are not available, it falls back to object attributes or callable attribute definitions.

Special options such as ``useSinusoidalDamageField`` and ``wavyCrack`` can overwrite the default damage assignment to impose prescribed spatial damage patterns.

.. raw:: html

   <h4><u>Notes</u></h4>

- Surface-related quantities are only evaluated when the particle is marked as non-interior.
- Plane strain enforces zero out-of-plane velocity and surface-position components where appropriate.
- Material directions may be written either as a 3-component vector or as a full 3x3 tensor, depending on the object definition.

=====

Validation of Particle Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before writing a particle, the script checks that all requested particle properties have been assigned valid values. If any required value is ``None``, the script prints an error message identifying the object and field, then exits immediately.

This provides a safeguard against incomplete object definitions and helps catch mistakes in custom geometry classes.

.. raw:: html

   <h4><u>Notes</u></h4>

- Required checks are performed only for fields requested in ``particleFileFields``.
- These checks are especially helpful when users define custom geometry objects with optional methods.

======

Particle Refinement and Writing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After determining the base particle properties, the script applies any material-specific refinement requested through ``particleRefinement``. Each parent particle may be subdivided into a regular grid of refined sub-particles. The refined particle coordinates and dimensions are computed, and each valid particle is written to the rank-local particle file.

The particle record always includes:

- particle ID
- particle position
- particle type

Additional fields are appended in the order defined by ``particleFieldOrder``.

The script also accumulates the total particle volume so that a final particle volume fraction can be reported.

.. raw:: html

   <h4><u>Notes</u></h4>

- Refinement is applied separately for each material, allowing selective local resolution increases.
- The written field order must remain consistent with the expected particle file header.
- ``RVector`` stores the particle half-dimensions in tensor form for GEOS-MPM.

=====

Parallel File Merge
~~~~~~~~~~~~~~~~~~~

Once all MPI ranks finish particle generation, the root rank merges the rank-local particle files into a single final particle file. During this step:

- a complete header line is written
- rank-local particles are appended in order
- particle IDs are renumbered globally
- temporary per-rank files are deleted

The script then reports:

- total number of particles
- total particle volume
- physical domain volume
- particle volume fraction
- average number of particles per rank

.. raw:: html

   <h4><u>Notes</u></h4>

- The final output file is named ``mpmParticleFile_<job>``.
- The merge step is performed only on rank 0.
- Particle IDs in the final file are sequential and global, regardless of the original rank-local numbering.

=====

Skipping Particle Generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If ``generateParticleFile`` is disabled, the script skips this entire section and prints a message indicating that particle file generation was not performed.

.. raw:: html

   <h4><u>Notes</u></h4>

- Particle generation is also automatically disabled for continuation runs when ``runContinuation`` is active.


:ref:`Back to top <particleFileWriter>`

.. raw:: html

   <hr style="height:3px; border:none; background-color:black;">

.. _run_configuration:


Run Configuration
-----------------

This section defines how simulations are configured for execution, including resource calculations, restart behavior, XML creation, batch script generation, job submission, and optional automated restart checks.


=====

.. _batch_run_configuration:

Batch Run Configuration
^^^^^^^^^^^^^^^^^^^^^^^

The script computes machine-dependent batch parameters for running GEOS-MPM jobs.

- The partition is set to ``pdebug`` when ``runDebug=True``, otherwise the user-defined partition is used.
- The total number of MPI cores is computed as:

  .. math::

     mCores = xpar \cdot ypar \cdot zpar

- The number of nodes is automatically determined based on the machine-specific cores per node.
- If no bank is specified, a default value is loaded from the user definitions file.
- The wall time is parsed and, for debug runs, capped at one hour.
- A restart-safe buffer time is computed to ensure sufficient time for writing restart files.
- The script estimates the simulation cost in core-hours and prints an approximate charge to the user’s allocation.

.. raw:: html

   <h4><u>Notes</u></h4>

- Users do not need to manually set ``mCores``; it is derived from domain partitioning.
- Node count is automatically adjusted for the target machine.
- Debug runs are constrained to short wall times for queue efficiency.
- Reported cost is an estimate and may vary depending on system accounting.

:ref:`Back to top <particleFileWriter>`

=====

.. _continuation_of_previous_job:


Continuation of Previous Job
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section prepares a simulation to restart from a previous GEOS-MPM run. If ``runContinuation`` is enabled, the script reads the restart job directory and restart cycle number from the user input file, verifies that the restart directory exists, and constructs the expected restart file name.

The corresponding restart directory and ``.root`` file are then copied into the current working directory so the new run can start from the saved state. The particle file from the previous job is also copied and renamed to match the current job, even though GEOS will read from the restart state rather than regenerate the particle distribution.

.. raw:: html

   <h4><u>Notes</u></h4>

- Continuation runs require both ``restartJobDir`` and ``restartCycleNum`` to be specified in the ``pfw`` dictionary.
- The restart file is named using the pattern ``mpm_<jobName>_restart_<cycle>``.
- Particle generation is typically disabled for continuation runs, since the restart file contains the necessary particle state.
- The string ``restartStr`` is later appended to the GEOS launch command so the solver starts from the restart state.

:ref:`Back to top <particleFileWriter>`

.. _create_input_file:

Create Input File
^^^^^^^^^^^^^^^^^

This section writes the GEOS XML input file that defines the background mesh, particle mesh, material regions, solver settings, constitutive models, events, numerical methods, and outputs for the simulation. Only rank 0 performs this step.

The script first determines which particle blocks and particle types belong to each material by examining the unique particle types present in each material's subregions. It then constructs:

- a particle block string for the particle mesh
- a particle type string describing the interpolation method for each block
- a particle region string assigning each block to a material
- a target-region string for the MPM solver

These strings are inserted into the XML template together with the domain bounds, grid resolution, periodicity, solver parameter string, material definitions, event schedule, and output settings. The final XML file is written as:

.. code-block:: text

   mpm_<jobname>.xml

.. raw:: html

   <h4><u>Notes</u></h4>

- Particle type IDs are mapped to GEOS particle interpolation names such as ``SinglePoint``, ``CPDI``, and ``CPTI``.
- Each material may contain multiple particle blocks if multiple subregions or particle types are present.
- The solver parameter string ``mpmSolverParameterString`` is inserted directly into the ``SolidMechanics_MPM`` block.
- If ``mpmEventsString`` is provided, it is included as a nested ``<MPMEvents>`` block in the solver definition.
- Output events include both VTK output and restart output, along with a halt event based on the available wall time.

:ref:`Back to top <particleFileWriter>`

.. _create_batch_script:

=====

Create Batch Script
^^^^^^^^^^^^^^^^^^^

This section creates a machine-specific batch script for launching the GEOS simulation. The script supports both Flux- and SLURM-based systems and writes the appropriate job script depending on the detected scheduler.

For Flux systems, the script writes a Flux batch file with the requested node count, wall time, bank, and queue, then launches GEOS using ``flux run``. For SLURM systems, it writes a standard ``sbatch`` script with the requested wall time, node count, task count, account, and partition, then launches GEOS using ``srun``.

In both cases, the batch script sets ``OMP_NUM_THREADS=1`` and passes the GEOS executable, XML input file, domain partitioning, and optional restart string to the run command.

.. raw:: html

   <h4><u>Notes</u></h4>

- The generated run script is named ``<timestamp>_runGEOS.sh``.
- On SLURM systems, ``#SBATCH --export=NONE`` is included to avoid issues caused by the modified MPI environment from ``mpi4py``.
- The GEOS executable is launched with the partition counts ``-x``, ``-y``, and ``-z`` so the solver uses the requested spatial decomposition.
- Restart runs include the ``-r`` flag pointing to the copied restart file.


:ref:`Back to top <particleFileWriter>`

=====

.. _submit_job:

Submit Job
^^^^^^^^^^

If ``mSubmitJobs`` is enabled, this section submits the generated GEOS batch script to the scheduler. The script uses:

- ``flux batch`` on Flux systems
- ``sbatch`` on SLURM systems

After submission, the scheduler output is parsed to extract the job ID, which is then printed for the user.

.. raw:: html

   <h4><u>Notes</u></h4>

- Job submission is optional and controlled by the ``mSubmitJobs`` flag.
- The parsed job ID is later reused if automatic restart checking is enabled.
- Submission behavior depends on the detected scheduler.

=====

.. _automatic_restart_checking:


Automatic Restart Checking
^^^^^^^^^^^^^^^^^^^^^^^^^^

If both ``mSubmitJobs`` and ``autoRestart`` are enabled, this section creates and submits a secondary batch script that runs after the GEOS job finishes. This follow-up job launches ``pfw_check.py`` to inspect the completed run and determine whether the simulation should be restarted.

The check job is submitted with a scheduler dependency so that it runs only after the GEOS job exits. On SLURM systems, this is done using:

.. code-block:: text

   --dependency=afterany:<jobID>

The generated script runs ``python3 pfw_check.py <inputFile> <jobID>`` and can be used to automate restart handling for long simulations.

.. raw:: html

   <h4><u>Notes</u></h4>

- Automatic restart checking is currently supported only for SLURM systems.
- If Flux is detected, the script prints a warning that ``autoRestart`` is not supported.
- The generated restart-check script is named ``<timestamp>_runCheck.sh``.
- This workflow is useful for simulations that may need multiple sequential submissions to reach the desired end time.

:ref:`Back to top <particleFileWriter>`