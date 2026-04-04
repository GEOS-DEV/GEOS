.. _Input Files MPM:

.. |below| raw:: html

   <span style="color:#999;">(below)</span>

.. |separate| raw:: html

   <span style="color:#999;">(separate page)</span>

############################################
MPM Input Files
############################################


.. toctree::
   :maxdepth: 1
   :hidden:


   pfw_input
   particleFileWriter
   pfw_geometryObjects
   pfw_analysis


Overview
========================
A few things to note. For new users, the Material Point Methods (MPM) is set up such that you define a grid on which you place particles.
The particles have assocaited material masses, momentum, internal and external forces associated with them. 
Later on, you'll see we define what is assigned to the particles in our 'pfw_input' files. 
After set up, we exectue our timestep:

#. Particle attributes are mapped to the closest node on the grid. 
#. Then, momentum equations associated with the grid are solved and the grid is updated. 
#. Attributes are updated on the grid and then mapped to the partcles. These attributes include particle position, velocity, volume, density, deformation gradent, stress, and internal variable (e.g. temperaure).
#. The grid is resent to the original state to avoid distortion.

In our simulations, we define various components of this process (e.g., length of time step, particle mass, force and/or stress on system)

The input files needed for running GEOS-MPM are the following. Details are in their respective links:

- **Environment Set Up**
  
  - :ref:`runClean_username.sh <runClean>` |below|
  - :ref:`userDefs_username.py <userDefs>` |below|

- **Simulation definition**
  
  - :doc:`pfw_input_simulationName.py <pfw_input>` |separate|
  - :doc:`pfw_geometryObjects.py <pfw_geometryObjects>` |separate|

- **Simulation Execution**
  
  - :doc:`particleFileWriter.py <particleFileWriter>` |separate|
  - :ref:`pfw_check.py <pfw_check>` |below|

- **Post-processing**
  
  - :doc:`pfw_analysis.py <pfw_analysis>` |separate|


**If you are not developing**, you will mostly 1) set up your simulations with the pfw_input...py file, 2) execute it with runClean...sh and need to set up your 3) userDefs....py file. 



.. raw:: html

   <hr style="height:4px; border:none; background-color:black;">


.. _runClean:

Environment Set Up: runClean_username.sh
========================
Overview
--------

This script is a SLURM-based launcher and preprocessing driver for GEOS-MPM runs. Its purpose is
to take one or more ``pfw_input_*.py`` files, create corresponding run directories, copy the
required support scripts, and execute ``particleFileWriter.py`` so that GEOS input files and
batch scripts are generated and optionally submitted.

The ``runClean`` script serves as a preprocessing and launch utility that:

- Automates setup of GEOS-MPM simulations
- Ensures consistent directory structure
- Copies required dependencies
- Executes preprocessing for one or more cases

The script supports running multiple input cases in sequence.

It is not the main simulation script, but rather a driver that prepares and initiates
GEOS-MPM runs.

SLURM Configuration
-------------------

The script is designed to run as a SLURM job:

.. code-block:: bash

   #!/bin/bash
   #SBATCH -t 00:30:00
   #SBATCH -N 1
   #SBATCH -A imcomp
   #SBATCH -p pdebug

- ``-t``: wall time limit (30 minutes)
- ``-N``: number of nodes
- ``-A``: account
- ``-p``: partition (debug queue)

This job prepares simulations rather than running the full GEOS solve.

Input File Specification
------------------------

The script processes one or more input files defined by:

.. code-block:: bash

   fileNames=(
       case1
       case2
   )

Each entry corresponds to an input file:

.. code-block:: text

   pfw_input_<fileName>.py

For example:

- ``case1`` → ``pfw_input_case1.py``
- ``case2`` → ``pfw_input_case2.py``

Directory Configuration
----------------------

Two directories must be specified:

.. code-block:: bash

   fileLocation='<PATH_TO_INPUT_FILES>'
   runLocation='<PATH_TO_RUN_LOCATION>'

- ``fileLocation``: location of input files and supporting scripts
- ``runLocation``: destination for simulation run directories

For each case, the script creates:

.. code-block:: text

   <runLocation>/<fileName>/

User-Specific Configuration
--------------------------

The script automatically detects the current user:

.. code-block:: bash

   userName="$(whoami)"

This is used to copy a user-specific configuration file:

.. code-block:: text

   userDefs_<username>.py

Execution Flow
--------------

For each ``fileName``, the script performs the following steps:

#. Validate that the file name is not empty.
#. Determine the number of processes to use.
#. Check whether the run directory already exists.
#. Optionally prompt for overwrite (interactive mode only).
#. Delete any existing run directory.
#. Create a new run directory.
#. Copy required input and support files.
#. Execute the particle file writer.

Directory Handling
------------------

If the run directory already exists:

- In interactive mode, the user is prompted before overwriting.
- In batch mode, the directory is overwritten without prompting.

The following command removes the existing directory:

.. code-block:: bash

   rm -rf $runLocation/$fileName/

.. warning::

   This operation permanently deletes all existing files in the run directory.

File Copy Operations
--------------------

The following files are copied into each run directory:

- ``pfw_input_<fileName>.py``
- ``particleFileWriter.py``
- ``pfw_check.py``
- ``pfw_geometryObjects.py``
- ``userDefs_<username>.py``

These files define the simulation input, preprocessing logic, geometry,
restart handling, and user-specific paths.

Execution of particleFileWriter.py
-----------------------------------

The script runs the particle file writer using either serial or parallel execution.

Serial execution:

.. code-block:: bash

   python3 particleFileWriter.py pfw_input_<fileName>

Parallel execution:

.. code-block:: bash

   srun -n <num_tasks> python3 particleFileWriter.py pfw_input_<fileName>

The number of tasks is determined by a command-line argument:

- No argument or ``1`` → serial execution
- Any other value → parallel execution with ``srun``

Purpose of particleFileWriter
-----------------------------

The ``particleFileWriter.py`` script:

- Reads the ``pfw_input`` file
- Generates particle distributions and geometry
- Writes GEOS-compatible input files (e.g., XML)
- Optionally creates and submits simulation batch jobs



:ref:`Back to top <Input Files MPM>`




.. raw:: html

   <hr style="height:4px; border:none; background-color:black;">

.. _userDefs:

Environment Set Up: userDefs_username.py
========================
Overview
--------

This file is a machine-specific Python configuration file used in GEOS-MPM workflows. It detects
the current computing environment and assigns runtime settings such as the ``geosx`` executable
path, default run directories, and allocation account.

The ``userDefs`` file:

- Detects the machine environment
- Sets the ``geosx`` executable path
- Defines default run and test directories
- Specifies the default bank or account

The script is intended to be copied and renamed for each user, allowing local configuration
without modifying the shared repository. It distinguishes between Lassen and other systems and
sets paths accordingly.

This file is not a simulation input or launcher script, but a support file used by preprocessing
and run scripts.

Using a separate ``userDefs_<username>.py`` file:

- Keeps machine-specific settings out of shared scripts
- Allows per-user configuration
- Improves portability across systems

It provides the machine-specific information required for preprocessing, job generation,
and execution scripts.



File Purpose
------------

The purpose of this file is to centralize local runtime configuration information that may vary
from one machine or user environment to another.

These settings are commonly used by other GEOS-MPM support scripts, such as preprocessing,
job submission, and restart utilities.

Template Usage
--------------

The file is intended to be used as a template.

A shared version should remain unmodified, while each user should create a personal copy named
according to their LC username:

.. code-block:: text

   userDefs_<username>.py

For example:

- ``userDefs_meitner78.py``
- ``userDefs_maric75.py``

This allows each user to define their own paths and defaults without committing personal changes
to the repository.

Platform Detection
------------------

The script imports Python's built-in ``platform`` module and checks the current hostname:

.. code-block:: python

   import platform
   lassen = 'lassen' in platform.node()

This creates a Boolean variable:

- ``True`` if the hostname contains ``lassen``
- ``False`` otherwise

This value is then used to select the appropriate machine-specific settings.

Machine-Specific Configuration
-----------------------------

The script uses a conditional block to define settings for two environments:

- Lassen
- non-Lassen systems

Lassen configuration:

.. code-block:: python

   if lassen:
     geosPath='/usr/WS1/meitner78/GEOS/build-lassen-gcc@8.3.1-relwithdebinfo/bin/geosx'
     testRunDirectory='/p/gpfs1/meitner78/geosxRuns/test/'
     defaultRunDirectory='/p/gpfs1/meitner78/geosxRuns/test/'
     defaultBank='cbronze'

Non-Lassen configuration:

.. code-block:: python

   else:
     geosPath='/usr/workspace/malic75/GEOS/build-quartz-gcc@12-release/bin/geosx'
     testRunDirectory='/p/lustre1/malic75/geosxRuns/test/'
     defaultRunDirectory='/p/lustre1/malic75/geosxRuns/test/'
     defaultBank='imcomp'

Defined Variables
-----------------

The script defines the following variables:

- ``geosPath``: full path to the ``geosx`` executable
- ``testRunDirectory``: default directory for test runs
- ``defaultRunDirectory``: default directory for general runs
- ``defaultBank``: default account or bank used in batch scripts

These variables allow the rest of the GEOS-MPM workflow to reference machine-specific values
without hardcoding them in multiple scripts.

Execution Logic
---------------

The script performs the following steps:

#. Imports the ``platform`` module.
#. Checks the hostname of the current machine.
#. Determines whether the current environment is Lassen.
#. Assigns the appropriate GEOS executable path.
#. Assigns the appropriate test and default run directories.
#. Assigns the default bank or account name.

Role in the Workflow
--------------------

This file is typically imported by other GEOS-MPM scripts that need access to local machine
settings.

For example, another script may use:

.. code-block:: python

   from userDefs_<username> import geosPath, defaultRunDirectory, defaultBank

This allows launch and preprocessing scripts to automatically use the correct executable,
filesystem path, and account for the current user and machine.



:ref:`Back to top <Input Files MPM>`



.. raw:: html

   <hr style="height:4px; border:none; background-color:black;">

.. _pfw_check:

Simulation Execution: pfw_check.py
===========================

Overview
--------

This script is a restart and job-monitoring utility for GEOS-MPM simulations. Its purpose is to
inspect SLURM output, determine whether a simulation completed successfully or needs to be
restarted, and automatically relaunch the job from the most recent restart file if necessary.

The ``pfw_check`` script:

- Monitors SLURM output files for job status
- Detects early termination due to time limits, node failure, or other issues
- Identifies the most recent restart file
- Generates and submits restart jobs
- Resubmits itself to continue monitoring until completion

It is not responsible for generating input files or launching initial simulations, but instead
provides automated restart and monitoring capabilities to ensure simulations run to completion
with minimal manual intervention.



Command-Line Usage
------------------

The script requires two arguments:

.. code-block:: bash

   python3 pfw_check.py <inputFile> <jobID>

- ``inputFile``: name of the ``pfw_input`` module (without ``.py``)
- ``jobID``: SLURM job ID of the simulation to monitor

If arguments are missing, the script will exit with a usage message.



User-Specific Configuration
--------------------------

The script automatically imports a user-specific configuration file:

.. code-block:: python

   userDefs_<username>.py

This file provides:

- ``geosPath``: path to the GEOS executable
- Default machine-specific settings

If the file is not found, the script will terminate with an error message.



Machine Detection
-----------------

The script detects the current machine using the hostname and assigns the number of cores per node:

.. code-block:: python

   machineList = {
     'lassen':44,
     'dane':112,
     'ruby':56,
     'rzhound':56,
     'tioga':64
   }

This information is used to compute the number of nodes required for the restart job:

.. code-block:: python

   mNodes = ceil(mCores / coresPerNode)



Input File Handling
------------------

The script imports the specified ``pfw_input`` file as a Python module:

.. code-block:: python

   job = importlib.import_module(inputFile)
   pfw = job.pfw

This allows access to simulation parameters such as:

- ``mCores``: number of MPI ranks
- ``mWallTime``: wall time limit
- ``mBank``: allocation account
- ``runDebug``: debug mode flag
- domain decomposition parameters (``xpar``, ``ypar``, ``zpar``)



Job Monitoring Logic
--------------------

The script reads the SLURM output file:

.. code-block:: text

   slurm-<jobID>.out

It scans the file in reverse to efficiently detect the final job status.

The following conditions are checked:

- ``Job complete`` → simulation finished successfully
- ``Job exited early`` → early termination
- ``TIME LIMIT`` → wall time exceeded
- ``NODE FAILURE`` → hardware failure

If none of these conditions are found, the job is assumed to have exited for an unknown reason.



Restart Detection
-----------------

If a restart is required, the script searches for the most recent restart directory:

.. code-block:: bash

   *_restart_########

It selects the latest valid restart file and verifies that it is usable.

If no restart file is found, the script exits with an error.



Restart Job Creation
--------------------

A new SLURM script is generated to restart the simulation:

.. code-block:: bash

   srun -n <mCores> <geosPath> -i <input.xml> -x <xpar> -y <ypar> -z <zpar> -r <restartFile>

The script sets:

- number of nodes (``mNodes``)
- number of cores (``mCores``)
- partition (``pbatch`` or ``pdebug``)
- allocation account (``mBank``)

The restart job is submitted using:

.. code-block:: bash

   sbatch <restart_script>



Self-Resubmission
-----------------

After submitting the restart job, the script creates and submits a follow-up job that reruns
``pfw_check.py`` after the restart completes:

.. code-block:: bash

   #SBATCH --dependency=afterany:<jobID>

This creates a loop where:

- the simulation runs
- ``pfw_check.py`` evaluates its status
- a restart is submitted if needed
- monitoring continues until completion



Execution Flow
--------------

The script performs the following steps:

#. Import user-specific configuration and simulation input.
#. Determine machine type and compute node allocation.
#. Read the SLURM output file in reverse.
#. Identify whether the job completed or failed.
#. If needed, locate the most recent restart file.
#. Generate a restart SLURM script.
#. Submit the restart job.
#. Submit a dependent job to rerun ``pfw_check.py``.




:ref:`Back to top <Input Files MPM>`