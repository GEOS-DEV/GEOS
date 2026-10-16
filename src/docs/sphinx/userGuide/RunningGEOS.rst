.. _RunningGEOS:

###############################################################################
Running GEOS
###############################################################################

GEOS is a command-line program. It reads an XML input file (see :ref:`InputFileAnatomy`),
runs the simulation and writes its results to disk.

.. note::
   After ``make install`` the executable is available as ``geos`` in the ``bin`` directory
   of the install tree. In the build tree it is called ``bin/geosx``; both names refer to the
   same program. This guide uses ``geos``.


Serial run
==========

.. code-block:: sh

   geos -i inputFiles/singlePhaseFlow/3D_10x10x10_compressible_smoke.xml

Paths inside the XML file (included files, meshes, tables) are resolved relative to the
location of the XML file, not to the directory you run from.
Output files are written to the current directory unless ``-o`` is given.
The base name of output files defaults to the input file name without ``.xml``,
and can be changed with ``-n``.


Parallel run
============

GEOS uses MPI. Launch it with the MPI launcher of your system (``mpirun``, ``mpiexec``, ``srun``, ...):

.. code-block:: sh

   mpirun -np 8 geos -i input.xml -x 2 -y 2 -z 2

How the mesh is split across MPI ranks depends on the mesh type:

* ``InternalMesh`` (and ``InternalWellbore``) is split into a Cartesian grid of
  ``x`` × ``y`` × ``z`` partitions given by ``-x``, ``-y`` and ``-z``. **Their product must equal the
  number of MPI ranks**, otherwise GEOS stops with
  *"The total number of processes = ... does not correspond to the total number of partitions = ..."*.
  Each direction must also have at least as many cells as partitions.
  Choose the split so that each partition is roughly cube-shaped.
* ``VTKMesh`` is partitioned automatically by a graph partitioner (``partitionMethod`` attribute of
  ``VTKMesh``, ParMETIS by default); ``-x``, ``-y``, ``-z`` are not needed.

See :ref:`SpatialPartitioning` for details.


Command-line options
====================

``geos --help`` prints the list below.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Description
   * - ``-?``, ``--help``
     - Print the usage message.
   * - ``-i``, ``--input <file>``
     - Input XML file (required). May be repeated; the files are merged.
   * - ``-r``, ``--restart <file>``
     - Restart from a restart file (see below).
   * - ``-x``, ``--xpartitions <n>``
     - Number of partitions in the x-direction.
   * - ``-y``, ``--ypartitions <n>``
     - Number of partitions in the y-direction.
   * - ``-z``, ``--zpartitions <n>``
     - Number of partitions in the z-direction.
   * - ``-s``, ``--schema <file>``
     - Write the XML schema of this build to ``<file>`` and exit (no ``-i`` needed).
   * - ``-v``, ``--validate-input``
     - Read and initialize the input, then exit without running the simulation.
   * - ``-b``, ``--use-nonblocking``
     - Use non-blocking MPI communication.
   * - ``-n``, ``--name <name>``
     - Problem name, used as the base name of output files. Defaults to the first input file name without ``.xml``.
   * - ``--suppress-pinned``
     - Do not use pinned memory for MPI communication buffers.
   * - ``-o``, ``--output <dir>``
     - Directory where output files are written (default: current directory).
   * - ``-t``, ``--timers <string>``
     - Caliper timer configuration, e.g. ``-t runtime-report``. See :doc:`/docs/sphinx/developerGuide/Contributing/Caliper`.
   * - ``--trace-data-migration``
     - Trace host-device data migration (GPU builds).
   * - ``-m``, ``--memory-usage <threshold>``
     - Print memory allocations of data repository members above the threshold.
   * - ``--pause-for <seconds>``
     - Pause before starting execution (e.g. to attach a debugger).
   * - ``-e``, ``--errorsOutput <file>``
     - Also write errors and warnings to a file (``.yaml`` supported).


Checking an input file
======================

``geos -i input.xml -v`` reads the input, builds the mesh and initializes all objects,
then exits. Use it to catch errors in a new input file before launching a long run.


Restarting a simulation
=======================

To make a simulation restartable, add a ``Restart`` output and trigger it from an event:

.. code-block:: xml

   <Events maxTime="1.0e7">
     <!-- solver and other output events ... -->
     <PeriodicEvent name="restarts" timeFrequency="1.0e6" target="/Outputs/restartOutput"/>
   </Events>

   <Outputs>
     <Restart name="restartOutput"/>
   </Outputs>

Each time the event fires, GEOS writes a file ``<problemName>_restart_<cycle>.root``
(``<cycle>`` has 9 digits) plus a directory of the same name containing one HDF5 file per MPI rank.
To restart, run the same input file and pass the restart name **without** the ``.root`` extension:

.. code-block:: sh

   mpirun -np 4 geos -i input.xml -r input_restart_000000024 -x 2 -y 1 -z 2

The simulation resumes at the time of the checkpoint. Use the same number of MPI ranks as the run
that wrote the restart, because each rank reads its own file.
An example is given in :ref:`TutorialCO2FieldCase`.


Reading the log
===============

GEOS writes its progress to the standard output: the input files read, a summary of the mesh
and regions, then the time and step size of every time step. If a step does not converge,
it is cut and retried (see :ref:`SolutionStrategy`).

The amount of information is controlled per object by its ``logLevel`` attribute
(``0``, the default, prints the least; higher levels add to lower ones). For example,
``logLevel="1"`` on a flow solver adds the Newton residuals, linear solver iterations and
time-step information. The meaning of each level is listed in the ``logLevel``
description of each element in the :doc:`/docs/sphinx/datastructure/CompleteXMLSchema`.
Statistics written to the log and to CSV files are described in :ref:`LogCsvOutputs`.


Common first-run problems
=========================

``No XML input file nor schema specified``
   ``-i`` is missing.

``The total number of processes ... does not correspond to the total number of partitions``
   With an ``InternalMesh``, the product of ``-x``, ``-y`` and ``-z`` must equal the number of MPI ranks.

An XML attribute or element is rejected
   Attribute and element names are case-sensitive and arrays are written in braces, e.g. ``{ 1, 2, 3 }``.
   Check the spelling against the :doc:`/docs/sphinx/datastructure/CompleteXMLSchema`.

A file referenced in the XML is not found
   Relative paths are resolved from the directory of the XML file that contains them.

The time step keeps being cut
   The nonlinear solver does not converge. Increase the solver ``logLevel`` to see the residuals,
   and see :ref:`SolutionStrategy` and :ref:`LinearSolvers`.
