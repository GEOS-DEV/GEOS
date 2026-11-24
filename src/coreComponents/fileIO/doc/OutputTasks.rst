.. _OutputTasks:

============
Outputs
============

This section describes how outputs are handled by GEOS.

The outputs are defined in an ``<Outputs>`` XML block. Each output is triggered by an event (see below).
All output files are written in the output directory, which is the current directory unless specified
with the ``-o`` command-line option. The following output types are available:

=================== ================================================================================================
Output              Purpose
=================== ================================================================================================
``VTK``             Plot files for ParaView_ or VisIT_ (recommended for visualization)
``Silo``            Plot files in the SILO_ format, mostly for VisIT_
``TimeHistory``     Time series of selected fields, collected by ``PackCollection`` tasks, written to HDF5_ files
``Restart``         Restart files, used to resume a simulation (see :ref:`RestartOutput`)
``Blueprint``       Conduit Blueprint plot files
``ChomboIO``        Coupling files for the Chombo code
``MemoryStats``     Memory usage statistics (log and CSV)
``Python``          No actual output: returns control to the Python interpreter when GEOS is driven by :ref:`pygeosxInterface`
=================== ================================================================================================

************************
Defining an output
************************

SILO Output
===========

The SILO output is defined through the ``<Silo>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml 

  <Outputs>
    <Silo name="siloOutput"/>
  </Outputs>
  
  
The parameter options are listed in the following table:

.. include:: /docs/sphinx/datastructure/Silo.rst


VTK Output
===========

The VTK output is defined through the ``<VTK>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml 

  <Outputs>
    <VTK name="vtkOutput"/>
  </Outputs>

The parameter options are listed in the following table:

.. include:: /docs/sphinx/datastructure/VTK.rst

Useful attributes include:

* ``plotLevel``: fields are written only if their plot level is lower or equal to this value (default: 1); increase it to output more fields.
* ``fieldNames`` and ``onlyPlotSpecifiedFieldNames``: list fields to output in addition to the default ones; with ``onlyPlotSpecifiedFieldNames="1"``, only the listed fields are written, which reduces the size of the output.
* ``outputRegionType``: output only the ``cell``, ``well``, ``surface`` or ``particle`` regions (default: ``all``).
* ``writeGhostCells``: also write ghost cells (for debugging parallel runs).
* ``format``: ``binary`` (default) or ``ascii``.

TimeHistory Output
==================

The TimeHistory output is defined through the ``<TimeHistory>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml

  <Outputs>
    <TimeHistory name="timeHistoryOutput" sources="{/Tasks/collectionTask}" filename="timeHistory" />
  </Outputs>

The parameter options are listed in the following table:

.. include:: /docs/sphinx/datastructure/TimeHistory.rst

In order to properly collect and output time history information the following steps must be accomplished:

#. Specify one or more collection tasks using the :ref:`TasksManager`.
#. Specify a `TimeHistory Output` using the collection task(s) as source(s).
#. Specify an event in the :ref:`EventManager` to trigger the collection task(s).
#. Specify an event in the :ref:`EventManager` to trigger the output.

Note: Currently if the collection and output events are triggered at the same simulation time, the one specified first will also trigger first. Thus in order to output time history for the current time in this case, always specify the time history collection events prior to the time history output events.

.. _RestartOutput:

Restart Output
==============

The ``<Restart>`` XML node writes restart files, which contain the complete state of the simulation:

.. code-block:: xml

  <Outputs>
    <Restart name="restartOutput"/>
  </Outputs>

Each time the output is triggered, GEOS writes a file named ``<problemName>_restart_<cycle>.root`` and a directory with
the same name containing one HDF5 file per MPI rank. The problem name is the name of the input file without its ``.xml``
extension, unless specified with the ``-n`` command-line option.

To resume a simulation, pass the restart file **without** its ``.root`` extension to the ``-r`` command-line option,
together with the original input file, and run with the same number of MPI ranks and the same partitioning:

.. code-block:: sh

  mpirun -np 4 geos -i myProblem.xml -x 2 -y 2 -r myProblem_restart_000000100

The simulation then resumes exactly where it left off in the event loop (see :ref:`EventManager`).

************************
Triggering the outputs
************************

The outputs can be triggered using the :ref:`EventManager`.
It is recommended to use a ``<PeriodicEvent>`` to output results with a defined frequency:

.. code-block:: xml

  <PeriodicEvent name="outputs"
                 timeFrequency="5000.0"
                 targetExactTimestep="1"
                 target="/Outputs/siloOutput" />

The keyword ``target`` has to match with the path of the output node, i.e. ``/Outputs/`` followed by its ``name``.

****************************
Visualisation of the outputs
****************************

We suggest the use of Paraview_ (or VisIT_) to visualize the VTK outputs, and MatPlotLib_ to plot the TimeHistory outputs.

Visualizing Silo outputs with VisIT
===================================

If the ``<Silo>`` XML node was defined, GEOS writes the results in a folder called ``siloFiles``.

In VisIT :

1. File > Open file...
2. On the right panel, browse to the ``siloFiles`` folder.
3. On the left panel, select the file(s) you want to visualize. Usually, one file is written according the
   frequency defined in the ``timeFrequency`` keyword of the Event that has triggered the output.
4. To load fields, use the "Add" button and browse to the fields you want to plot.
5. To plot fields, use the "Draw" button.

Please consult the VisIT_ documentation for further explanations on its usage.


Visualizing VTK outputs with VisIT
===================================

If the ``<VTK>`` XML node was defined, GEOS writes the results in a folder named after the ``plotFileRoot`` attribute (by default, the ``name`` of the ``<VTK>`` node).
For problems with multiple active regions (e.g. :ref:`TutorialHydraulicFractureWithAdvancedXML`), additional work may be required to ensure that vtk files can be read by VisIt.
One option is to use the ``outputRegionType`` attribute to output separate sets of files per region.  For example:

.. code-block:: xml

  <Problem>
    <Events>
      <!-- Use nested events to trigger both vtk outputs -->
      <PeriodicEvent
        name="outputs"
        timeFrequency="2.0"
        targetExactTimestep="0">
        <PeriodicEvent
          name="outputs_cell"
          target="/Outputs/vtkOutput_cell"/>
        <PeriodicEvent
          name="outputs_surface"
          target="/Outputs/vtkOutput_surface"/>
      </PeriodicEvent>
    </Events>
    
    <Outputs>
      <VTK
        name="vtkOutput_cell"
        outputRegionType="cell"
        plotFileRoot="vtk_cell"/>

      <VTK
        name="vtkOutput_surface"
        outputRegionType="surface"
        plotFileRoot="vtk_surface"/>
    </Outputs>
  </Problem>


In VisIt, the prepared files can be visualized as follows:

1. Click the Open icon (or select File/Open from the top menu).
2. Enter the vtk output folder in the left-hand panel (these will be separate if you use the ``outputRegionType`` option above).
3. Select the desired ``*.vtm`` group in the right-hand panel (by double-clicking or selecting OK).
4. (Repeat for any additional datasets/regions.)
5. In the main window, select the desired dataset from the Active Sources dropdown menu.
6. Click on the Add icon, and select a parameter to visualize from the mesh, pseudocolor, or vector plot options
7. At some point you may be prompted to create a Database Correlation.  If you select Yes, then VisIt will create a new time slider, which can be used to change the time state while keeping the various aspects of the simulation synchronized.
8. Finally, click the Draw button and use the time slider to visualize the output.

Please consult the VisIT_ documentation for further explanations on its usage.


Visualizing VTK outputs with Paraview
=====================================

If the ``<VTK>`` XML node was defined, GEOS writes a folder and a ``.pvd`` file named after the ``plotFileRoot``
attribute (by default, the ``name`` of the ``<VTK>`` node).

The ``.pvd`` file references one dataset per output time, written according to the frequency of the Event that has
triggered the output. Each dataset contains one ``.vtu`` file per MPI process and per region.

To visualize the results:

1. Open the ``.pvd`` file in ParaView (File > Open) and click *Apply*. All output times are loaded, and can be browsed with the time controls.
2. Select the field to display in the coloring drop-down menu (e.g. ``pressure``, ``displacement``).
3. Use the *Block Selector* or the *Extract Block* filter to show only some regions (e.g. the reservoir or the wells).
4. Use the *Clip* or *Slice* filters to look inside a 3D model, and *Plot Over Line* or *Plot Data Over Time* to extract curves.

If a field is missing, increase ``plotLevel`` or add it to ``fieldNames`` (see above).

Visualizing TimeHistory outputs with MatPlotLib
===============================================

If the ``<TimeHistory>`` XML node was defined, GEOS writes a file named after the string defined
in the ``filename`` keyword and formatted as specified by the string defined in the ``format``
keyword (only HDF5_ is currently supported).

The TimeHistory file contains the collected time history information from each specified time history collector.
This information includes datasets for the time itself, any metadata sets describing index association with specified
collection sets, and the time history information itself.

It is recommended to use MatPlotLib_ and format-specific accessors (like H5PY for HDF5_) to access and easily plot the
time history data.

.. _SILO: https://wci.llnl.gov/simulation/computer-codes/silo
.. _VTK: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5
.. _VisIT: https://wci.llnl.gov/simulation/computer-codes/visit/downloads
.. _Paraview: https://www.paraview.org/
.. _MatPlotLib: https://matplotlib.org/
