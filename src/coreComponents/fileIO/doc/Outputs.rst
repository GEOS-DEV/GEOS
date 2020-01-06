:orphan:

============
Outputs
============

This section describes how outputs are handled by GEOSX

The outputs are defined in a ``<Outputs>`` XML block.

There are two available formats to output the results of a simulation: SILO_ and VTK_.

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

.. include:: /coreComponents/fileIO/schema/docs/Silo.rst


VTK Output
===========

The VTK output is defined through the ``<VTK>`` XML node (subnode of ``<Outputs> XML block``) as shown here:

.. code-block:: xml 

  <Outputs>
    <VTK name="vtkOutput"/>
  </Outputs>

The parameter options are listed in the following table:

.. include:: /coreComponents/fileIO/schema/docs/VTK.rst

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

The keyword ``target`` has to match with the name of the ``<Silo>`` or ``<VTK>`` node.

****************************
Visualisation of the outputs
****************************

We suggest the use of VisIT_ and Paraview_ to visualize the outputs.

Visualizing Silo outputs with VisIT
===================================

If the ``<Silo>`` XML node was defined, GEOSX writes the results in a folder called ``siloFiles``.

In VisIT : 

1. File > Open file...
2. On the right panel, browse to the ``siloFiles`` folder.
3. On the left panel, select the file(s) you want to visualize. Usually, one file is written according the
   frequency defined in the ``timeFrequency`` keyword of the Event that has triggered the output.
4. To load fields, use the "Add" button and browse to the fields you want to plot.
5. To plot fields, use the "Draw" button.

Please consult the VisIT_ documentation for further explanations on its usage.

Visualizing VTK outputs with Paraview
=====================================

If the ``<VTK>`` XML node was defined, GEOSX writes a folder and a ``.pvd`` file named after the string defined
in ``name`` keyword.

The ``.pvd`` file contains references to the ``.pvtu`` files. One ``.pvtu`` file is output according the frequency defined in the ``timeFrequency`` keyword of the Event that has triggered the output.

One ``.pvtu`` contains references to ``.vtu`` files. There is as much ``.vtu`` file as there were MPI processes
used for the computation.

All these files can be opened with paraview. To have the whole results for every output time steps, you can
open the ``.pvd`` file.

.. _SILO: https://wci.llnl.gov/simulation/computer-codes/silo
.. _VTK: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
.. _VisIT: https://wci.llnl.gov/simulation/computer-codes/visit/downloads
.. _Paraview: https://www.paraview.org/
