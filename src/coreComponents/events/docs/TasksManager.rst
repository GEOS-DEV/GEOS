.. _TasksManager:

Tasks Manager
=============

The GEOS tasks manager allows a user to specify tasks to be executed. These tasks are compatible targets for the :ref:`EventManager`.

The tasks manager is configured via the ``Tasks`` block in an input .xml file, i.e.:

.. code-block:: xml

   <Tasks>
     <PackCollection name="historyCollection" objectPath="nodeManager" fieldName="velocity" />
   </Tasks>


Tasks Manager Configuration
---------------------------

Task
***************************
The children of the Tasks block define different Tasks to be triggered by events specified in the :ref:`EventManager` during the execution of the simulation.
The available task types are summarized below; the complete list of attributes for each one is given in the :doc:`Datastructure Index </docs/sphinx/datastructure/CompleteXMLSchema>`.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Task
     - Purpose
   * - ``PackCollection``
     - Collect time-history data from a field (on all or a subset of objects) for output by a ``TimeHistory`` output (see below and :ref:`OutputTasks`).
   * - ``CompositionalMultiphaseStatistics``, ``SinglePhaseStatistics``, ``SolidMechanicsStatistics``, ``SourceFluxStatistics``
     - Compute and report statistics (min/average/max of pressure, phase volumes, displacements, injected/produced amounts, ...) in the log and/or CSV files (see :ref:`LogCsvOutputs`).
   * - ``*Initialization`` (e.g. ``SinglePhasePoromechanicsInitialization``, ``SolidMechanicsLagrangianFEMInitialization``, ``HydrofractureInitialization``)
     - Run a solver once to bring the model to equilibrium with the initial conditions before the main simulation starts (see :ref:`PoromechanicalInitialization`).
   * - ``SolidMechanicsStateReset``
     - Reset the displacements (and velocities) of a solid mechanics solver and optionally enable/disable inelastic behavior, typically after an initialization step.
   * - ``FieldApplicator``
     - Apply field specifications to elements created during the simulation (e.g. fracture elements created by the ``SurfaceGenerator``).
   * - ``CellToCellDataCollection``
     - Output cell-to-cell connection data (currently the transmissibilities) computed by a flow solver.
   * - ``PVTDriver``, ``RelpermDriver``, ``TriaxialDriver``
     - Run a constitutive model on its own, without a mesh, to test or calibrate it (see :ref:`PVTDriver`, :ref:`RelpermDriver`, :ref:`TriaxialDriver`).

.. include:: /docs/sphinx/datastructure/Tasks.rst
    :start-line: 3
    
PackCollection
***************************
The ``PackCollection`` Task is used to collect time history information from fields. Either the entire field or specified named sets of indices in the field can be collected.

.. include:: /docs/sphinx/datastructure/PackCollection.rst
    :start-line: 3

Note: The time history information collected via this task is buffered internally until it is output by a linked TimeHistory Output.


***************************
Triggering the Tasks
***************************
Tasks can be triggered using the :ref:`EventManager`.
Recurring tasks should use a ``<PeriodicEvent>`` and one-time tasks should use a ``<SoloEvent>``:

.. code-block:: xml

  <PeriodicEvent name="historyCollectEvent"
                 timeFrequency="1.0"
                 targetExactTimestep="1"
                 target="/Tasks/historyCollection" />

The keyword ``target`` has to match the ``name`` of a Task specified as a child of the ``<Tasks>`` block.

