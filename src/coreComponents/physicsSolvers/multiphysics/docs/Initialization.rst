.. _PoromechanicalInitialization:

Poromechanical Models Initialization
===========================================

This section demonstrates the initialization of poromechanical models in GEOS: before the main simulation starts,
an initialization task (e.g. ``SinglePhasePoromechanicsInitialization``, see :ref:`TasksManager`) runs the coupled solver once
to compute the stress field in equilibrium with the initial pressure and the gravity or user-defined stress gradients.


.. toctree::
   :maxdepth: 1

   gravityInducedStressInitialization/Example

   userTableStressInitialization/Example
