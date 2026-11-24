.. _UserGuide:

###############################################################################
User Guide
###############################################################################

This guide describes how to set up and run GEOS simulations.
If you are new to GEOS, start with the :ref:`Tutorials`, then read
:ref:`Capabilities`, :ref:`InputFileAnatomy` and :ref:`RunningGEOS`.
The remaining chapters follow the order in which an input file is usually written:
mesh, materials, physics and numerics, boundary conditions, then time stepping and outputs.
The complete list of XML elements and attributes is in the
:doc:`/docs/sphinx/datastructure/CompleteXMLSchema`.

.. toctree::
   :maxdepth: 2

   Capabilities

   InputFileAnatomy

   RunningGEOS

   /coreComponents/fileIO/doc/InputXMLFiles.rst

   /coreComponents/mesh/docs/Mesh

   /coreComponents/constitutive/docs/Constitutive

   /coreComponents/physicsSolvers/PhysicsSolvers

   /coreComponents/discretizationMethods/docs/NumericalMethodsManager

   /coreComponents/linearAlgebra/docs/LinearSolvers.rst

   /coreComponents/fieldSpecification/docs/FieldSpecification

   /coreComponents/functions/docs/FunctionManager

   /coreComponents/events/docs/EventManager

   /coreComponents/events/docs/TasksManager

   /coreComponents/physicsSolvers/multiphysics/docs/Initialization.rst

   /coreComponents/fileIO/doc/Index.rst

   /coreComponents/mesh/mpiCommunications/docs/SpatialPartition.rst

   /coreComponents/constitutiveDrivers/docs/ConstitutiveDrivers

   /pygeosx/pygeosx_documentation.rst
