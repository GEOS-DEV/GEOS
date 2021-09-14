################################################################################
Constitutive models
################################################################################

In GEOSX, all constitituive models (e.g, those defining fluid and rock properties)
are implemented in the namespace ``constitutive`` and derive from a common base class,
``ConstitutiveBase``. All objects are owned and handled by the ``ConstitutiveManager``.

Standalone models
======================================================
Standalone models are use to implement constitutive laws such as

  - mechanical material models (i.e.g, linear elasticity, plasticity, etc.)
  - PVT fluid behaviors
  - relative permeability relationships
  - porosity and permeability dependencies on state variables
  - contact laws.

Storage, allocation and update of properties
-----------------------------------------------------
Each constitutive model owns, as member variables, ``LvArray::Array`` containers
that hold the properties (or fields) handled and, their derivatives w.r.t to the
fields needed to update each property. Each property is stored as an array in which the
first dimension represents the elementIndex and the second dimension the index of the
integration point. Obviously these dimensions are determined by the size of the grid
and by the type of of discretization method chosen. Vectorial and tensorial will
also have an additional dimension to identify their components. Similarly, for multiphase
fluid models an additional dimension will be added to properties defined for each phase
and or component.
As such, for example, a single phase fluid model in which density and viscosity are functions of the fluid pressure will have
the following members

.. literalinclude:: ../SingleFluidBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

Compound constitutive models
========================================================
