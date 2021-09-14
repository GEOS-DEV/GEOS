.. _constitutiveModels:

################################################################################
Constitutive models in GEOSX
################################################################################

In GEOSX, all constitituive models (e.g, those defining fluid and rock properties)
are implemented in the namespace ``constitutive`` and derive from a common base class,
``ConstitutiveBase``. All objects are owned and handled by the ``ConstitutiveManager``.

Standalone models
======================================================
Standalone constitutive models implement constitutive laws such as

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
integration point. These dimensions are determined by the number of elemetns of the
subregion on which each constitutive model is registered, and by the type of discretization
method chosen. Additionally, vectorial and tensorial will also have an additional dimension to identify
their components. Similarly, for multiphase fluid models an additional dimension will be
added to properties defined for each phase and/or component.
Thus, for example, a single phase fluid model in which density and viscosity are
functions of the fluid pressure will have the following members

.. literalinclude:: /coreComponents/constitutive/fluid/SingleFluidBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

Resizing of all fields of the constitutive models is performed during the initialization phase by
the ``ConstitutiveManger`` through the call ``ConstitutiveManger::hangConstitutiveRelation``
which sets the appropriate subRegion as the parent Group of each constitutive model object.
Additionally, it also resizes all fields based on the size of the subregion and on the number of quadrature
points on it, by calling ``CONSTITUTIVE_MODEL::allocateConstitutiveData``. For the
single phase fluid example used before this call is

.. literalinclude:: /coreComponents/constitutive/fluid/SingleFluidBase.cpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_00
   :end-before: //END_SPHINX_INCLUDE_00

Any property (or field) stored on a constitutive model must be updated within a computational
kernel to ensure that, in GPU runs, `host` and `device` memory are properly synced and to ensure that the
updates are performed on `device`. Additionally, some properties are may be updated
within finite element kernels of specific physics (e.g., stress in a mechanics kernel). Consequently,
for each constitutive model class, a corresponding `nameOfTheModelUpdates` which, since it only contains
``LvArray::arrayView`` containers to the data, can be capture by value inside compuational kernels.
For example, for the single phase fluid model the 

.. literalinclude:: /coreComponents/constitutive/fluid/SingleFluidBase.hpp
   :language: c++
   :start-after: //START_SPHINX_INCLUDE_01
   :end-before: //END_SPHINX_INCLUDE_01

Compound models
========================================================
