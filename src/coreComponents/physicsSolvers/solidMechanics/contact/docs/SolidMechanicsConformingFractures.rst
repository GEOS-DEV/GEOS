.. _SolidMechanicsConformingFractures:

###########################################
Solid mechanics conforming fractures solver
###########################################

Introduction
============

The conforming fractures solvers model the mechanical behavior of pre-existing fractures and faults
whose surfaces are explicitly represented in the mesh (i.e., the mesh conforms to the fracture surfaces).
The fracture surfaces are represented by a ``SurfaceElementRegion``, and the frictional contact between
the two faces of the fracture is enforced with Lagrange multipliers (tractions) as additional unknowns.
The contact behavior is described by a contact constitutive model, such as ``Coulomb``.
The governing equations are described in :ref:`ContactMechanicsSolver`.

Three variants are available:

* ``SolidMechanicsLagrangeContact``: Lagrange multipliers, with a stabilization built from a finite-volume stencil (``stabilizationName``),
* ``SolidMechanicsLagrangeContactBubbleStab``: Lagrange multipliers with a bubble-function stabilization,
* ``SolidMechanicsAugmentedLagrangianContact``: augmented Lagrangian formulation.

Worked examples are given in :ref:`ExamplesingleFracCompression`, :ref:`ExampleTFrac`, and :ref:`faultVerification`.

Parameters
=========================

The solver is specified as a child of the ``Solvers`` block, for example ``<SolidMechanicsLagrangeContact>``.
The following attributes are supported in the input block for `SolidMechanicsLagrangeContact`:

.. include:: /docs/sphinx/datastructure/SolidMechanicsLagrangeContact.rst

The following data are allocated and used by the solver:

.. include:: /docs/sphinx/datastructure/SolidMechanicsLagrangeContact_other.rst

