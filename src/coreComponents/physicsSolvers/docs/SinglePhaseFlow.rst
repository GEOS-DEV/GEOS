.. _SinglePhaseFlow:

#####################################
Single phase flow FV solver
#####################################

Overview
=========================

This is a cell-centered Finite Volume solver for compressible single-phase flow in porous media.
Fluid pressure as the primary solution variable.
Darcy's law is used to calculate fluid velocity from pressure gradient.
The solver currently only supports Dirichlet-type boundary conditions applied on cells or faces.

The following mass balance equation is solved in the domain:

.. math::
   \frac{\partial}{\partial t}(\phi\rho) + \boldsymbol{\nabla} \cdot (\rho\boldsymbol{u}) + q = 0

where

.. math::
   \boldsymbol{u} = -\frac{1}{\mu}\boldsymbol{k}(\nabla p - \rho \boldsymbol{g})

and :math:`\phi` is porosity, :math:`\rho` is fluid density, :math:`\mu` is fluid viscosity,
:math:`\boldsymbol{k}` is the permeability tensor, :math:`\boldsymbol{g}` is the gravity vector,
and :math:`q` is the source function (currently not supported).

Usage
=========================

The solver is enabled by adding a ``<SinglePhaseFlow>`` node in the Solvers section.
Like any solver, time stepping is driven by events, see :ref:`EventManager`.

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/SinglePhaseFlow.rst

In particular:

* ``discretization`` must point to a Finite Volume flux approximation scheme defined in the Numerical Methods section of the input file (see :ref:`FiniteVolume`)
* ``fluidName`` must point to a single phase fluid model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``solidName`` must point to a solid mechanics model defined in the Constitutive section of the input file (see :ref:`Constitutive`)
* ``targetRegions`` attribute is currently not supported, the solver is always applied to all regions.

Primary solution field label is ``pressure``.
Initial conditions must be prescribed on this field in every region, and boundary conditions
must be presribed on this field on cell or face sets of interest.

In addition, the solver declares a scalar field named ``referencePorosity`` and a vector field
named ``permeability``, that contains principal values of the symmetric rank-2 permeability tensor
(tensor axis are assumed aligned with the global coordinate system).
These fields must be populated via :ref:`FieldSpecification` section and ``permeability`` should
be supplied as the value of ``coefficientName`` attribute of the flux approximation scheme used.


Input example
=========================

.. code-block:: xml

  <Solvers
    gravityVector="0.0,0.0,-9.81">

    <SinglePhaseFlow name="SinglePhaseFlow"
                     verboseLevel="3"
                     gravityFlag="1"
                     discretization="singlePhaseTPFA"
                     fluidName="water"
                     solidName="rock"
                     targetRegions="Region2">
      <SystemSolverParameters name="SystemSolverParameters"
                              krylovTol="1.0e-10"
                              newtonTol="1.0e-6"
                              maxIterNewton="8"/>
    </SinglePhaseFlow>
  </Solvers>
