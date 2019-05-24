.. _SinglePhaseFlow:

#####################################
Single phase flow FV solver
#####################################

Overview
=========================

This is a cell-centered Finite Volume solver for compressible single-phase flow in porous media.
Fluid pressure as the primary solution variable.
Darcy's law is used to calculate fluid velocity from pressure gradient.
The solver currently only supports Dirichlet-type boundary conditions (BC) applied on cells or faces and Neumann no-flow type BC.

The following mass balance equation is solved in the domain:

.. math::
   \frac{\partial}{\partial t}(\phi\rho) + \boldsymbol{\nabla} \cdot (\rho\boldsymbol{u}) + q = 0,

where

.. math::
   \boldsymbol{u} = -\frac{1}{\mu}\boldsymbol{k}(\nabla p - \rho \boldsymbol{g})

and :math:`\phi` is porosity, :math:`\rho` is fluid density, :math:`\mu` is fluid viscosity,
:math:`\boldsymbol{k}` is the permeability tensor, :math:`\boldsymbol{g}` is the gravity vector,
and :math:`q` is the source function (currently not supported). The details on the computation of the density and the viscosity are given in :ref:`CompressibleSinglePhaseFluid`.

When the entire pore space is filled by a single phase, we can substitute the Darcy's law into the mass balance equation to obtain the single phase flow equation

.. math::
   \frac{\partial}{\partial t}(\phi\rho) - \boldsymbol{\nabla} \cdot \frac{\rho \boldsymbol{k}}{\mu} (\nabla p - \gamma \nabla z) + q = 0,

with :math:`\gamma \nabla z= \rho \boldsymbol{g}`.


Space discretization with Finite volume mehtod
=========================
Let :math:`\Omega \subset \mathbb{R}^n, \, n =1,2,3` be an open set defining the computational domain. We consider :math:`\Omega` meshed by element such that :math:`\Omega = \cup_{i}V_i` and integrate the single phase flow equation, described above, over each element :math:`V_i`:

.. math::
   \int_{V_i} \frac{\partial}{\partial t}(\phi\rho) dV - \int_{V_i} \boldsymbol{\nabla} \cdot \frac{\rho \boldsymbol{k}}{\mu} (\nabla p - \gamma \nabla z) dV + \int_{V_i} q dV  = 0.

Applying the divergence theorem to the second term lead to

.. math::
  \int_{V_i} \frac{\partial}{\partial t}(\phi\rho)_i - \oint_{S_i} \nabla \left(\frac{\rho \boldsymbol{k}}{\mu}(\nabla p -\gamma \nabla z)\right) \cdot \boldsymbol{n} dS + \int_{V_i} q dV  = 0.

where :math:`S_i` represents the surface area of the element :math:`V_i` and :math:`\boldsymbol{n}` is a outward unit vector normal to the surface.

For the flux term, the (static) transmissibility is currently computed with a Two-Point Flux Approximation (TPFA) as described in :ref:`FiniteVolume`.

The pressure-dependent mobility :math:`\lambda = \frac{\rho}{\mu}` at the interface is approximated using a first-order upwinding on the sign of the potential difference.

Time discretization
=========================
Let :math:`t_0 < t_1 < \cdots < t_N=T` be a grid discretization of the time interval :math:`[t_0,T], \, t_0, T \in \mathbb{R}^+`. We use the backward Euler (fully implicit) method to integrate the single phase flow equation between two grid points :math:`t_n` and :math:`t_{n+1}, \, n< N` to obtain the residual equation:

.. math::
   \int_{V_i} \frac{(\phi\rho)_i^{n+1} - (\phi\rho)_i^n}{\Delta t} - \oint_{S_i} \nabla \left(\frac{\rho \boldsymbol{k}}{\mu}(\nabla p -\gamma \nabla z)\right)^{n+1} \cdot \boldsymbol{n} dS + \int_{V_i} q^{n+1} dV = 0

where :math:`\Delta t = t_{n+1}-t_n` is the time-step. The expression of this residual equation and its derivative are used to form a linear system, which is solved via the solver package.

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
must be prescribed on this field on cell or face sets of interest.

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
