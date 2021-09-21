.. _SinglePhasePoroelasticSolver:

#################################
Single-Phase Poromechanics Solver
#################################

Introduction
===========================================

This section describes the single-phase poromechanics model, which is used to simulate fully coupled processes of three-dimensional deformation and single-phase fluid flow through the void space in a porous medium.
This single-phase poromechanics solver is an example of multiphysics solver obtained coupling a quasi-static solid mechanics solver with a single-phase flow solvers.
The governing equations are solved using a fully-implicit time-marching scheme (backward Euler).
The spatial discretization relies on a mixed formulation combing continuous Galerkin finite elements with either a two-point flux approximation finite volume methor or a hybrid mimetic finite difference method.

Theory
=========================

Governing Equations
--------------------------

Let :math:`\Omega` denote the three-dimensional domain occupied bu the porous mediu. 
Let :math:`t` denote time, which belongs to the time interval of interest :math:`\mathcal{I} = (0, t_{\max}]`.
The requisite set of governing equations consists of two conservation laws, one for the linear momentum balance and one for the mass balance of the fluid phase.
Assuming quasi-static mechanical behavior, the strong form of the single-phase poromechanical initial/boundary value problem (IBVP) consists of finding the solid displacement vector :math:`\boldsymbol{u}` and the fluid pore pressure :math:`p` such that

.. math::
     & -\nabla \cdot \boldsymbol{\sigma} = \rho \boldsymbol{g} 
     && \text{in } \Omega \times \mathcal{I}
     && \text{(quasi-static linear momentum balance)}

     & \dot{m}_f + \nabla \cdot \boldsymbol{w} = q
     && \text{in } \Omega \times \mathcal{I}
     && \text{(fluid phase mass balance)}

where 

- :math:`\boldsymbol{\sigma}  = ( \boldsymbol{\sigma}^{\prime} - b p \boldsymbol{1})` is the total Cauchy stress tensor, with :math:`\boldsymbol{\sigma}^{\prime} = \mathbb{C} : \nabla^s \boldsymbol{u}` the effective stress, :math:`\mathbb{C}` the rank-4 drained tensor of tangent moduli, :math:`b` the Biot coefficient, and :math:`\boldsymbol{1}` the rank-2 identity tensor;
- :math:`\rho \boldsymbol{g}` is a body force due to self-weigth of the mixture, with :math:`\rho = [(1 - \phi) \rho_s + \phi \rho_f]` the density of the mixture, :math:`\phi` the porosity, :math:`\rho_s` the density of the solid phase, :math:`\rho_f` the density of the fluid phase, and :math:`\boldsymbol{g}` the gravity vector;
- :math:`m_f = (\phi \rho_f)` is the fluid mass content, namely the fluid mass per unit reference volume;
- :math:`\boldsymbol{w} = - \rho_f \frac{1}{\mu} \boldsymbol{\kappa} \cdot \nabla \Phi_f` is th fluid phase mass flux, with :math:`\mu` the fluid viscosity, :math:`\boldsymbol{\kappa}` the absolute permeability tensor, and :math:`\Phi_f = ( p - \rho_f \boldsymbol{g} \cdot \boldsymbol{x})` the fluid phase potential, and :math:`\boldsymbol{x}` the position vector in :math:`\mathbb{R}^3`;
- :math:`q` is a mass source/sink per unit volum term for the fluid phase;
- :math:`\nabla`, :math:`\nabla^s`, and :math:`\nabla \cdot` are the gradient, symmetric gradient, and divergence operator, respectively;
- the superposed dot, :math:`\dot{(\bullet)}`, indicates the derivative of quantitiy :math:`(\bullet)` with respect to time.


Parameters
===========================================

The poroelasticity model is implemented as a main solver listed in
``<Solvers>`` block of the input XML file that calls both SolidMechanicsLagrangianSSLE and SinglePhaseFlow solvers.
In the main solver, it requires the specification of solidSolverName, fluidSolverName, and couplingTypeOption.

The following attributes are supported:

.. include:: /coreComponents/schema/docs/SinglePhasePoromechanics.rst

* ``couplingTypeOption``: defines the coupling scheme.

The solid constitutive model used here is PoroLinearElasticIsotropic, which derives from ElasticIsotropic and includes an additional parameter: Biot's coefficient. The fluid constitutive model is the same as SinglePhaseFlow solver. For the parameter setup of each individual solver, please refer to the guideline of the specific solver.

An example of a valid XML block for the constitutive model is given here:

.. literalinclude:: ../integratedTests/PoroElastic_Terzaghi_FIM.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_CONSTITUTIVE -->
  :end-before: <!-- SPHINX_POROELASTIC_CONSTITUTIVE_END -->

Example
===========================================

.. literalinclude:: ../integratedTests/PoroElastic_Terzaghi_FIM.xml
  :language: xml
  :start-after: <!-- SPHINX_POROELASTIC_SOLVER -->
  :end-before: <!-- SPHINX_POROELASTIC_SOLVER_END -->
