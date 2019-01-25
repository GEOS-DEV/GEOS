##########################
Coupled poroelastic solver
##########################

Poroelasticity model
===========================================

This section describes the use of the poroelasticity models implemented in GEOSX.

Currently, the poroelasticity model is implemented using a "fixed-stress" coupling, which is a sequential method that solves the flow
problem first and then geomechanics problem with the rate of volumetric (mean) total stress variation frozen. 

In our model, the geomechanics (elasticity) equation is expressed in terms of the total stress :math:`\mathbf{\sigma}`:

.. math::
    \nabla \mathbf{\sigma} + \rho_b \mathbf{g} = 0

where it relates to effective stress :math:`\mathbf{\sigma\prime}` and pore pressure :math:`p` through Biot's coefficient :math:`b`:

.. math::
    \mathbf{\sigma} = \mathbf{\sigma\prime} - b p\mathbf{I}

The fluid mass conservation equation is expressed in terms of pore pressure and volumetric (mean) total stress:

.. math::
    \left( \frac{1}{M} + \frac{b^2}{K_{dr}} \right) \frac{\partial p}{\partial t} + \frac{b}{K_{dr}} \frac{\partial \sigma_v}{\partial t} + \nabla \cdot \mathbf{v}_f = f

where :math:`M` is the Biot's modulus and :math:`K_{dr}` is the drained bulk modulus. 

Unlike the conventional reservoir model that uses Lagrange’s porosity, in the coupled geomechanics and flow model, Euler’s porosity :math:`\phi` is adopted so the porosity variation is derived as:

.. math::
    \partial \phi = \left( \frac{b-\phi}{K_s}\right) \partial p + \left( b-\phi \right) \partial \epsilon_v

where :math:`K_{s}` is the bulk modulus of the solid grain and :math:`\epsilon_v` is the volumetric strain.  

For more details, refer to `Kim(2010)
<https://pangea.stanford.edu/ERE/pdf/pereports/PhD/Kim10.pdf>`_.

Model parameters
===========================================

The poroelasticity model is implemented as a main solver listed in
``<Solvers>`` block of the input XML file that calls both SolidMechanics_LagrangianFEM and SinglePhaseFlow solvers. 
In the main solver, it requires the specification of solidSolverName, fluidSolverName, and couplingTypeOption.

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/Poroelastic.rst

* ``couplingTypeOption``: defines the coupling scheme. Currently, only "FixedStress" is implemented.

For the parameter setup of each individual solver, please refer to the guideline of the specific solver.

Input example
===========================================

.. code-block:: xml

  <Solvers
    gravityVector="0.0,0.0,-9.81">

    <Poroelastic name="poroSolve" 
                 solidSolverName="lagsolve"
                 fluidSolverName="SinglePhaseFlow"
                 couplingTypeOption="FixedStress"
                 verboseLevel="1"
                 targetRegions="Region2">
      <SystemSolverParameters name="SystemSolverParameters"
                              maxIterNewton="40"
                              verbosityFlag="2"/>
    </Poroelastic>

  
    <SolidMechanics_LagrangianFEM name="lagsolve" 
                                  timeIntegrationOption="QuasiStatic"
                                  verboseLevel="1"
                                  >
      <SystemSolverParameters name="SystemSolverParameters"
                              krylovTol="1.0e-10"
                              newtonTol="1.0e-3"
                              maxIterNewton="10"
                              verbosityFlag="2"/>
    </SolidMechanics_LagrangianFEM>
  
    <SinglePhaseFlow name="SinglePhaseFlow"
                     verboseLevel="1"
                     gravityFlag="1"
                     discretization="singlePhaseTPFA"
                     fluidName="water"
                     solidName="shale">
      <SystemSolverParameters name="SystemSolverParameters"
                              krylovTol="1.0e-10"
                              newtonTol="1.0e-6"
                              maxIterNewton="8"/>
    </SinglePhaseFlow>
  </Solvers>