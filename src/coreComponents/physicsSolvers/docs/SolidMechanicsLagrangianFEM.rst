.. _SolidMechanicsLagrangianFEM:

#####################################
Lagrangian solid mechanics FEM solver
#####################################

Overview
=========================

Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/SolidMechanics_LagrangianFEM.rst

Input example
=========================

.. code-block:: xml

  <Solvers>
    <SolidMechanics_LagrangianFEM name="lagsolve"
                                  timeIntegrationOption="QuasiStatic"
                                  discretization="FE1"
                                  targetRegions="Region2">
      <SystemSolverParameters name="solverParams0"
                              useMLPrecond="1"
                              scalingOption="0"
                              krylovTol="1.0e-8"
                              newtonTol="1.0e-4"
                              maxIterNewton="8"
                              verbosityFlag="0"/>
     </SolidMechanics_LagrangianFEM>
  </Solvers>