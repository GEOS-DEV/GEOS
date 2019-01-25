#######################################
Compositional multiphase flow FV solver
#######################################

Overview
=========================

Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/CompositionalMultiphaseFlow.rst

Input example
=========================

.. code-block:: xml

  <Solvers
    gravityVector="0.0,0.0,-9.81">

    <CompositionalMultiphaseFlow name="compflow"
                                 verboseLevel="1"
                                 gravityFlag="1"
                                 discretization="fluidTPFA"
                                 fluidName="fluid1"
                                 solidName="rock"
                                 relPermName="relperm"
                                 temperature="297.15"
                                 useMass="0"
                                 targetRegions="Region2">
      <SystemSolverParameters name="SystemSolverParameters"
                              krylovTol="1.0e-10"
                              newtonTol="1.0e-6"
                              maxIterNewton="15"
                              useDirectSolver="1"
                              solverType="Klu"
                              ilut_fill="0"
                              ilut_drop="0"/>
    </CompositionalMultiphaseFlow>
  </Solvers>