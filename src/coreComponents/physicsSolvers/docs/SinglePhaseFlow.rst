#####################################
Single phase flow FV solver
#####################################

Overview
=========================

Usage
=========================

The following attributes are supported:

.. include:: /coreComponents/fileIO/schema/docs/SinglePhaseFlow.rst


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
