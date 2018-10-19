########################
Solver Parameters
########################

The following tables provides a comprehensive list of solver parameters
controlling the linear and nonlinear behavior. To keep things organized,
some parameters are nested as ``topLevel.subLevel.subSubLevel``.

***************************
Nonlinear Solver Parameters
***************************

================================= ====== ========== =================================================================
Parameter                         Type   Default    Definition [Options]
================================= ====== ========== =================================================================
``nonlinearSolver.type``          string ``newton`` Nonlinear solver type ``[newton, backtrackingNewton, appleyard]``
``nonlinearSolver.tolerance``     real   1e-5       Relative convergence tolerance
``nonlinearSolver.maxIterations`` int    25         Max nonlinear iterations before timestep cut
================================= ====== ========== =================================================================

************************
Linear Solver Parameters
************************

========================================= ====== ========= ====================================================
Parameter                                 Type   Default   Definition ``[Options]``
========================================= ====== ========= ====================================================
``linearSolver.type``                     string ``gmres`` Linear solver type ``[direct, cg, gmres, bicgstab]``
``linearSolver.tolerance``                real   1e-6      Relative convergence tolerance
``linearSolver.useAdaptiveTolerance``     bool   false     Use Eisenstat-Walker adaptive tolerance
``linearSolver.maxIterations``            int    300       Max linear iterations
``linearSolver.scaling.useRowScaling``    bool   false     Apply row scaling before solving
``linearSolver.scaling.useRowColScaling`` bool   false     Apply row & column scaling before solving
``linearSolver.gmres.maxRestart``         int    300       Restart basis after ``maxRestart`` iterations
========================================= ====== ========= ====================================================

**************************
Preconditioner Parameters
**************************

======================================= ====== =============== ==================================================================================
Parameter                               Type   Default         Definition ``[Options]``
======================================= ====== =============== ==================================================================================
``preconditioner.type``                 string ``ilut``        Preconditioner type ``[none, iluk, ilut, amg, userDefined]``
``preconditioner.blockSize``            int    1               Some algorithms support block-variants for matrices with dense sub-blocks.
``preconditioner.ilu.fill``             int    0               Sparsity pattern fill factor for ILUK, ILUT
``preconditioner.ilu.threshold``        real   0.0             Dropping threshold for ILUT
``preconditioner.amg.maxLevels``        int    20              Max number of multigrid levels
``preconditioner.amg.cycleType``        string ``V``           Multigrid cycle type ``[V,W]``
``preconditioner.amg.smootherType``     string ``gaussSeidel`` Smoother to apply within AMG cycle ``[jacobi, gaussSeidel, chebyshev, ilut, ...]``
``preconditioner.amg.coarseType``       string ``direct``      Treatment of coarsest level solve ``[direct,smoother]``
``preconditioner.amg.numSweeps``        int    2               Number of smoothing sweeps
``preconditioner.amg.symmetricProblem`` bool   true            Optimize settings for symmetric or nonsymmetric problems
``preconditioner.amg.nullSpaceType``    string ``constant``    Null space to use ``[constant,rigidBody]``
``preconditioner.sai.placeHolder``      tbd    tbd             Placeholder anticipating Sparse Approximate Inverse preconditioners
======================================= ====== =============== ==================================================================================
