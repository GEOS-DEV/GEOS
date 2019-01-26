################################################################################
Linear solver parameters
################################################################################

************
Introduction
************

Any physics solver relying on standard finite element and finite volume techniques requires the solution of algebraic linear systems, which are obtained upon linearization and discretization of the governing equations, of the form:

.. math::

  \mathsf{A} \mathsf{x} = \mathsf{b}

with a :math:`\mathsf{A}` a square sparse matrix, :math:`\mathsf{x}` the solution vector, and :math:`\mathsf{b}` the right-hand side.
For example, in a classical linear elastostatics problem :math:`\mathsf{A}` is the stiffness matrix, and :math:`\mathsf{x}` and :math:`\mathsf{b}` are the displacement and nodal force vectors, respectively.


This solution stage represents the most computationally expensive portion of a typical simulation.
Solution algorithms generally belong to two families of methods: direct methods and iterative methods.
In GEOSX both options are made available wrapping around well-established open-source linear algebra libraries, namely
`HYPRE <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__,
`PETSC <https://www.mcs.anl.gov/petsc/>`__,
`SuperLU <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__, and
`Trilinos <https://trilinos.github.io/>`__.

**************
Direct methods
**************

The major advantages are their reliability, robustness, and ease of implementation.
However, they have large memory requirements and exhibit poor scalability.  Direct methods should be used in a prototyping stage, for example when developing a new formulation or algorithm, when the dimension of the problem, namely the size of matrix :math:`\mathsf{A}`, is small.
Irrespective of the selected direct solver implementation, the stages can be idenitified:

(1) **Setup Stage**: the matrix is first analyzed and then factorized
(#) **Solve Stage**: the solution to the linear systems involving the factorized matrix is computed
(#) **Finalize Stage**: the systems involving the factorized matrix have been solved and the direct solver lifetime ends

The default option in GEOSX relies on `SuperLU <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__, a general purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations, that is called taking advantage of the interface provided in `HYPRE <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__.

******************
Iterative  methods
******************

Iterative methods are the method of choice for large, three‐dimensional applications.
However, an appropriate preconditioner is essential for enabling the solution and achieving competitive performance.

*******
Summary
*******

The following table summarizes the available input parameters for the linear solver.

+--------------------------+----------+-----------+------------------------------------------------+
| Parameter                | Type     | Default   | Definition                                     |
+==========================+==========+===========+================================================+
|                          |          |           |                                                |
| ``type``                 | string   | ``gmres`` | Linear solver type.                            |
|                          |          |           | Valide values of this parameter are:           |
|                          |          |           |                                                |
|                          |          |           | * ``direct``:                                  |
|                          |          |           |   Direct solver                                |
|                          |          |           |                                                |
|                          |          |           | * ``cg``:                                      |
|                          |          |           |   Conjugate gradient method.                   |
|                          |          |           |   The coefficient matrix :math:`\mathsf{A}`    |
|                          |          |           |   must be symmetric positive definite          |
|                          |          |           |                                                |
|                          |          |           | * ``gmres``:                                   |
|                          |          |           |   Generalized minimum residual method.         |
|                          |          |           |   The coefficient matrix :math:`\mathsf{A}`    |
|                          |          |           |   can be indefinite or nonsymmetric            |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+
|                          |          |           |                                                |
| ``tolerance``            | real     | 1e-6      | Relative convergence tolerance of the          |
|                          |          |           | iterative method when using ``cg``,            |
|                          |          |           | ``gmres``, or ``bicgstab``.                    |
|                          |          |           | If the method converges, the iterative         |
|                          |          |           | solution :math:`\mathsf{x}_k` is such that     |
|                          |          |           | the relative residual norm satisfies           |
|                          |          |           | :math:`|| \mathsf{b}` -                        |
|                          |          |           | :math:`\mathsf{A} \mathsf{x}_k ||_2` <         |
|                          |          |           | ``tolerance`` * :math:`|| \mathsf{b} ||_2`     |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+
|                          |          |           |                                                |
| ``useAdaptiveTolerance`` | bool     | false     | Use Eisenstat-Walker adaptive tolerance        |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+
|                          |          |           |                                                |
| ``maxIterations``        | int      | 300       | Maximum number of linear iterations of the     |
|                          |          |           | iterative method when using ``cg``, ``gmres``, |
|                          |          |           | or ``bicgstab``.                               |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+
|                          |          |           |                                                |
| ``useRowScaling``        | bool     | false     | Apply row scaling before solving               |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+
|                          |          |           |                                                |
| ``useRowColScaling``     | bool     | false     | Apply row & column scaling before solving      |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+
|                          |          |           |                                                |
| ``maxRestart``           | int      | int       | Restart basis after ``maxRestart`` iterations  |
|                          |          |           | of the iterative method when using ``gmres``   |
|                          |          |           |                                                |
+--------------------------+----------+-----------+------------------------------------------------+


.. warning::

   The preconditioner parameter list is incomplete

The following table summarizes the available input parameters for the preconditioner.

+--------------------------+----------+--------------+------------------------------------------------+
| Parameter                | Type     | Default      | Definition                                     |
+==========================+==========+==============+================================================+
|                          |          |              |                                                |
| ``type``                 | string   | ``ilut``     | Preconditioner type.                           |
|                          |          |              | Valide values of this parameter are:           |
|                          |          |              |                                                |
|                          |          |              | * ``none``:                                    |
|                          |          |              |   No preconditioner                            |
|                          |          |              |                                                |
|                          |          |              | * ``ilut``:                                    |
|                          |          |              |   Dual threshold incomplete LU factorization   |
|                          |          |              |                                                |
|                          |          |              | * ``amg``:                                     |
|                          |          |              |   Algebraic multigrid                          |
|                          |          |              |                                                |
|                          |          |              | * ``ai``:                                      |
|                          |          |              |   Approximate inverse                          |
|                          |          |              |                                                |
|                          |          |              | * ``userDefined``:                             |
|                          |          |              |   The user is responsible for providing: (i) a |
|                          |          |              |   setup method that constructs the             |
|                          |          |              |   preconditioner,  and (ii) an apply method    |
|                          |          |              |   that performs the action of the              |
|                          |          |              |   preconditioning operator on a given vector   |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``ilut.fill``            | int      | 0            | Sparsity pattern fill factor for ``ilut``      |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``ilut.threshold``       | real     | 0.0          | Dropping threshold for ``ilut``                |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.maxLevels``        | int      | 20           | Max number of multigrid levels                 |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.cycleType``        | string   | ``V``        | Multigrid cycle type.                          |
|                          |          |              | Valide values of this parameter are:           |
|                          |          |              |                                                |
|                          |          |              | * ``V``:                                       |
|                          |          |              |   V-cycle                                      |
|                          |          |              |                                                |
|                          |          |              | * ``W``:                                       |
|                          |          |              |   W-cycle                                      |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.smootherType``     | string   | ``GS``       | Smoother to apply within AMG cycle             |
|                          |          |              | Valide values of this parameter are:           |
|                          |          |              |                                                |
|                          |          |              | * ``J``:                                       |
|                          |          |              |   Jacobi                                       |
|                          |          |              |                                                |
|                          |          |              | * ``GS``:                                      |
|                          |          |              |   Gauss-Seidel                                 |
|                          |          |              |                                                |
|                          |          |              | * ``ilut``:                                    |
|                          |          |              |   Dual threshold incomplete LU factorization   |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.coarseType``       | string   | ``direct``   | Solver used on the coarsest level              |
|                          |          |              | Valide values of this parameter are:           |
|                          |          |              |                                                |
|                          |          |              | * ``direct``:                                  |
|                          |          |              |   Direct solver                                |
|                          |          |              |                                                |
|                          |          |              | * ``smoother``:                                |
|                          |          |              |   Smoother used as iterative solver            |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.numSweeps``        | int      | 2            | Number of smoothing sweeps                     |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.symmetricProblem`` | bool     | true         | Optimizing setting for symmetric or            |
|                          |          |              | nonsymmetric problem                           |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
|                          |          |              |                                                |
| ``amg.nullSpaceType``    | string   | ``constant`` | Null space to use.                             |
|                          |          |              | Valide values of this parameter are:           |
|                          |          |              |                                                |
|                          |          |              | * ``constant``:                                |
|                          |          |              |   Constant null space vector                   |
|                          |          |              |                                                |
|                          |          |              | * ``RB``:                                      |
|                          |          |              |   Rigid body modes                             |
|                          |          |              |                                                |
+--------------------------+----------+--------------+------------------------------------------------+
