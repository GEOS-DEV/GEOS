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
Irrespective of the selected direct solver implementation (let's call it ``solver``), three methods are needed

* ``solverSetup``: analyze and factorize the matrix
* ``solverSolve``: compute the solution to the linear system
* ``solverDestroy``: called whenever all the systems involving its matrix have been solved and the solver lifetime ends

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

================================= ====== ========= =============================
Parameter                         Type   Default   Definition
--------------------------------- ------ --------- -----------------------------
``linearSolver.type``             string ``gmres`` Linear solver type.  Valid values  of this parameter are:

                                                   - ``direct``: Direct solver

                                                   - ``cg``: Conjugate gradient method.  The coefficient matrix :math:`\mathsf{A}` must be symmetric positive definite

                                                   - ``gmres``: Generalized minimum residual method.  The coefficient matrix :math:`\mathsf{A}` can be indefinite or nonsymmetric

                                                   - ``bicgstab``: Biconjugate gradient stabilized method.  The coefficient matrix :math:`\mathsf{A}` can be indefinite or nonsymmetric

``linearSolver.tolerance``        real   1e-6      Relative convergence tolerance of the iterative method when using ``cg``, ``gmres``, or ``bicgstab``. If the method converges, the iterative solution :math:`\mathsf{x}_k` is such that the relative residual norm satisfies :math:`|| \mathsf{b} - \mathsf{A} * \mathsf{x}_k ||_2` < ``linearSolver.tolerance`` * :math:`|| \mathsf{b} ||_2`

``linearSolver.maxIterations``    int    300       Maximum number of linear iterations of the iterative method when using ``cg``, ``gmres``, or ``bicgstab``.

``linearSolver.gmres.maxRestart`` int    300       Restart basis after ``maxRestart`` iterations of the iterative method when using ``gmres``
================================= ====== ========= =============================
