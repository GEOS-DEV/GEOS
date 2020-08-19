.. _LinearSolvers:

################################################################################
Linear Solvers
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

The major advantages are their reliability, robustness, and ease of use.
However, they have large memory requirements and exhibit poor scalability.  Direct methods should be used in a prototyping stage, for example when developing a new formulation or algorithm, when the dimension of the problem, namely the size of matrix :math:`\mathsf{A}`, is small.
Irrespective of the selected direct solver implementation, three stages can be idenitified:

(1) **Setup Stage**: the matrix is first analyzed and then factorized
(#) **Solve Stage**: the solution to the linear systems involving the factorized matrix is computed
(#) **Finalize Stage**: the systems involving the factorized matrix have been solved and the direct solver lifetime ends

The default option in GEOSX relies on `SuperLU <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__, a general purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations, that is called taking advantage of the interface provided in `HYPRE <https://computation.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__.

******************
Iterative  methods
******************

As the problem size (number of computational cells) increases, global iterative solution strategies are the method of choice---typically nonsymmetric Krylov solvers.
Because of the possible poor conditioning of :math:`\mathsf{A}`, preconditioning is essential to solve such systems efficiently.
*''Preconditioning is simply a means of transforming the original linear system into one which has the same solution, but which is likely to be easier to solve with an iterative solver''* [Saad (2003)].

The design of a robust and efficient preconditioner is based on a trade-off between two competing objectives:

* **Robustness**: reducing the number of iterations needed by the preconditioned solver to achieve convergence;
* **Efficiency**: limiting the time required to construct and apply the preconditioner.

Assuming a preconditioning matrix :math:`\mathsf{M}` is available, three standard approaches are used to apply the preconditioner:

(1) **Left preconditioning**: the preconditioned system is :math:`\mathsf{M}^{-1} \mathsf{A} \mathsf{x} = \mathsf{M}^{-1} \mathsf{b}`
(#) **Right preconditioning**: the preconditioned system is :math:`\mathsf{A} \mathsf{M}^{-1} \mathsf{y} = \mathsf{b}`, with :math:`\mathsf{x} = \mathsf{M}^{-1} \mathsf{y}`
(#) **Split preconditioning**: the preconditioned system is :math:`\mathsf{M}^{-1}_L \mathsf{A} \mathsf{M}^{-1}_R \mathsf{y} = \mathsf{M}^{-1}_L \mathsf{b}`, with :math:`\mathsf{x} = \mathsf{M}^{-1}_R \mathsf{y}`


*******
Summary
*******

The following table summarizes the available input parameters for the linear solver.

.. include:: ../../fileIO/schema/docs/LinearSolverParameters.rst
