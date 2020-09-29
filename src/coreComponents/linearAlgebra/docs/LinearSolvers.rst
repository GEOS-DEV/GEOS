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
`PETSc <https://www.mcs.anl.gov/petsc/>`__,
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

***************************
Preconditioner descriptions
***************************

This section provides a brief description of the available preconditioners.

* **None**: no preconditioning is used, i.e., :math:`\mathsf{M}^{-1} = \mathsf{I}`.
* **Jacobi**: diagonal scaling preconditioning, with :math:`\mathsf{M}^{-1} = \mathsf{D}^{-1}`, with :math:`\mathsf{D}` the matrix diagonal.
  Further details can be found in:

  - `HYPRE documentation <https://hypre.readthedocs.io/en/latest/api-sol-parcsr.html>`__,
  - `PETSc documentation <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCJACOBI.html>`__,
  - `Trilinos documentation <https://docs.trilinos.org/dev/packages/ifpack/doc/html/classIfpack__PointRelaxation.html>`__.

* **ILUK**: incomplete LU factorization with fill level k of the original matrix: :math:`\mathsf{M}^{-1} = \mathsf{U}^{-1} \mathsf{L}^{-1}`.
  Further details can be found in:

  - `HYPRE documentation <https://hypre.readthedocs.io/en/latest/solvers-hypre-ilu.html>`__,
  - `PETSc documentation <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCILU.html>`__,
  - `Trilinos documentation <https://docs.trilinos.org/dev/packages/ifpack/doc/html/classIfpack__ILU.html>`__.

* **ILUT**: a dual threshold incomplete LU factorization: :math:`\mathsf{M}^{-1} = \mathsf{U}^{-1} \mathsf{L}^{-1}`.
  Further details can be found in:

  - `HYPRE documentation <https://hypre.readthedocs.io/en/latest/solvers-hypre-ilu.html>`__,
  - not yet available through PETSc interface,
  - `Trilinos documentation <https://docs.trilinos.org/dev/packages/ifpack/doc/html/classIfpack__ILUT.html>`__.

* **ICC**: incomplete Cholesky factorization of a symmetric positive definite matrix: :math:`\mathsf{M}^{-1} = \mathsf{L}^{-T} \mathsf{L}^{-1}`.
  Further details can be found in:

  - not yet available through *hypre* interface,
  - `PETSc documentation <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCICC.html>`__,
  - `Trilinos documentation <https://docs.trilinos.org/dev/packages/ifpack/doc/html/classIfpack__IC.html>`__.

* **AMG**: algebraic multigrid (can be classical or aggregation-based according to the specific package).
  Further details can be found in:

  - `HYPRE documentation <https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html>`__,
  - `PETSc documentation <https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCGAMG.html>`__,
  - `Trilinos documentation <https://docs.trilinos.org/dev/packages/ml/doc/html/index.html>`__.

* **MGR**: multigrid reduction. Available through *hypre* interface only. Specific documentation coming soon.
  Further details can be found in `MGR documentation <https://hypre.readthedocs.io/en/latest/solvers-mgr.html>`__.

* **Block**: custom preconditioner designed for a 2 x 2 block matrix.

********************
Block preconditioner
********************

This framework allows the user to design a block preconditioner for a 2 x 2 block matrix. The key component is the Schur complement
:math:`\mathsf{S} = \mathsf{A}^{11} - \mathsf{A}^{10} \mathsf{\widetilde{A}}_{00}^{-1} \mathsf{A}_{01}` computation, that requires
an approximation of the leading block. Currently, available options for :math:`\mathsf{\widetilde{A}}_{00}^{-1}` are:

* diagonal with diagonal values (essentially, a Jacobi preconditioner);
* diagonal with row sums as values (e.g., used for CPR-like preconditioners).

Once the Schur complement is computed, to properly define the block preconditioner we need:

* the preconditioner for :math:`\mathsf{A}_{00}` (any of the above listed single-matrix preconditioner);
* the preconditioner for :math:`\mathsf{S}` (any of the above listed single-matrix preconditioner);
* the application strategy. This can be:

  * diagonal: none of the coupling terms is used;
  * upper triangular: only the upper triangular coupling term is used;
  * lower-upper triangular: both coupling terms are used.

Moreover, a block scaling is available. Feasible options are:

* none: keep the original scaling;
* Frobenius norm: equilibrate Frobenius norm of the diagonal blocks;
* user provided.
