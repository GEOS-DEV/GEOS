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

where :math:`\mathsf{A}` is the square sparse matrix, :math:`\mathsf{x}` the solution vector, and :math:`\mathsf{b}` the right-hand side.
For example, in a classical linear elastostatics problem :math:`\mathsf{A}` is the stiffness matrix, and :math:`\mathsf{x}` and :math:`\mathsf{b}` are the displacement and nodal force vectors, respectively.


This solution stage represents the most computationally expensive portion of a typical simulation.
Solution algorithms generally belong to two families of methods: direct methods and iterative methods.
In GEOS both options are made available wrapping around well-established open-source linear algebra libraries, namely
`HYPRE <https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__,
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

The default option in GEOS relies on `SuperLU <http://crd-legacy.lbl.gov/~xiaoye/SuperLU/>`__, a general purpose library for the direct solution of large, sparse, nonsymmetric systems of linear equations, that is called taking advantage of the interface provided in `HYPRE <https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__.

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

.. include:: /docs/sphinx/datastructure/LinearSolverParameters.rst

************************
HYPRE through hypredrive
************************

When GEOS is built with both ``HYPRE`` and ``HYPREDRV`` support, the standard HYPRE-based
linear solver path can delegate iterative solves to `hypredrive <https://github.com/LLNL/hypredrive>`__.
This provides a YAML-driven way to configure HYPRE solvers and preconditioners while keeping
GEOS backward compatible with the existing ``LinearSolverParameters`` input block.

At runtime, GEOS still owns the assembled linear system objects and passes the following data
to hypredrive in library mode:

* the HYPRE matrix,
* the right-hand side and solution vectors,
* the degree-of-freedom map when it is available.

The hypredrive input only controls solver and preconditioner configuration. Matrix assembly,
vector ownership, and all higher-level GEOS solver logic remain unchanged.

Build requirements
==================

This feature is optional. It is available only when GEOS is configured with:

* ``ENABLE_HYPRE=ON``,
* ``ENABLE_HYPREDRV=ON``,
* ``HYPREDRV_DIR`` pointing to a hypredrive installation or package configuration directory.

If GEOS is built without hypredrive support, setting ``hypredriveInputFile`` in the XML input
is an error.

When enabled, hypredrive owns HYPRE runtime initialization and finalization. No additional user
input is required for this.

Two configuration modes
=======================

GEOS can use hypredrive in two different ways.

1. **Authoritative YAML file**: the user supplies a hypredrive YAML file through
   ``LinearSolverParameters/hypredriveInputFile``.
2. **Generated fallback YAML**: if no file is supplied, GEOS translates the supported
   ``LinearSolverParameters`` options into an in-memory hypredrive YAML document.

In both cases, the switch is transparent to the calling physics solver.

Authoritative YAML file
=======================

The most flexible workflow is to provide an explicit hypredrive YAML file:

.. code-block:: xml

  <SinglePhasePoromechanics
    name="PoroelasticitySolver"
    solidSolverName="LinearElasticitySolver"
    flowSolverName="SinglePhaseFlowSolver"
    targetRegions="{ Domain }">
    <LinearSolverParameters
      solverType="gmres"
      preconditionerType="mgr"
      logLevel="2"
      hypredriveInputFile="inputFiles/linearSolvers/flow-mgr.yml"/>
  </SinglePhasePoromechanics>

When ``hypredriveInputFile`` is provided, the file is **authoritative**:

* GEOS passes the file path to hypredrive as-is.
* GEOS does **not** merge, patch, or override the YAML ``solver`` or ``preconditioner`` blocks.
* GEOS only injects the runtime matrix, vectors, and dof map.
* Any ``preconditioner.reuse`` settings in that YAML file are preserved and become the source of
  truth for reuse policy.

This mode is recommended when:

* the desired hypredrive configuration uses options that are not exposed through GEOS XML,
* the user wants exact control of solver and preconditioner YAML,
* the user wants to reproduce a standalone hypredrive configuration inside GEOS.

For readability it is still a good idea to keep the XML ``solverType`` and
``preconditionerType`` consistent with the supplied YAML file, but the YAML file is the source
of truth for hypredrive solver options.

When hypredrive is active, GEOS keeps one ``HYPREDRV_t`` handle alive per
``HypredriveSolver`` and reuses that handle across compatible setup/solve cycles. GEOS only
recreates the handle when the authoritative/generated YAML changes, when the system-setup
timestamp changes, or when the matrix/dof-layout signature changes. This persistent handle is
what allows hypredrive's own preconditioner-reuse logic to work across repeated GEOS linear
solves.

Generated YAML from GEOS input
==============================

If ``hypredriveInputFile`` is left empty, GEOS attempts to generate an equivalent hypredrive
input from the existing ``LinearSolverParameters`` values. This preserves the familiar GEOS XML
workflow while letting hypredrive manage the actual HYPRE solver objects.

In this generated path, GEOS currently keeps ``preconditioner.reuse`` disabled. Reuse is
available today through the authoritative-YAML mode, where the user can write the reuse block
directly in hypredrive YAML.

The generated path currently supports the following iterative solver types:

* ``cg`` (translated to hypredrive ``pcg``),
* ``gmres``,
* ``fgmres``,
* ``bicgstab``.

The generated path currently supports the following preconditioners:

* ``none``,
* ``amg``,
* ``mgr``,
* ``iluk``,
* ``ilut``.

For generated AMG YAML, GEOS translates GEOS smoother and coarse-grid options to the canonical
hypredrive names. Important examples are:

* ``fgs`` -> ``forward-hl1gs``,
* ``bgs`` -> ``backward-hl1gs``,
* ``l1jacobi`` -> ``l1-jacobi``,
* ``sgs`` -> ``hsgs``,
* ``l1sgs`` -> ``l1-hsgs``,
* ``direct`` coarse solve -> ``ge``.

For example, a generated AMG configuration passed to hypredrive has the form:

.. code-block:: yaml

  solver:
    gmres:
      max_iter: 300
      relative_tol: 1.0e-8
      krylov_dim: 200
  preconditioner:
    amg:
      tolerance: 0.0
      max_iter: 1
      relaxation:
        down_type: forward-hl1gs
        up_type: forward-hl1gs
        coarse_type: ge

If GEOS cannot represent the requested configuration as valid hypredrive YAML, it automatically
falls back to the legacy GEOS HYPRE setup path. This keeps existing inputs working even when
they use options that are not yet translated to hypredrive.

Common cases that currently stay on the legacy HYPRE path are:

* unsupported generated solver or preconditioner combinations,
* ``krylovAdaptiveTol`` enabled,
* AMG options that are not translated by the generated YAML builder.

When exact control is required, prefer the authoritative YAML-file mode.

Generated MGR YAML and symbolic dof labels
==========================================

The generated hypredrive MGR path mirrors the MGR strategies already implemented in GEOS.
In addition to the ``solver`` and ``preconditioner`` sections, GEOS emits a
``linear_system.dof_labels`` block so that MGR levels can refer to symbolic names instead of raw
integer component ids.

The labels are derived from the GEOS degree-of-freedom manager:

* field names are taken in registration order,
* names are sanitized to YAML-friendly tokens,
* each component is suffixed with its component index.

For a displacement-pressure system, the generated YAML may look like:

.. code-block:: yaml

  linear_system:
    dof_labels:
      u_0: 0
      u_1: 1
      u_2: 2
      p: 3
  preconditioner:
    mgr:
      level:
        0:
          f_dofs: [u_0, u_1, u_2]

This makes the generated MGR configuration much easier to inspect and compare with standalone
hypredrive YAML files.

Logging, debugging, and fallback control
========================================

The linear-solver ``logLevel`` can be used to inspect what GEOS is passing to hypredrive:

* if ``logLevel < 1``, no additional hypredrive input is printed,
* if ``logLevel >= 1`` and ``hypredriveInputFile`` is set, GEOS logs a delimited YAML block with a comment identifying the authoritative file and, when readable from GEOS, the authoritative YAML contents,
* if ``logLevel >= 1`` and GEOS generates the YAML itself, GEOS logs the full generated YAML in the same delimited format,
* GEOS emits this hypredrive input dump once per hypredrive-using solver during startup, after solver initialization has finalized the effective configuration and before the ``Import fields`` log section,
* for those solvers, GEOS suppresses the usual GEOS linear-solver parameter table and logs only the hypredrive YAML block,
* the hypredrive YAML does not repeat on later setups,
* when hypredrive statistics printing is enabled, the destroy-time ``STATISTICS SUMMARY`` banner uses ``general.name`` when present in authoritative YAML, and GEOS labels generated-fallback objects with the owning solver name automatically,
* internally, GEOS drives hypredrive through handle-scoped ``system``, ``timestep``, and
  ``newton`` annotations inside the hypredrive adapter rather than from the generic physics
  solver layer.

This is especially useful when validating a new solver recipe or comparing the generated GEOS
configuration with a standalone hypredrive input file.

For regression studies and side-by-side comparisons, hypredrive can be disabled at runtime by
setting the environment variable:

.. code-block:: bash

  export GEOS_HYPREDRV_FORCE_LEGACY=1

When this variable is set, GEOS bypasses hypredrive and uses the legacy HYPRE implementation even
if GEOS was built with hypredrive support.

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

* **MGR**: multigrid reduction. Available through *hypre* interface only.
  Further details can be found in `MGR documentation <https://hypre.readthedocs.io/en/latest/solvers-mgr.html>`__, also see section below.

* **Block**: custom preconditioner designed for a 2 x 2 block matrix.

************************
HYPRE MGR Preconditioner
************************

MGR stands for multigrid reduction, a multigrid method that uses the interpolation, restriction operators, and the Galerkin triple product, to reduce a linear system to a smaller one, similar to a Schur complement approach. As such, it is designed to target block linear systems resulting from discretizations of multiphysics problems. GEOS uses MGR through an implementation in `HYPRE <https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__. More information regarding MGR can be found `here <https://hypre.readthedocs.io/en/latest/solvers-mgr.html>`__. Currently, MGR strategies are implemented for hydraulic fracturing, singlephase and multiphase poromechanics, singlephase poromechanics with fractures, compositional flow with and without wells. More multiphysics solvers with MGR will be enabled in the future.

To use MGR for a specific block system, several components need to be specified.

(1) **The number of reduction levels and the coarse points (corresponding to fields) for each level**. For example, for single-phase hydraulic fracturing, there are two fields, i.e. displacement and fluid pressure, a two-level MGR strategy can be used with the fluid pressure being the coarse degrees of freedom.
(#) **Interpolation/restriction operators and the coarse-grid computation strategy**. A simple but effective strategy is to use Jacobi diagonalization for interpolation and injection for restriction. For most cases, a Galerkin coarse grid strategy can be used, but for special cases such as poroelastic, a non-Galerkin approach is preferable.
(#) **Global smoother**. Depending on the problem, a global relaxation step could be beneficial. Some options include ILU(k), (block) Jacobi, (block) Gauss-Seidel.
(#) **Solvers for F-relaxation and coarse-grid correction**. These solvers should be chosen carefully for MGR to be effective. The choice of these solvers should correspond to the properties of the blocks specified by the C- and F-points. For example, if the :math:`\mathsf{A}_{FF}` block is hyperbolic, a Jacobi smoother is sufficient while for an elliptic operator an AMG V-cycle might be required. For the single-phase hydraulic fracturing case, an AMG V-cycle is needed for both F-relaxation and coarse-grid correction.


Note that these are only general guidelines for designing a working MGR recipe. For complicated multiphysics problems, experimentation with different numbers of levels, choices of C- and F-points, and smoothers/solvers, etc., is typically needed to find the best strategy. Currently, these options are only available to developers. We are working on exposing these functionalities to the users in future releases.

********************
Block preconditioner
********************

This framework allows the user to design a block preconditioner for a 2 x 2 block matrix. The key component is the Schur complement
:math:`\mathsf{S} = \mathsf{A}_{11} - \mathsf{A}_{10} \mathsf{\widetilde{A}}_{00}^{-1} \mathsf{A}_{01}` computation, that requires
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

********************
Adaptive tolerance
********************

This feature is available for iterative solvers and can be enabled using `krylovAdaptiveTol` flag in `LinearSolverParameters`. It follows the Eisenstat-Walker inexact Newton approach described in [Eisenstat and Walker 1996]. The key idea is to relax the linear solver tolerance at the beginning of the nonlinear iterations loop and tighten it when getting closer to the final solution. The initial tolerance is defined by `krylovWeakestTol` and starting from second nonlinear iteration the tolerance is chosen using the following steps:

- compute the current to previous nonlinear norm ratio: :math:`\mathsf{nr} = \mathsf{min}( \mathsf{norm}^{curr} / \mathsf{norm}^{prev}, 1.0 )`
- estimate the new linear solver tolerance: :math:`\mathsf{tol}_{new} = \mathsf{\gamma} \cdot \mathsf{nr}^{ax}`
- compute a safeguard to avoid too sharp tolerance reduction: :math:`\mathsf{tol}_{alt} = \mathsf{tol}_{old}^{2}` (the bound is the quadratic reduction with respect to the previous tolerance value)
- apply safeguards and compute the final tolerance: :math:`\mathsf{tol} = \mathsf{max}( \mathsf{tol}_{new}, \mathsf{tol}_{alt} )`, :math:`\mathsf{tol} = \mathsf{min}( \mathsf{tol}_{max}, \mathsf{max}( \mathsf{tol}_{min}, \mathsf{tol} ) )`

Here :math:`\mathsf{\gamma}` is the forcing term, :math:`ax` is the adaptivity exponent, :math:`\mathsf{tol}_{min}` and :math:`\mathsf{tol}_{max}` are prescribed tolerance bounds (defined by `krylovStrongestTol` and `krylovWeakestTol`, respectively).
