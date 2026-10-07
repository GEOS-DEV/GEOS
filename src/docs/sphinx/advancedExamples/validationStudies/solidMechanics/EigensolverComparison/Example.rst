.. _AdvancedExampleEigensolverComparison:


#####################################################################
Free Block and Free-Free Beam: Arnoldi versus LOBPCG Eigensolvers
#####################################################################

------------------------------------------------------------------
Problem description
------------------------------------------------------------------

The ``Modal`` time integration option of the :ref:`SolidMechanicsLagrangianFEM` solver offers two eigensolvers: ``arnoldi`` and ``lobpcg``.
This example compares them on two free bodies.
It shows the effect of the deflation of the rigid-body modes.
The theory of both solvers is in :ref:`SolidMechanicsLagrangianFEM`.

The first problem is a free steel block of 2 m by 1 m by 1 m.
It has 32 by 16 by 16 trilinear hexahedra, that is 28 611 unknowns.
The solver computes 16 modes: the six rigid-body modes and ten elastic modes.
The mesh is isotropic and the multigrid preconditioner is effective.

The second problem is the slender free-free beam of the example :ref:`AdvancedExampleFreeFreeBeamModes`, with 7 575 unknowns.
The multigrid preconditioner is less effective for this problem, because of the aspect ratio of the beam.
The solver computes 10 modes, and 20 modes in one run.

**Input files**

The block uses two GEOS xml files that differ by the eigensolver:

.. code-block:: console

  inputFiles/solidMechanics/modalFreeBlock_base.xml
  inputFiles/solidMechanics/modalFreeBlock_arnoldi.xml
  inputFiles/solidMechanics/modalFreeBlock_lobpcg.xml

The beam uses the files of the example :ref:`AdvancedExampleFreeFreeBeamModes`.
The other combinations are obtained by changing the attributes ``modalSolverType`` and ``modalDeflateRigidBodyModes``.
The script :download:`run_comparison.sh <run_comparison.sh>` makes all the runs.
The script :download:`summarize_comparison.py <summarize_comparison.py>` writes the table and the convergence histories from the logs.
The script :download:`EigensolverComparison.py <EigensolverComparison.py>` makes the figure.
The results are in :download:`EigensolverComparison.csv <EigensolverComparison.csv>` and :download:`ConvergenceHistory.csv <ConvergenceHistory.csv>`.

------------------------------------------------------------------
Solver
------------------------------------------------------------------

The ``lobpcg`` input file of the block selects the eigensolver and the deflation of the rigid-body modes:

.. literalinclude:: ../../../../../../../inputFiles/solidMechanics/modalFreeBlock_lobpcg.xml
    :language: xml
    :start-after: <!-- SPHINX_MODAL_BLOCK_SOLVER -->
    :end-before: <!-- SPHINX_MODAL_BLOCK_SOLVER_END -->

The ``LinearSolverParameters`` block gives the preconditioner.
``arnoldi`` uses the conjugate gradient solver with this preconditioner to solve to ``krylovTol``.
``lobpcg`` only applies one cycle of the preconditioner.
The shift is -500 Hz, which is below the first mode.

------------------------------------------------------------------
Results
------------------------------------------------------------------

All the runs give the same frequencies to the printed digits.
For example, the first elastic mode of the block is at 715.9340 Hz, and the one of the beam at 10.8928 Hz.
The table gives the cost of each run.
The operator applications are the linear solves for ``arnoldi`` and the preconditioner applications for ``lobpcg``.
The times are for one core of a workstation and only give the order of magnitude.

.. csv-table:: Cost of the eigensolvers
   :file: EigensolverComparison.csv
   :header-rows: 1

The figure shows the convergence history of ``lobpcg`` and the time of the combinations.
A hatched bar is a run that did not converge.

.. plot:: docs/sphinx/advancedExamples/validationStudies/solidMechanics/EigensolverComparison/EigensolverComparison.py

**Discussion**

- ``arnoldi`` converges for every case.
  Its cost is the number of solves times the cost of a solve.
  The default completeness check keeps the converged modes and expands the Krylov space of a random vector, until the best pair that is not kept has converged.
  It adds 38 solves for the block and 9 for the beam with 10 modes, and it is a large part of the cost of the runs of the block.
  With deflation of the rigid modes, the iteration has no repeated eigenvalue to find.
  The runs with deflation keep the check and are cheaper than the default runs by 1 % for the block and by 29 % for the beam.
  If the multiplicity of the eigenvalues is known, ``modalCompletenessCheck="0"`` removes this cost.
- ``lobpcg`` without deflation needs 1440 preconditioner applications for the block, against 810 with deflation.
  It does not converge for the beam within 400 iterations.
  The rigid-body modes have a rounding noise in :math:`\mathbf{K}\mathbf{x}`.
  The noise keeps them active and disturbs the other modes.
- ``lobpcg`` with deflation converges in 80 iterations for the block and 134 for the beam.
  It is 3.0 times faster than ``arnoldi`` with deflation for the block, and 8.9 times faster for the beam.
  For the 20 modes of the beam, it needs 2.9 s, against 13.8 s for ``arnoldi`` with the default settings.
- The advantage of ``lobpcg`` depends on the quality of the preconditioner.
  The block is well preconditioned, and a solve to a tight tolerance costs about 35 iterations of the conjugate gradient method.
  The beam needs about 270 iterations for each solve, so the accurate solves of ``arnoldi`` are expensive and the single cycles of ``lobpcg`` are cheap.

**Recommendation**

Use ``arnoldi`` as the default, because it is robust and finds the modes closest to any shift.
For a free body, use ``modalDeflateRigidBodyModes="1"`` with both solvers.
Use ``lobpcg`` for a free body when only the lowest modes are needed and the preconditioner is effective.
Increase ``modalSubspaceSize`` above ``modalNumModes`` if a group of repeated eigenvalues is cut at the last mode.

------------------------------------------------------------------
To go further
------------------------------------------------------------------

**References**

- [Knyazev2001]_ introduces the LOBPCG method.
- [Duersch2018]_ and [HetmaniukLehoucq2006]_ describe robust implementations.
- [StathopoulosWu2002]_ gives the orthonormalization used.
- [Stewart2001]_, [WuSimon2000]_ and [LehoucqSorensenYang1998]_ describe the Krylov-Schur and the Arnoldi methods.

**Feedback on this example**

For any feedback on this example, please submit a `GitHub issue on the project's GitHub page <https://github.com/GEOS-DEV/GEOS/issues>`_.
