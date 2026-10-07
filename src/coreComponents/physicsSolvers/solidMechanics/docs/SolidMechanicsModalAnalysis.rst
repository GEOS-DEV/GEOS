.. _SolidMechanicsModalAnalysis:

#####################################
Solid Mechanics Modal Analysis Solver
#####################################

Introduction
============
The `SolidMechanicsModalAnalysis` solver computes the vibration modes of a structure: the free vibrations about the current state, found as the generalized eigenvalue problem :math:`\mathbf{K} \boldsymbol{\phi} = \lambda \mathbf{M} \boldsymbol{\phi}`.
It derives from the :ref:`SolidMechanicsLagrangianFEM` solver and uses its finite element discretization, material models and boundary conditions.
Unlike this solver, it has no time integration: one execution computes the modes.

Theory
======


At each application, the solver does not advance in time. It computes the free vibration modes of the structure about the current state.
The tangent stiffness matrix :math:`\mathbf{K}` and the mass matrix :math:`\mathbf{M}` define the generalized eigenvalue problem

.. math::
   \mathbf{K} \boldsymbol{\phi}_k = \lambda_k \mathbf{M} \boldsymbol{\phi}_k,
   \qquad \lambda_k = \omega_k^2, \qquad f_k = \frac{\omega_k}{2\pi},

where :math:`f_k` is the frequency of mode :math:`k`.
Both matrices are symmetric and positive semi-definite.
``modalMassType`` selects the mass matrix.
The default, ``lumped``, is diagonal: the nodal mass is the row sum of the consistent mass matrix, which is the quantity used by the explicit dynamics option.
The ``consistent`` mass is the exact matrix :math:`M_{ab} = \int_\Omega \rho N_a N_b \, d\Omega`.
It is available for first-order tetrahedra, where it is :math:`\rho V (1 + \delta_{ab}) / 20` for each pair of nodes :math:`a, b` and each direction.
The modes are normalized so that :math:`\boldsymbol{\phi}_k^T \mathbf{M} \boldsymbol{\phi}_k = 1`.
Then :math:`\boldsymbol{\phi}_k^T \mathbf{K} \boldsymbol{\phi}_k = \lambda_k`, and two modes with different eigenvalues are orthogonal for :math:`\mathbf{M}` and for :math:`\mathbf{K}` [Bathe1996]_ [Parlett1998]_.

Free body and constraints
-------------------------

Displacement boundary conditions are treated as homogeneous.
The constrained degrees of freedom are removed symmetrically from :math:`\mathbf{K}` and :math:`\mathbf{M}`, and the mode shapes are zero there.
If no displacement boundary condition is imposed, the structure is a free body.
Then :math:`\mathbf{K}` is singular and the six rigid-body modes are the first six eigenpairs, with :math:`\lambda = 0`.
They are three translations and three rotations.
The participation factors of the six rigid modes satisfy :math:`\sum_{k=1}^{6} \Gamma_{k,d}^2 = m_{tot}` for each direction :math:`d`, where :math:`m_{tot}` is the total mass.

With ``modalVerifyFreeBody="1"``, the solver checks a free body strictly.
It requires no displacement boundary condition and at least seven modes.
Before the eigensolve, it checks that the six analytical rigid vectors are in the null space of :math:`\mathbf{K}`, with a stiffness residual relative to :math:`\|\mathbf{K}\|_\infty`, and that they are independent in the :math:`\mathbf{M}` inner product.
It solves the eigenproblem with a tolerance 1000 times smaller than ``modalTolerance``.
After it, it checks that exactly six eigenvalues are zero, relative to the largest one, and that the residuals of the eigenpairs of the original pencil, the :math:`\mathbf{M}`-orthonormality and the :math:`\mathbf{K}`-orthogonality of the modes are within ``modalTolerance``.
If a check fails, the solver stops with an error.

Spectral shift
--------------

The modes closest to a spectral shift :math:`\sigma` are computed.
The shift is given as a signed frequency by ``modalShiftFrequency``, and :math:`\sigma = \mathrm{sign}(f)\,(2 \pi f)^2`.
A negative frequency gives a negative shift :math:`\sigma = -\alpha`.
The shifted matrix :math:`\mathbf{K} - \sigma \mathbf{M} = \mathbf{K} + \alpha \mathbf{M}` is then symmetric positive definite, even for a free body.
Its linear solver is the one given in the ``LinearSolverParameters`` block of the solver.
The multigrid preconditioner of an elastic problem uses the rigid-body modes as near-null space.

Two eigensolvers are available, with the ``modalSolverType`` attribute.
Both are accessed through a common interface, so a new method can be added with one class and one entry in the factory.

Arnoldi eigensolver
-------------------

The ``arnoldi`` eigensolver is a block Krylov-Schur method with spectral transformation [Stewart2001]_ [WuSimon2000]_.
It is the shift-and-invert approach of ARPACK [LehoucqSorensenYang1998]_ [EricssonRuhe1980]_.
It applies the operator

.. math::
   \mathbf{T} = (\mathbf{K} - \sigma \mathbf{M})^{-1} \mathbf{M},

which is self-adjoint for the :math:`\mathbf{M}` inner product.
Its eigenvalues are :math:`\theta = 1 / (\lambda - \sigma)`, so the modes closest to the shift have the largest :math:`|\theta|`.
The solver builds an :math:`\mathbf{M}`-orthonormal basis :math:`V_m` of the Krylov space of :math:`\mathbf{T}`, with the relation

.. math::
   \mathbf{T} V_m = V_m H_m + V_{next} B,

and extracts the Ritz pairs :math:`(\theta_i, V_m y_i)` from the small symmetric matrix :math:`H_m`.
Each new basis vector costs one linear solve with :math:`\mathbf{K} - \sigma \mathbf{M}`.
The basis is restarted with the Krylov-Schur technique: the solver keeps the best Ritz vectors and the residual block.
The default basis has the larger of twice the number of modes and 20 vectors, set by ``modalSubspaceSize``.
For a single vector this is the Lanczos method with thick restart.
The implementation keeps the full matrix :math:`H_m` and :math:`\mathbf{M}`-orthogonalizes every vector twice, so that it does not depend on a tridiagonal structure.

The residual of a Ritz pair satisfies

.. math::
   \| \mathbf{T} x - \theta x \|_{\mathbf{M}} = |\theta| \, \| (\mathbf{K} - \sigma \mathbf{M})^{-1} r \|_{\mathbf{M}},
   \qquad r = \mathbf{K} x - \lambda \mathbf{M} x,

for an :math:`\mathbf{M}`-normalized :math:`x`.
The convergence test :math:`\| \mathbf{T} x - \theta x \|_{\mathbf{M}} \leq \mathrm{tol}\,|\theta|` is therefore a test on :math:`\| (\mathbf{K} - \sigma \mathbf{M})^{-1} r \|_{\mathbf{M}}`.
Because every solve must be accurate, set ``krylovTol`` at least two orders of magnitude below ``modalTolerance``.

A single Krylov vector cannot find a repeated eigenvalue.
The six rigid-body modes of a free body are such a case.
Two remedies are available.
The default ``modalCompletenessCheck`` keeps the converged modes and builds the Krylov space of a random vector that is :math:`\mathbf{M}`-orthogonal to them.
That space has a direction in every eigenspace that is left, so a missed copy appears among the wanted Ritz pairs.
The check ends only when the best Ritz pair that is not kept has converged, because a Ritz value that has not converged can hide a missed copy when the Ritz values are clustered.
It repeats while it finds new modes.
The alternative is ``modalBlockSize`` at least equal to the multiplicity of the eigenvalues.
Neither is a proof that no mode is missed.
The check is a strong heuristic, because the largest eigenvalues of the complement appear first in the Krylov space.

LOBPCG eigensolver
------------------

The ``lobpcg`` eigensolver is the locally optimal block preconditioned conjugate gradient method [Knyazev2001]_.
It computes the lowest eigenvalues of the pencil :math:`(\mathbf{K}, \mathbf{M})`.
It never solves a linear system accurately: it applies a preconditioner :math:`\mathbf{T} \approx (\mathbf{K} - \sigma \mathbf{M})^{-1}`, which is one cycle of the preconditioner of the ``LinearSolverParameters`` block (for example algebraic multigrid).
For a direct linear solver, the exact solve is the preconditioner.
The shift only enters through this preconditioner.
It must be at or below the first eigenvalue, so that :math:`\mathbf{K} - \sigma \mathbf{M}` is positive definite.

The block :math:`X = [x_1, \dots, x_n]` of :math:`\mathbf{M}`-orthonormal iterates is improved at each iteration in three steps.

1. The residuals :math:`r_i = \mathbf{K} x_i - \lambda_i \mathbf{M} x_i` are computed, and preconditioned: :math:`w_i = \mathbf{T} r_i`.
2. The search space :math:`S = [X, W, P]` is formed with the iterates, the preconditioned residuals and the previous search directions :math:`P`.
3. The Rayleigh-Ritz problem is solved on :math:`S`.

The Rayleigh-Ritz problem is made well-posed by the SVQB algorithm [StathopoulosWu2002]_.
The Gram matrix :math:`G = S^T \mathbf{M} S` is scaled to a unit diagonal, :math:`D G D = U \Sigma U^T`, and the directions with a negligible eigenvalue are dropped.
The matrix :math:`C = D U \Sigma^{-1/2}` satisfies :math:`C^T G C = I`.
The solver then diagonalizes the small matrix

.. math::
   C^T \left( S^T \mathbf{K} S \right) C \, y = \theta \, y,

and computes the new iterates :math:`X^{+} = S C y` for the lowest :math:`n` Ritz values.
The new search directions :math:`P^{+}` are the part of :math:`X^{+}` that comes from :math:`W` and :math:`P`.
The use of :math:`P` gives the method the three-term recurrence of the conjugate gradient method.
The orthonormalization of the basis follows the robust implementations of [HetmaniukLehoucq2006]_ and [Duersch2018]_.

The implementation has these properties:

- The error estimate is :math:`\| \mathbf{T} r_i \|_{\mathbf{M}}`, for :math:`\mathbf{M}`-normalized iterates.
  It estimates the quantity that the Arnoldi test uses.
  The plain residual :math:`\| r_i \|` cannot be used for the eigenvalues close to zero, because it has a floor: the rounding error of :math:`\mathbf{K} x`.
- The converged iterates stay in the Rayleigh-Ritz space but do not make new directions (soft locking).
- The images :math:`\mathbf{K} X` and :math:`\mathbf{M} X` are recomputed at each iteration, because the update by coefficients loses accuracy when the directions become nearly dependent.
- If the search space has fewer independent directions than modes, random vectors are added.
- ``modalSubspaceSize`` larger than ``modalNumModes`` adds guard vectors.
  They keep the method from cutting a group of repeated eigenvalues at the last requested mode.

In each iteration, the preconditioner and :math:`\mathbf{M}` are applied to the residual of every iterate, converged or not, because the convergence test uses the preconditioned residual of all of them.
:math:`\mathbf{K}` is applied to the residual of every active vector, and :math:`\mathbf{K}` and :math:`\mathbf{M}` to the new iterates.
The dense work is :math:`O(n^2)` vector operations.
The memory is about fifteen blocks of :math:`n` vectors.

Deflation of the rigid-body modes
---------------------------------

The rigid-body modes of a free structure are known.
With ``modalDeflateRigidBodyModes="1"``, the solver computes the three translations and the three rotations from the node coordinates.
It :math:`\mathbf{M}`-orthonormalizes them into a matrix :math:`Y`, and projects every new vector with :math:`z \leftarrow z - Y Y^T \mathbf{M} z`.
Both eigensolvers then compute the other modes in the :math:`\mathbf{M}`-orthogonal complement of the rigid modes.
The rigid modes are the first modes of the result and count in ``modalNumModes``.
Their eigenvalue is the Rayleigh quotient of the analytical vector, which is zero up to rounding.
This option needs a structure without displacement boundary conditions.

Deflation has two effects.
It removes the repeated eigenvalue, so the Arnoldi iteration does not need the completeness check.
It also removes the rounding noise of :math:`\mathbf{K} x` for the rigid modes.
This noise keeps the rigid modes in the active set of LOBPCG and prevents the convergence of the other modes.
Deflation is therefore advised with ``lobpcg``.

Choice of the eigensolver
-------------------------

The cost of ``arnoldi`` is the number of operator applications, about two to six times the number of modes, times the cost of an accurate linear solve.
The cost of ``lobpcg`` is the number of iterations, typically 80 to 200, times the number of modes, times the cost of one preconditioner application.
It is faster when one preconditioner cycle is much cheaper than a solve to a tight tolerance, which is the usual case.
The examples :ref:`AdvancedExampleEigensolverComparison` and :ref:`AdvancedExampleFreeFreeBeamModes` give measured costs.
In summary:

- ``arnoldi`` is the robust default.
  It finds the modes closest to any shift, including a shift inside the spectrum with a direct solver, and its behavior only depends on the quality of the linear solver.
- ``lobpcg`` with ``modalDeflateRigidBodyModes="1"`` is about three to nine times faster than ``arnoldi`` with deflation on the cases of the examples.
  Without deflation it can stagnate.
  It finds the lowest modes only.
  Its convergence depends on the quality of the preconditioner.

The convergence test of both solvers is on :math:`\| (\mathbf{K} - \sigma \mathbf{M})^{-1} \mathbf{r} \|_{\mathbf{M}} \leq` ``modalTolerance``.
The log reports the relative residual :math:`\|\mathbf{r}\|_2 / (|\lambda - \sigma| \, \|\mathbf{M}\boldsymbol{\phi}\|_2)`.
For the rigid-body modes, this residual is limited by the rounding error of :math:`\mathbf{K}\boldsymbol{\phi}` and can be larger than the tolerance.

Results and limits
------------------

The results are:

- A table in the log with the frequency, the eigenvalue, the residual and the participation factors :math:`\Gamma_{k,d} = \boldsymbol{\phi}_k^T \mathbf{M} \mathbf{e}_d` of each mode in the three directions :math:`d`.
- The arrays ``modalEigenvalues``, ``modalFrequencies``, ``modalResiduals`` and ``modalParticipationFactors`` of the solver, which are saved in the restart files.
  A frequency has the sign of its eigenvalue, so a rigid mode with a small negative eigenvalue has a small negative frequency.
- The nodal fields ``modeShape1``, ``modeShape2``, ..., one for each mode, which a VTK output writes.
  They are not saved in the restart files: after a restart, the frequencies and factors are available, and the mode shapes are zero until the modal analysis runs again.

The following limits apply.
The stiffness is the tangent stiffness at the current state, without the geometric (pre-stress) stiffness.
The consistent mass is only available for first-order tetrahedra.
Contact, damping and body-force or traction loads are not used.
The solver is standalone: a coupled solver cannot drive it.

The modal analysis runs on the CPU, CUDA and HIP backends, and it does not use unified memory.
The vectors, the matrices and the preconditioner stay in the memory space of the backend during the eigensolve.
The dot products of the orthogonalization and of the Rayleigh-Ritz steps are computed in batches by device kernels, and only the small matrix of results is copied to the host.
Linear combinations of many vectors are done by one fused kernel.
The small dense eigenproblems of the projected matrices are solved on the host.

The examples :ref:`AdvancedExampleFreeFreeBeamModes` and :ref:`AdvancedExampleEigensolverComparison` verify this solver against the Euler-Bernoulli beam theory and compare the two eigensolvers.

.. code-block:: xml

   <SolidMechanicsModalAnalysis name="solid"
                                discretization="FE1"
                                targetRegions="{ Region }"
                                modalNumModes="10"
                                modalShiftFrequency="-100"
                                modalTolerance="1e-8">
     <LinearSolverParameters solverType="cg"
                             preconditionerType="amg"
                             krylovTol="1e-10"/>
   </SolidMechanicsModalAnalysis>

   <Events maxTime="1">
     <PeriodicEvent name="modalAnalysis" forceDt="1" target="/Solvers/solid"/>
   </Events>

Modal analysis references
-------------------------

.. [Bathe1996] K.-J. Bathe. *Finite Element Procedures*. Prentice Hall, 1996.

.. [Parlett1998] B. N. Parlett. *The Symmetric Eigenvalue Problem*. SIAM, Classics in Applied Mathematics 20, 1998.

.. [LehoucqSorensenYang1998] R. B. Lehoucq, D. C. Sorensen & C. Yang. *ARPACK Users' Guide: Solution of Large-Scale Eigenvalue Problems with Implicitly Restarted Arnoldi Methods*. SIAM, 1998.

.. [EricssonRuhe1980] T. Ericsson & A. Ruhe. "The spectral transformation Lanczos method for the numerical solution of large sparse generalized symmetric eigenvalue problems". *Math. Comp.* 35(152), 1251-1268, 1980.

.. [WuSimon2000] K. Wu & H. Simon. "Thick-restart Lanczos method for large symmetric eigenvalue problems". *SIAM J. Matrix Anal. Appl.* 22(2), 602-616, 2000.

.. [Stewart2001] G. W. Stewart. "A Krylov-Schur algorithm for large eigenproblems". *SIAM J. Matrix Anal. Appl.* 23(3), 601-614, 2001.

.. [Knyazev2001] A. V. Knyazev. "Toward the optimal preconditioned eigensolver: locally optimal block preconditioned conjugate gradient method". *SIAM J. Sci. Comput.* 23(2), 517-541, 2001.

.. [StathopoulosWu2002] A. Stathopoulos & K. Wu. "A block orthogonalization procedure with constant synchronization requirements". *SIAM J. Sci. Comput.* 23(6), 2165-2182, 2002.

.. [HetmaniukLehoucq2006] U. Hetmaniuk & R. Lehoucq. "Basis selection in LOBPCG". *J. Comput. Phys.* 218(1), 324-332, 2006.

.. [Duersch2018] J. A. Duersch, M. Shao, C. Yang & M. Gu. "A robust and efficient implementation of LOBPCG". *SIAM J. Sci. Comput.* 40(5), C655-C676, 2018.

.. [Blevins2001] R. D. Blevins. *Formulas for Natural Frequency and Mode Shape*. Krieger, 2001.

.. [Han1999] S. M. Han, H. Benaroya & T. Wei. "Dynamics of transversely vibrating beams using four engineering theories". *J. Sound Vib.* 225(5), 935-988, 1999.

.. [Cowper1966] G. R. Cowper. "The shear coefficient in Timoshenko's beam theory". *J. Appl. Mech.* 33(2), 335-340, 1966.

.. [TimoshenkoGoodier1970] S. P. Timoshenko & J. N. Goodier. *Theory of Elasticity*. 3rd edition, McGraw-Hill, 1970.

Parameters
=========================

The `SolidMechanicsModalAnalysis` is specified by the title of the subblock of the `Solvers` block.
It has the attributes of the :ref:`SolidMechanicsLagrangianFEM` solver, and the following attributes:

.. include:: /docs/sphinx/datastructure/SolidMechanicsModalAnalysis.rst

The following data are allocated and used by the solver:

.. include:: /docs/sphinx/datastructure/SolidMechanicsModalAnalysis_other.rst
