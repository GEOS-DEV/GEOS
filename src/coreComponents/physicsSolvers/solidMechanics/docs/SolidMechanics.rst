.. _SolidMechanicsLagrangianFEM:

#####################################
Solid Mechanics Solver
#####################################

List of Symbols
===================

.. math::
   i,j,k &\equiv \text {indices over spatial dimensions} \notag \\
   a,b,c &\equiv \text {indices over nodes} \notag \\
   l &\equiv \text {indices over volumetric elements} \notag \\
   q,r,s &\equiv \text {indices over faces} \notag \\
   n &\equiv \text {indices over time} \notag \\
   kiter &\equiv \text {iteration count for non-linear solution scheme} \notag \\
   \Omega &\equiv \text {Volume of continuum body} \notag \\
   \Omega_{crack} &\equiv \text {Volume of open crack} \notag \\
   \Gamma &\equiv \text {External surface of } \Omega \notag \\
   \Gamma_t &\equiv \text {External surface where tractions are applied} \notag \\
   \Gamma_u &\equiv \text {External surface where kinematics are specified} \notag \\
   \Gamma_{crack} &\equiv \text {entire surface of crack} \notag \\
   \Gamma_{cohesive} &\equiv \text {surface of crack subject to cohesive tractions} \notag \\
   \eta_0 &\equiv \text {set of all nodes} \notag \\
   \eta_f & \equiv \text {set of all nodes on flow mesh}  \notag \\
   m &\equiv \text{mass} \notag \\
   \kappa_k &\equiv \text{all elements connected to element k} \notag \\
   \phi & \equiv \text { porosity} \notag \\
   p_f & \equiv \text { fluid pressure} \notag \\
   \mathbf{u} & \equiv \text { displacement} \notag \\
   \mathbf{q} & \equiv \text { volumetric flow rate} \notag \\
   \mathbf{T} & \equiv \text { Cauchy stress} \notag \\
   \rho & \equiv \text { density in the current configuration} \notag \\
   \mathbf{x}& \equiv \text { current position} \notag \\
   \mathbf{w}&\equiv \text { aperture, or gap vector} \notag


Introduction
============
The `SolidMechanicsLagrangianFEM` solver applies a Continuous Galerkin finite element method to solve the linear momentum balance equation.
The primary variable is the displacement field which is discretized at the nodes.

Theory
=========================

Governing Equations
--------------------------

The `SolidMechanicsLagrangianFEM` solves the equations of motion as given by

.. math::
   T_{ij,j} + \rho(b_{i}-\ddot{x}_{i}) = 0,

which is a 3-dimensional expression for the well known expression of Newtons Second Law (:math:`F = m a`).
These equations of motion are discretized using the Finite Element Method,
which leads to a discrete set of residual equations:

.. math::
   (R_{solid})_{ai}=\int\limits_{\Gamma_t} \Phi_a t_i   dA  - \int\limits_\Omega \Phi_{a,j} T_{ij}   dV +\int\limits_\Omega \Phi_a \rho(b_{i}-\Phi_b\ddot{x}_{ib})  dV = 0

Quasi-Static Time Integration
-----------------------------
The Quasi-Static time integration option solves the equation of motion after removing the inertial term, which is expressed by

.. math::
   T_{ij,j} + \rho b_{i} = 0,

which is essentially a way to express the equation for static equilibrium (:math:`\Sigma F=0`).
Thus, selection of the Quasi-Static option will yield a solution where the sum of all forces at a given node is equal to zero.
The resulting finite element discretized set of residual equations are expressed as

.. math::
   (R_{solid})_{ai}=\int\limits_{\Gamma_t} \Phi_a t_i   dA  - \int\limits_\Omega \Phi_{a,j} T_{ij}   dV + \int\limits_\Omega \Phi_a \rho b_{i}  dV = 0,

Taking the derivative of these residual equations wrt. the primary variable (displacement) yields

.. math::
    \pderiv{(R_{solid}^e)_{ai}}{u_{bj}} &=
            - \int\limits_{\Omega^e} \Phi_{a,k} \frac{\partial T_{ik}}{\partial u_{bj}}   dV,

And finally, the expression for the residual equation and derivative are used to express a non-linear system of equations

.. math::
   \left. \left(\pderiv{(R_{solid}^e)_{ai}}{u_{bj}} \right)\right|^{n+1}_{kiter}
   \left( \left. \left({u}_{bj} \right) \right|^{n+1}_{{kiter}+1} - \left. \left({u}_{bj} \right) \right|^{n+1}_{kiter} \right)
   = - (R_{solid})_{ai}|^{n+1}_{kiter} ,

which are solved via the solver package.

Implicit Dynamics Time Integration (Newmark Method)
---------------------------------------------------
For implicit dynamic time integration, we use an implementation of the classical Newmark method.
This update method can be posed in terms of a simple SDOF spring/dashpot/mass model.
In the following, :math:`M` represents the mass, :math:`C` represent the damping of the dashpot, :math:`K`
represents the spring stiffness, and :math:`F` represents some external load.

.. math::
   M a^{n+1} + C v^{n+1} + K u^{n+1} &= F_{n+1},  \\

and a series of update equations for the velocity and displacement at a point:

.. math::
   u^{n+1} &= u^n + v^{n+1/2} \Delta t,  \\
   u^{n+1} &= u^n + \left( v^{n} + \inv{2} \left[ (1-2\beta) a^n + 2\beta a^{n+1} \right] \Delta t \right) \Delta t, \\
   v^{n+1} &= v^n + \left[(1-\gamma) a^n + \gamma a^{n+1} \right] \Delta t.

As intermediate quantities we can form an estimate (predictor) for the end of step displacement and midstep velocity by
assuming zero end-of-step acceleration.

.. math::
   \tilde{u}^{n+1} &= u^n + \left( v^{n} + \inv{2}  (1-2\beta) a^n  \Delta t \right) \Delta t = u^n + \hat{\tilde{u}}\\
   \tilde{v}^{n+1} &= v^n + (1-\gamma) a^n  \Delta t =  v^n + \hat{\tilde{v}}

This gives the end of step displacement and velocity in terms of the predictor with a correction for the end
step acceleration.

.. math::
   u^{n+1} &= \tilde{u}^{n+1} + \beta a^{n+1} \Delta t^2 \\
   v^{n+1} &= \tilde{v}^{n+1} + \gamma a^{n+1} \Delta t

The acceleration and velocity may now be expressed in terms of displacement, and ultimately in terms
of the incremental displacement.

.. math::
   a^{n+1} &= \frac{1}{\beta \Delta t^2} \left(u^{n+1} - \tilde{u}^{n+1} \right)  = \frac{1}{\beta \Delta t^2} \left( \hat{u} - \hat{\tilde{u}} \right) \\
   v^{n+1} &= \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(u^{n+1} - \tilde{u}^{n+1} \right) = \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u} - \hat{\tilde{u}} \right)

plugging these into equation of motion for the SDOF system gives:

.. math::
   M \left(\frac{1}{\beta \Delta t^2} \left(\hat{u} - \hat{\tilde{u}} \right)\right) + C \left( \tilde{v}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u} - \hat{\tilde{u}} \right) \right) + K u^{n+1} &= F_{n+1}  \\

Finally, we assume Rayliegh damping for the dashpot.

.. math::
   C = a_{mass} M + a_{stiff} K

Of course we know that we intend to model a system of equations with many DOF.
Thus the representation for the mass, spring and dashpot can be replaced by our finite element discretized
equation of motion.
We may express the system in context of a nonlinear residual problem

.. math::
    (R_{solid}^e)_{ai} &=
        \int\limits_{\Gamma_t^e} \Phi_a t_i   dA  \\
        &- \int\limits_{\Omega^e} \Phi_{a,j} \left(T_{ij}^{n+1}+  a_{stiff} \left(\pderiv{T_{ij}^{n+1}}{\hat{u}_{bk}} \right)_{elastic} \left( \tilde{v}_{bk}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u}_{bk} - \hat{\tilde{u}}_{bk} \right) \right) \right)  dV \notag \\
        &+\int\limits_{\Omega^e} \Phi_a \rho \left(b_{i}- \Phi_b  \left( a_{mass} \left( \tilde{v}_{bi}^{n+1} + \frac{\gamma}{\beta \Delta t} \left(\hat{u}_{bi} - \hat{\tilde{u}}_{bi} \right) \right) + \frac{1}{\beta \Delta t^2}  \left( \hat{u}_{bi} - \hat{\tilde{u}}_{bi} \right) \right) \right)  dV ,\notag \\
    \pderiv{(R_{solid}^e)_{ai}}{\hat{u}_{bj}} &=
        - \int\limits_{\Omega^e} \Phi_{a,k} \left(\pderiv{T_{ik}^{n+1}}{\hat{u}_{bj}}+  a_{stiff} \frac{\gamma}{\beta \Delta t} \left(\pderiv{T_{ik}^{n+1}}{\hat{u}_{bj}} \right)_{elastic} \right)   dV \notag \\
        &- \left( \frac{\gamma a_{mass}}{\beta \Delta t} + \frac{1}{\beta \Delta t^2}  \right) \int\limits_{\Omega^e} \rho \Phi_a \Phi_c    \pderiv{ \hat{u}_{ci} }{\hat{u}_{bj}}dV .

Again, the expression for the residual equation and derivative are used to express a non-linear system of equations

.. math::
   \left. \left(\pderiv{(R_{solid}^e)_{ai}}{u_{bj}} \right)\right|^{n+1}_{kiter}
   \left( \left. \left({u}_{bj} \right) \right|^{n+1}_{{kiter}+1} - \left. \left({u}_{bj} \right) \right|^{n+1}_{kiter} \right)
   = - (R_{solid})_{ai}|^{n+1}_{kiter} ,

which are solved via the solver package. Note that the derivatives involving :math:`u` and :math:`\hat{u}` are interchangable,
as are differences between the non-linear iterations.

Explicit Dynamics Time Integration  (Special Implementation of Newmark Method with \gamma=0.5, \beta=0)
-------------------------------------------------------------------------------------------------------
For the Newmark Method, if \gamma=0.5, \beta=0, and the inertial term contains a diagonalized "mass matrix",
the update equations may be carried out without the solution of a system of equations.
In this case, the update equations simplify to a non-iterative update algorithm.

First the mid-step velocity and end-of-step displacements are calculated through the update equations

.. math::
   \tensor{v}^{n+1/2} &= \tensor{v}^{n} +  \tensor{a}^n \left( \frac{\Delta t}{2} \right), \text{ and} \\
   \tensor{u}^{n+1} &= \tensor{u}^n + \tensor{v}^{n+1/2} \Delta t.

Then the residual equation/s are calculated, and acceleration at the end-of-step is calculated via

.. math::
   \left( \tensor{M} + \frac{\Delta t}{2} \tensor{C} \right) \tensor{a}^{n+1} &=  \tensor{F}_{n+1} - \tensor{C} v^{n+1/2} - \tensor{K} u^{n+1} .

Note that the mass matrix must be diagonal, and damping term may not include the stiffness based damping
coefficient for this method, otherwise the above equation will require a system solve.
Finally, the end-of-step velocities are calculated from the end of step acceleration:

.. math::
   \tensor{v}^{n+1} &= \tensor{v}^{n+1/2} + \tensor{a}^{n+1} \left( \frac{\Delta t}{2} \right).

Note that the velocities may be stored at the midstep, resulting one less kinematic update.
This approach is typically referred to as the "Leapfrog" method.
However, in GEOS we do not offer this option since it can cause some confusion that results from the
storage of state at different points in time.


Modal Analysis (Vibration Modes)
--------------------------------

With ``timeIntegrationOption="Modal"``, the solver does not advance in time.
At each application of the solver, it computes the free vibration modes of the structure about the current state.
The tangent stiffness matrix :math:`\mathbf{K}` and the lumped mass matrix :math:`\mathbf{M}` define the generalized eigenvalue problem

.. math::
   \mathbf{K} \boldsymbol{\phi}_k = \lambda_k \mathbf{M} \boldsymbol{\phi}_k,
   \qquad \lambda_k = \omega_k^2, \qquad f_k = \frac{\omega_k}{2\pi},

where :math:`f_k` is the frequency of mode :math:`k`.
Both matrices are symmetric and positive semi-definite.
The mass matrix is diagonal.
The nodal mass is the row sum of the consistent mass matrix, which is the quantity used by the explicit dynamics option.
The modes are normalized so that :math:`\boldsymbol{\phi}_k^T \mathbf{M} \boldsymbol{\phi}_k = 1`.
Then :math:`\boldsymbol{\phi}_k^T \mathbf{K} \boldsymbol{\phi}_k = \lambda_k`, and two modes with different eigenvalues are orthogonal for :math:`\mathbf{M}` and for :math:`\mathbf{K}` [Bathe1996]_ [Parlett1998]_.

Free body and constraints
~~~~~~~~~~~~~~~~~~~~~~~~~

Displacement boundary conditions are treated as homogeneous.
The constrained degrees of freedom are removed symmetrically from :math:`\mathbf{K}` and :math:`\mathbf{M}`, and the mode shapes are zero there.
If no displacement boundary condition is imposed, the structure is a free body.
Then :math:`\mathbf{K}` is singular and the six rigid-body modes are the first six eigenpairs, with :math:`\lambda = 0`.
They are three translations and three rotations.
The participation factors of the six rigid modes satisfy :math:`\sum_{k=1}^{6} \Gamma_{k,d}^2 = m_{tot}` for each direction :math:`d`, where :math:`m_{tot}` is the total mass.

Spectral shift
~~~~~~~~~~~~~~

The modes closest to a spectral shift :math:`\sigma` are computed.
The shift is given as a signed frequency by ``modalShiftFrequency``, and :math:`\sigma = \mathrm{sign}(f)\,(2 \pi f)^2`.
A negative frequency gives a negative shift :math:`\sigma = -\alpha`.
The shifted matrix :math:`\mathbf{K} - \sigma \mathbf{M} = \mathbf{K} + \alpha \mathbf{M}` is then symmetric positive definite, even for a free body.
Its linear solver is the one given in the ``LinearSolverParameters`` block of the solver.
The multigrid preconditioner of an elastic problem uses the rigid-body modes as near-null space.

Two eigensolvers are available, with the ``modalSolverType`` attribute.
Both are accessed through a common interface, so a new method can be added with one class and one entry in the factory.

Arnoldi eigensolver
~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~

The results are:

- A table in the log with the frequency, the eigenvalue, the residual and the participation factors :math:`\Gamma_{k,d} = \boldsymbol{\phi}_k^T \mathbf{M} \mathbf{e}_d` of each mode in the three directions :math:`d`.
- The arrays ``modalEigenvalues``, ``modalFrequencies``, ``modalResiduals`` and ``modalParticipationFactors`` of the solver, which are saved in the restart files.
  A frequency has the sign of its eigenvalue, so a rigid mode with a small negative eigenvalue has a small negative frequency.
- The nodal fields ``modeShape1``, ``modeShape2``, ..., one for each mode, which a VTK output writes.
  They are not saved in the restart files: after a restart, the frequencies and factors are available, and the mode shapes are zero until the modal analysis runs again.

The following limits apply.
The stiffness is the tangent stiffness at the current state, without the geometric (pre-stress) stiffness.
Only the lumped mass is available.
Contact, damping and body-force or traction loads are not used.
The ``Modal`` option is only available to a standalone solver, not to a solver that a coupled solver drives.

The modal analysis runs on the CPU, CUDA and HIP backends, and it does not use unified memory.
The vectors, the matrices and the preconditioner stay in the memory space of the backend during the eigensolve.
The dot products of the orthogonalization and of the Rayleigh-Ritz steps are computed in batches by device kernels, and only the small matrix of results is copied to the host.
Linear combinations of many vectors are done by one fused kernel.
The small dense eigenproblems of the projected matrices are solved on the host.

The examples :ref:`AdvancedExampleFreeFreeBeamModes` and :ref:`AdvancedExampleEigensolverComparison` verify this option against the Euler-Bernoulli beam theory and compare the two eigensolvers.

.. code-block:: xml

   <SolidMechanicsLagrangianFEM name="solid"
                                discretization="FE1"
                                targetRegions="{ Region }"
                                timeIntegrationOption="Modal"
                                modalNumModes="10"
                                modalShiftFrequency="-100"
                                modalTolerance="1e-8">
     <LinearSolverParameters solverType="cg"
                             preconditionerType="amg"
                             krylovTol="1e-10"/>
   </SolidMechanicsLagrangianFEM>

   <Events maxTime="1">
     <PeriodicEvent name="modalAnalysis" forceDt="1" target="/Solvers/solid"/>
   </Events>

Modal analysis references
~~~~~~~~~~~~~~~~~~~~~~~~~

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

In the preceding XML block, The `SolidMechanicsLagrangianFEM` is specified by the title of the subblock of the `Solvers` block.
The following attributes are supported in the input block for `SolidMechanicsLagrangianFEM`:

.. include:: /docs/sphinx/datastructure/SolidMechanicsLagrangianFEM.rst

The following data are allocated and used by the solver:

.. include:: /docs/sphinx/datastructure/SolidMechanicsLagrangianFEM_other.rst

Example
=========================

An example of a valid XML block is given here:

.. literalinclude:: ../../../../../inputFiles/solidMechanics/sedov_finiteStrain_smoke.xml
  :language: xml
  :start-after: <!-- SPHINX_SOLID_MECHANICS_SOLVER -->
  :end-before: <!-- SPHINX_SOLID_MECHANICS_SOLVER_END -->
