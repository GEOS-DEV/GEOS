
.. _SolutionStrategy:

#####################################
Solution strategy
#####################################

All  physics solvers share a common solution strategy for nonlinear time-dependent
problems. Here, we briefly describe the nonlinear solver and the timestepping
strategy employed.

Nonlinear solver
=================================
At each time-step, the nonlinear system of discrete residual equations, i.e.

.. math::
    r(x) = 0

is solved by employing the Newton-Raphson method. Here, :math:`x` is the vector of
primary unknowns. Thus, each physics solver is responsible for assembling
the Jacobian matrix :math:`J` containing the analytical derivatives of the residual
vector :math:`r` with respect to the primary variables. Then, at each Newton's iteration
:math:`\nu`, the following linear system is solved

.. math::
    J^{\nu} \delta x^{\nu+1} = -r^{\nu},

where, :math:`\delta x^{\nu+1}` is the Newton's update. This linear system can be
solved with a variety of different linear solvers described in :doc:`/coreComponents/linearAlgebra/docs/LinearSolvers`.
The Newton update, :math:`\delta x^{\nu+1}` is then applied to the primary variables:

..  math::
  x^{\nu+1} = x^{\nu} + \delta x^{\nu+1}.

This procedure is repeated until convergence.

Timestepping strategy
==================================




Parameters
============================

.. include:: /coreComponents/fileIO/schema/docs/NonLinearSolverParameters.rst

Example
================================
