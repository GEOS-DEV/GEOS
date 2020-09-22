/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TrilinosTpetraSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRASOLVER_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRASOLVER_HPP_

#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"

namespace geosx
{

class DofManager;
class TpetraVector;
class TpetraMatrix;

/**
 * @brief Wrapper for Trilinos/Tpetra-based direct and interative linear solvers.
 */
class TrilinosTpetraSolver
{
public:

  /**
   * @brief Solver constructor, with parameter list reference.
   * @param[in] parameters structure containing linear solver parameters
   */
  TrilinosTpetraSolver( LinearSolverParameters parameters );

  /**
   * @brief Virtual destructor.
   *
   */
  ~TrilinosTpetraSolver();

  /**
   * @brief Solve a linear system Ax=b with an iterative solver.
   * @param[in,out] mat the matrix
   * @param[in,out] sol the solution
   * @param[in,out] rhs the right-hand side
   * @param dofManager the Degree-of-Freedom manager associated with matrix
   */
  void solve( TpetraMatrix & mat,
              TpetraVector & sol,
              TpetraVector & rhs,
              DofManager const * const dofManager = nullptr );

  /**
   * @brief Get the result of previous solve.
   * @return struct with last solve stats
   */
  LinearSolverResult const & result()
  {
    return m_result;
  }

private:

  LinearSolverParameters m_parameters;
  LinearSolverResult m_result;

  void solve_direct( TpetraMatrix & mat,
                     TpetraVector & sol,
                     TpetraVector & rhs );

  void solve_krylov( TpetraMatrix & mat,
                     TpetraVector & sol,
                     TpetraVector & rhs );
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSTPETRASOLVER_HPP_
