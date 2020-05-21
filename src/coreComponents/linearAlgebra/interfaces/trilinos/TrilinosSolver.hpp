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
 * @file TrilinosSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSSOLVER_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

class EpetraVector;
class EpetraMatrix;
class LinearSolverParameters;

/**
 * @class TrilinosSolver
 * @brief This class creates and provides basic support for AztecOO, Amesos and ML libraries.
 */
class TrilinosSolver
{
public:

  /**
   * @brief Solver constructor, with parameter list reference
   *
   * @param[in] parameters structure containing linear solver parameters
   */
  TrilinosSolver( LinearSolverParameters const & parameters );

  /**
   * @brief Virtual destructor.
   *
   */
  ~TrilinosSolver();

  /**
   * @brief Solve system with an iterative solver.
   * @param[in,out] mat the matrix
   * @param[in,out] sol the solution
   * @param[in,out] rhs the right-hand side
   *
   * Solve Ax=b with A an EpetraMatrix, x and b EpetraVector.
   */
  void solve( EpetraMatrix & mat,
              EpetraVector & sol,
              EpetraVector & rhs );

  /**
   * @brief Number of krylov iterations to convergence.
   *
   * @note Value is meaningless if a direct solver is called and will return 1
   */
  integer iterations();

  /**
   * @brief Relative residual reduction.
   *
   * If the solve is successful, this value should be less than the target krylov
   * tolerance.  If the solver stagnates, however, it may be higher.
   *
   * @return Reduction value
   * @note Value is meaningless if a direct solver is called and will return machine precision;
   */
  real64  reduction();

  /**
   * @brief Setup time (in seconds) for preconditioners and/or direct factorizations
   * @return Setup time
   */
  real64  setupTime();

  /**
   * @brief Solve time (in seconds) exclusive of setup costs
   * @return Solve time
   */
  real64  solveTime();

  /**
   * @brief Total time (in seconds), the sum of setupTime() and solveTime()
   * @return Total time
   */
  real64  totalTime();

private:

  LinearSolverParameters const & m_parameters;

  void solve_direct( EpetraMatrix & mat,
                     EpetraVector & sol,
                     EpetraVector & rhs );

  void solve_krylov( EpetraMatrix & mat,
                     EpetraVector & sol,
                     EpetraVector & rhs );

  integer m_iterations;
  real64 m_reduction;
  real64 m_setupTime;
  real64 m_solveTime;
};

} // end geosx namespace

#endif /* TRILINOSSOLVER_HPP_ */
