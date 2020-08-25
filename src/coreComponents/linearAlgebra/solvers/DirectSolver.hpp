/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_DIRECTSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_DIRECTSOLVER_HPP_

#include "linearAlgebra/utilities/BlockVectorView.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"
#include "linearAlgebra/utilities/BlockOperatorView.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"

namespace geosx
{

/**
 * @brief Base class for Direct solvers
 * @tparam LAI interface handled by this solver
 */
template< typename LAI >
class DirectSolver
{
public:

  /// Alias of linear algebra parallel matrix
  using Matrix = typename LAI::ParallelMatrix;

  /// Alias of linear algebra parallel vector
  using Vector = typename LAI::ParallelVector;

  /**
   * @brief Constructor.
   * @param [in] params linear solver parameters.
   */
  DirectSolver( LinearSolverParameters const params );

  /**
   * @brief Default destructor
   */
  ~DirectSolver() = default;

  /**
   * @brief Solve preconditioned system
   * @param [in] matrix system matrix.
   * @param [in] b system right hand side.
   * @param [out] x system solution.
   */
  void solve( Matrix & matrix,
              Vector & b,
              Vector & x ) const;

  /**
   * @brief Get result of a linear solve.
   * @return struct containing status and various statistics of the last solve
   */
  LinearSolverResult const & result() const
  {
    return m_result;
  }

protected:

  /// solver verbosity level
  LinearSolverParameters m_params;

  /// results of a solve
  mutable LinearSolverResult m_result;
};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_DIRECTSOLVER_HPP_
