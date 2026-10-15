/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file LinearSolverBase.hpp
 */
#ifndef GEOS_LINEARSOLVERBASE_HPP
#define GEOS_LINEARSOLVERBASE_HPP

#include "PreconditionerBase.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"
#include "common/Stopwatch.hpp"

#include <string>

namespace geos
{

/**
 * @brief Execution metadata associated with a linear-solver invocation.
 *
 * This context lets solver backends annotate per-solve setup, logging, and
 * statistics with the timestep/newton state that produced the linear system.
 */
struct LinearSolverExecutionContext
{
  /// Human-readable name of the owning solver.
  std::string solverName;
  /// Current outer cycle index, or `-1` when unavailable.
  integer cycleNumber = -1;
  /// Current timestep-attempt index, or `-1` when unavailable.
  integer timeStepAttempt = -1;
  /// Current configuration-attempt index, or `-1` when unavailable.
  integer configurationAttempt = -1;
  /// Current nonlinear-iteration index, or `-1` when unavailable.
  integer nonlinearIteration = -1;
  /// Timestamp of the last system assembly/setup used by this solve.
  Timestamp systemSetupTimestamp = 0;
};

/**
 * @brief Simple interface for linear solvers that allows to extract solution results.
 * @tparam LAI linear algebra interface type
 */
template< typename LAI >
class LinearSolverBase : public PreconditionerBase< LAI >
{
public:

  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param params solver parameters
   */
  explicit LinearSolverBase( LinearSolverParameters params )
    : m_params( std::move( params ) )
  {}

  /**
   * @brief Solve preconditioned system
   * @param [in] rhs system right hand side.
   * @param [in,out] sol system solution (input = initial guess, output = solution).
   */
  virtual void solve( Vector const & rhs, Vector & sol ) const = 0;

  /**
   * @brief Provide contextual execution metadata for backends that need it.
   * @param context The current linear-solve execution context.
   */
  virtual void setExecutionContext( LinearSolverExecutionContext const & context )
  {
    GEOS_UNUSED_VAR( context );
  }

  /**
   * @brief Supply near-null-space vectors used to construct physics-aware preconditioners.
   * @param nearNullKernel Distributed vectors spanning the near null space.
   */
  virtual void setNearNullKernel( arrayView1d< Vector const > const & nearNullKernel )
  {
    GEOS_UNUSED_VAR( nearNullKernel );
  }

  /**
   * @brief @return parameters of the solver.
   */
  LinearSolverParameters const & parameters() const
  {
    return m_params;
  }

  /**
   * @brief @return result of the most recent solve.
   */
  LinearSolverResult const & result() const
  {
    return m_result;
  }

protected:

  /// Parameters for the solver
  LinearSolverParameters m_params;

  /// Result of most recent solve (status, timings)
  mutable LinearSolverResult m_result;
};

}

#endif //GEOS_LINEARSOLVERBASE_HPP
