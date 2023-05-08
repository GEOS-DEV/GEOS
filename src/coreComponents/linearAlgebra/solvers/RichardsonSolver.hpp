/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file RichardsonSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_RICHARDSONSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_RICHARDSONSOLVER_HPP_

#include "linearAlgebra/solvers/KrylovSolver.hpp"

namespace geos
{

/**
 * @brief Implements right-preconditioned modified Richardson iteration.
 * @tparam VECTOR type of vectors this solver operates on
 * @note Richardson is not a Krylov subspace method, but
 *       for convenience inherits from KrylovSolver; that
 *       class should really be renamed to IterativeSolver.
 */
template< typename VECTOR >
class RichardsonSolver final : public KrylovSolver< VECTOR >
{
public:

  /// Alias for the base type
  using Base = KrylovSolver< VECTOR >;

  /// Alias for the vector type
  using Vector = typename Base::Vector;

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Solver object constructor.
   * @param[in] params  parameters for the solver
   * @param[in] matrix  reference to the system matrix
   * @param[in] precond reference to the preconditioning operator
   */
  RichardsonSolver( LinearSolverParameters params,
                    LinearOperator< Vector > const & matrix,
                    LinearOperator< Vector > const & precond );

  ///@}

  /**
   * @name KrylovSolver interface
   */
  ///@{

  /**
   * @brief Solve preconditioned system
   * @param [in] b system right hand side.
   * @param [inout] x system solution (input = initial guess, output = solution).
   */
  virtual void solve( Vector const & b, Vector & x ) const override;

  virtual string methodName() const override
  {
    return "Richardson";
  };

  ///@}

protected:

  /// Alias for vector type that can be used for temporaries
  using VectorTemp = typename KrylovSolver< VECTOR >::VectorTemp;

  using Base::m_params;
  using Base::m_operator;
  using Base::m_precond;
  using Base::m_residualNorms;
  using Base::m_result;
  using Base::createTempVector;
  using Base::logProgress;
  using Base::logResult;

private:

  real64 m_omega;
};

} // geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_RICHARDSONSOLVER_HPP_
