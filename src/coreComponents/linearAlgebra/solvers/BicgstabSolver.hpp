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
 * @file BicgstabSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_BICGSTABSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_BICGSTABSOLVER_HPP_

#include "linearAlgebra/solvers/KrylovSolver.hpp"

namespace geosx
{

/**
 * @brief This class implements Bi-Conjugate Gradient Stabilized method
 *        for monolithic and block linear operators.
 * @tparam VECTOR type of vectors this solver operates on.
 * @note  The notation is consistent with "Iterative Methods for
 *        Linear and Non-Linear Equations" from C.T. Kelley (1995)
 *        and "Iterative Methods for Sparse Linear Systems"
 *        from Y. Saad (2003).
 */
template< typename VECTOR >
class BicgstabSolver : public KrylovSolver< VECTOR >
{
public:

  /// Alias for base type
  using Base = KrylovSolver< VECTOR >;

  /// Alias for template parameter
  using Vector = typename Base::Vector;

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Constructor.
   * @param [in] params parameters for the solver
   * @param [in] A reference to the system matrix.
   * @param [in] M reference to the preconditioning operator.
   */
  BicgstabSolver( LinearSolverParameters params,
                  LinearOperator< Vector > const & A,
                  LinearOperator< Vector > const & M );

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
  virtual void solve( Vector const & b, Vector & x ) const override final;

  virtual string methodName() const override final
  {
    return "BiCGStab";
  };

  ///@}

protected:

  /// Alias for vector type that can be used for temporaries
  using VectorTemp = typename KrylovSolver< VECTOR >::VectorTemp;

  using Base::m_params;
  using Base::m_operator;
  using Base::m_precond;
  using Base::m_result;
  using Base::m_residualNorms;
  using Base::createTempVector;
  using Base::logProgress;
  using Base::logResult;

};

} // namespace GEOSX

#endif /*GEOSX_LINEARALGEBRA_SOLVERS_BICGSTABSOLVER_HPP_ */
