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
 * @file GMRESsolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_GMRESSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_GMRESSOLVER_HPP_

#include "linearAlgebra/solvers/KrylovSolver.hpp"

namespace geosx
{

/**
 * @brief This class implements Generalized Minimized RESidual method
 *        (right-preconditioned) for monolithic and block linear operators.
 * @tparam VECTOR type of vectors this solver operates on.
 * @note  The notation is consistent with "Iterative Methods for
 *        Linear and Non-Linear Equations" from C.T. Kelley (1995)
 *        and "Iterative Methods for Sparse Linear Systems"
 *        from Y. Saad (2003).
 */
template< typename VECTOR >
class GMRESsolver : public KrylovSolver< VECTOR >
{
public:

  using Base = KrylovSolver< VECTOR >;
  using Vector = typename Base::Vector;

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Solver object constructor.
   */
  GMRESsolver( LinearOperator< Vector > const & A,
               LinearOperator< Vector > const & M,
               real64 const tolerance,
               localIndex const maxIterations,
               integer const verbosity = 0,
               localIndex const maxRestart = 100 );

  /**
   * @brief Virtual destructor.
   */
  virtual ~GMRESsolver() override;

  ///@}

  /**
   * @name KrylovSolver interface
   */
  ///@{

  virtual void solve( Vector const & b, Vector & x ) const override final;

  virtual string methodName() const override final
  {
    return "GMRES";
  };

  ///@}

protected:

  using VectorTemp = typename KrylovSolver< VECTOR >::VectorTemp;

  using Base::m_operator;
  using Base::m_precond;
  using Base::m_tolerance;
  using Base::m_maxIterations;
  using Base::m_logLevel;
  using Base::m_result;
  using Base::m_residualNorms;
  using Base::createTempVector;
  using Base::logProgress;
  using Base::logResult;

  /// Number of iterations needed to restart GMRES
  localIndex m_maxRestart;

  /// Storage for Krylov subspace vectors
  array1d< VectorTemp > m_kspace;

  /// Flag indicating whether kspace vectors have been created
  bool m_kspaceInitialized;
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_GMRESSOLVER_HPP_
