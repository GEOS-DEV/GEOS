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
 * @file BiCGSTABsolver.hpp
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
class BiCGSTABsolver : public KrylovSolver< VECTOR >
{
public:

  using Vector = typename KrylovSolver< VECTOR >::Vector;

  //! @name Constructor/Destructor Methods
  //@{

  /**
   * @brief Solver object constructor.
   */
  BiCGSTABsolver( LinearOperator< Vector > const & A,
                  LinearOperator< Vector > const & M,
                  real64 const tolerance,
                  localIndex const maxIterations,
                  integer const verbosity = 0 );

  /**
   * @brief Virtual destructor.
   */
  virtual ~BiCGSTABsolver() override;

  //@}

  //! @name KrylovSolver interface
  //@{

  virtual void solve( Vector const & b, Vector & x ) const override final;

  //@}

protected:

  using VectorTemp = typename KrylovSolver< VECTOR >::VectorTemp;

  using KrylovSolver< VECTOR >::m_operator;
  using KrylovSolver< VECTOR >::m_precond;
  using KrylovSolver< VECTOR >::m_tolerance;
  using KrylovSolver< VECTOR >::m_maxIterations;
  using KrylovSolver< VECTOR >::m_verbosity;
  using KrylovSolver< VECTOR >::m_numIterations;
  using KrylovSolver< VECTOR >::m_residualNormVector;
  using KrylovSolver< VECTOR >::m_convergenceFlag;

};

} // namespace GEOSX

#endif /*GEOSX_LINEARALGEBRA_SOLVERS_BICGSTABSOLVER_HPP_ */
