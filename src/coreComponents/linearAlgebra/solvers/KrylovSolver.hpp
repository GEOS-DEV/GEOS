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

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_

#include "linearAlgebra/utilities/BlockVectorView.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"
#include "linearAlgebra/utilities/BlockOperatorView.hpp"

namespace geosx
{

template< typename Vector > class LinearOperator;
template< typename Vector > class BlockVectorView;

/**
 * @brief Base class for Krylov solvers
 * @tparam VECTOR type of vector handled by this solver
 */
template< typename VECTOR >
class KrylovSolver : public LinearOperator< VECTOR >
{
public:

  /// Base type
  using Base = LinearOperator< VECTOR >;

  /// Alias for template parameter
  using Vector = typename Base::Vector;

  /**
   * @brief Constructor
   */
  KrylovSolver( LinearOperator< Vector > const & A,
                LinearOperator< Vector > const & M,
                real64 const tolerance,
                localIndex const maxIterations,
                integer const verbosity )
    : Base(),
    m_operator( A ),
    m_precond( M ),
    m_tolerance( tolerance ),
    m_maxIterations( maxIterations ),
    m_verbosity( verbosity )
  {
    GEOSX_ERROR_IF_LT_MSG( m_maxIterations, 0, "Krylov solver: max number of iteration must be non-negative." );
    GEOSX_LAI_ASSERT_EQ( m_operator.numGlobalRows(), m_precond.numGlobalRows() );
    GEOSX_LAI_ASSERT_EQ( m_operator.numGlobalCols(), m_precond.numGlobalCols() );
  }

  /**
   * @brief Virtual destructor
   */
  virtual ~KrylovSolver() override = default;

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with CG
   * using monolithic GEOSX matrices.
   *
   * @param A system matrix.
   * @param x system solution (input = initial guess, output = solution).
   * @param b system right hand side.
   * @param M preconditioner.
   */
  virtual void solve( Vector const & b, Vector & x ) const = 0;

  virtual void apply( Vector const & src, Vector & dst ) const override final
  {
    solve( src, dst );
  }

  virtual globalIndex numGlobalRows() const override final
  {
    return m_operator.numGlobalRows();
  }

  virtual globalIndex numGlobalCols() const override final
  {
    return m_operator.numGlobalCols();
  }

  localIndex numIterations( ) const { return m_numIterations; };

  arrayView1d< const real64 > residualNormVector( ) const { return m_residualNormVector.toViewConst(); };

  bool convergenceFlag( ) const { return m_convergenceFlag; };

private:

  // Helper struct to get temporary stored vector type from vector view types
  template< typename VEC >
  struct VectorStorageHelper
  {
    using type = VEC;
  };

  template< typename VEC >
  struct VectorStorageHelper< BlockVectorView< VEC > >
  {
    using type = BlockVector< VEC >;
  };

protected:

  /// Alias for vector type that can be used for temporaries
  using VectorTemp = typename VectorStorageHelper< VECTOR >::type;

  /// reference to the operator to be solved
  LinearOperator< Vector > const & m_operator;

  /// reference to the preconditioning operator
  LinearOperator< Vector > const & m_precond;

  /// relative residual norm reduction tolerance
  real64 m_tolerance;

  /// maximum number of Krylov iterations
  localIndex m_maxIterations;

  /// solver verbosity level
  integer m_verbosity;

  /// actual number if Krylov iterations
  mutable localIndex m_numIterations = 0;

  /// residual norm vector
  mutable array1d< real64 > m_residualNormVector;

  /// convergence flag
  mutable bool m_convergenceFlag;

};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_
