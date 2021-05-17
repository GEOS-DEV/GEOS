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
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"

namespace geosx
{

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
   * @brief Factory method for instantiating Krylov solver objects.
   * @param parameters solver parameters
   * @param matrix linear operator to solve (can be a matrix or a matrix-free operator)
   * @param precond preconditioning operator (must be set up by the user prior to calling solve()/apply())
   * @return an owning pointer to the newly instantiated solver
   */
  static std::unique_ptr< KrylovSolver< VECTOR > > create( LinearSolverParameters const & parameters,
                                                           LinearOperator< VECTOR > const & matrix,
                                                           LinearOperator< VECTOR > const & precond );

  /**
   * @brief Constructor.
   * @param [in] params parameters solver parameters
   * @param [in] matrix reference to the system matrix
   * @param [in] precond reference to the preconditioning operator
   */
  KrylovSolver( LinearSolverParameters params,
                LinearOperator< Vector > const & matrix,
                LinearOperator< Vector > const & precond );

  /**
   * @brief Virtual destructor
   */
  virtual ~KrylovSolver() override = default;

  /**
   * @brief Solve preconditioned system
   * @param [in] b system right hand side.
   * @param [inout] x system solution (input = initial guess, output = solution).
   */
  virtual void solve( Vector const & b, Vector & x ) const = 0;


  /**
   * @brief Apply operator to a vector.
   *
   * @param src Input vector (src).
   * @param dst Output vector (dst).
   */
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

  virtual localIndex numLocalRows() const override final
  {
    return m_operator.numLocalRows();
  }

  virtual localIndex numLocalCols() const override final
  {
    return m_operator.numLocalCols();
  }

  virtual MPI_Comm getComm() const override final
  {
    return m_operator.getComm();
  }

  /**
   * @brief @return parameters of the solver.
   */
  LinearSolverParameters const & parameters() const
  {
    return m_params;
  }

  /**
   * @brief @return the result of a linear solve.
   */
  LinearSolverResult const & result() const
  {
    return m_result;
  }

  /**
   * @brief Get convergence history of a linear solve.
   * @return array containing residual norms of every iteration (including initial)
   */
  arrayView1d< real64 const > history() const
  {
    return m_residualNorms;
  }

  /**
   * @brief Get name of the Krylov subspace method.
   * @return the abbreviated name of the method
   */
  virtual string methodName() const = 0;

private:

  ///@cond DO_NOT_DOCUMENT

  template< typename VEC >
  struct VectorStorageHelper
  {
    using type = VEC;

    static VEC createFrom( VEC const & src )
    {
      VEC v;
      v.createWithLocalSize( src.localSize(), src.getComm() );
      return v;
    }
  };

  template< typename VEC >
  struct VectorStorageHelper< BlockVectorView< VEC > >
  {
    using type = BlockVector< VEC >;

    static BlockVector< VEC > createFrom( BlockVectorView< VEC > const & src )
    {
      BlockVector< VEC > v( src.blockSize() );
      for( localIndex i = 0; i < src.blockSize(); ++i )
      {
        v.block( i ).createWithLocalSize( src.block( i ).localSize(), src.block( i ).getComm() );
      }
      return v;
    }
  };

  ///@endcond DO_NOT_DOCUMENT

protected:

  /// Alias for vector type that can be used for temporaries
  using VectorTemp = typename VectorStorageHelper< VECTOR >::type;

  /**
   * @brief Helper function to create temporary vectors based on a source vector.
   * @param src the source vector, whose size and parallel distribution will be used
   * @return the new vector
   *
   * The main purpose is to deal with BlockVector/View/Wrapper hierarchy.
   */
  static VectorTemp createTempVector( Vector const & src )
  {
    return VectorStorageHelper< VECTOR >::createFrom( src );
  }

  /**
   * @brief Output iteration progress (called by implementations).
   * @param iter  current iteration number
   * @param rnorm current residual norm
   */
  void logProgress( localIndex const iter, real64 const rnorm ) const
  {
    m_residualNorms[iter] = rnorm;
    if( m_params.logLevel >= 2 )
    {
      GEOSX_LOG_RANK_0( methodName() << " iteration " << iter << ": residual = " << rnorm );
    }
  }

  /**
   * @brief Output convergence result (called by implementations).
   */
  void logResult() const
  {
    if( m_params.logLevel >= 1 )
    {
      GEOSX_LOG_RANK_0( methodName() << ' ' <<
                        ( m_result.success() ? "converged" : "failed to converge" ) <<
                        " in " << m_result.numIterations << " iterations " <<
                        "(" << m_result.solveTime << " s)" );
    }
  }

  /// parameters of the solver
  LinearSolverParameters m_params;

  /// reference to the operator to be solved
  LinearOperator< Vector > const & m_operator;

  /// reference to the preconditioning operator
  LinearOperator< Vector > const & m_precond;

  /// results of a solve
  mutable LinearSolverResult m_result;

  /// Absolute residual norms at each iteration (if available)
  mutable array1d< real64 > m_residualNorms;
};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVSOLVER_HPP_
