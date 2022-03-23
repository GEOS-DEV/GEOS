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
 * @file ExpBicgstabSolver.cpp
 *
 */

#include "ExpBicgstabSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "common/LinearOperator.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geosx
{


template< typename VECTOR >
ExpBicgstabSolver< VECTOR >::ExpBicgstabSolver( LinearSolverParameters params,
                                                LinearOperator< Vector > const & A,
                                                LinearOperator< Vector > const & M )
  : KrylovSolver< VECTOR >( std::move( params ), A, M )
{}

template< typename VECTOR >
void ExpBicgstabSolver< VECTOR >::solve( Vector const & b,
                                         Vector & x ) const
{
  Stopwatch watch;

  // Define vectors
  VectorTemp r( x );

  // Compute initial rk
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Define temporary vectors
  VectorTemp S  = createTempVector( r );
  VectorTemp Y  = createTempVector( r );
  VectorTemp Q  = createTempVector( r );
  VectorTemp Q2 = createTempVector( r );
  VectorTemp P2 = createTempVector( r );
  VectorTemp R2 = createTempVector( r );
  VectorTemp S2 = createTempVector( r );
  VectorTemp W  = createTempVector( r );
  VectorTemp Z  = createTempVector( r );
  VectorTemp W2 = createTempVector( r );
  VectorTemp Z2 = createTempVector( r );
  VectorTemp T  = createTempVector( r );
  VectorTemp V  = createTempVector( r );

  // Define vectors
  VectorTemp r0( r );

  // Define scalars and reinitialize some
  real64 rho = rnorm0 * rnorm0; // same as r.dot( r0 );
  
  m_precond.apply( r, R2 );
  m_operator.apply( R2, W );
  
  real64 d2 = W.dot( r0 );
    
  m_precond.apply( W, W2 );
  m_operator.apply( W2, T );

  real64 alpha = rho / d2;
  real64 beta = 0.0;
  real64 omega;
  real64 d1, d3;

  // Initialize iteration state
  m_result.status = LinearSolverResult::Status::NotConverged;
  m_residualNorms.clear();

  integer & k = m_result.numIterations;
  for( k = 0; k <= m_params.krylov.maxIterations; ++k )
  {
    real64 const rnorm = r.norm2();
    m_residualNorms.emplace_back( rnorm );
    logProgress();

    // Convergence check on ||rk||/||b||
    if( rnorm <= absTol )
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    // 
    if ( k == 0 )
    {
      R2.copy( P2 );
      W.copy( S );
      W2.copy( S2 );
      T.copy( Z );
    }
    else
    {
      P2.axpby( 1.0, R2, beta );    // TODO
      P2.axpy( -beta * omega, S2 ); // add method axpbypcz
      S.axpby( 1.0, W, beta );
      S.axpy( -beta * omega, Z );
      S2.axpby( 1.0, W2, beta );
      S2.axpy( -beta * omega, Z2 );
      Z.axpby( 1.0, T, beta );
      Z.axpy( -beta * omega, V );
    }
    
    r.copy( Q );
    Q.axpy ( -alpha, x );
    R2.copy( Q2 );
    Q2.axpy( -alpha, S2 );
    W.copy( Y );
    Y.axpy( -alpha, Z );
    
    d1 = Q.dot( Y );
    d2 = Y.dot( Y );

    if ( isZero( d2 ) )
    {
      d1 = Q.dot( Q );  
      GEOSX_KRYLOV_BREAKDOWN_IFNOT_ZERO( d1 )
      x.axpy( alpha, P2 );
      ++k;
      m_residualNorms.emplace_back( 0.0 );
      logProgress();
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    omega = d1 / d2;

    x.axpy( alpha, P2 );
    x.axpy( omega, Q2 );

    Q.copy( r );
    r.axpy ( -omega, Y );
    W2.copy( R2 );
    R2.axpy( -alpha, Z2 );
    R2.axpby( 1.0, Q2, -omega );
    T.copy( W );
    W.axpy( -alpha, V );
    W.axpby( 1.0, Y, -omega );

    real64 rho_old = rho;

    rho = r.dot( r0 );
    d1 = S.dot( r0 );
    d2 = W.dot( r0 );
    d3 = Z.dot( r0 );

    m_precond.apply( W, W2 );
    m_operator.apply( W2, T );

    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( rho_old );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( omega );
    beta = ( rho / rho_old ) * ( alpha / omega );
    alpha = rho / ( d2 + beta *d1 - beta * omega * d3 );
  }


  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class ExpBicgstabSolver< TrilinosInterface::ParallelVector >;
template class ExpBicgstabSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class ExpBicgstabSolver< HypreInterface::ParallelVector >;
template class ExpBicgstabSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class ExpBicgstabSolver< PetscInterface::ParallelVector >;
template class ExpBicgstabSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

} //namespace geosx
