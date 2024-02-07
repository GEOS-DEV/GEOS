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
 * @file CGsolver.cpp
 *
 */


#include "CgSolver.hpp"

#include "common/Stopwatch.hpp"
#include "common/TimingMacros.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "common/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geos
{

// BEGIN_RST_NARRATIVE CGsolver.rst
// ==============================
// CG Solver
// ==============================
// Implementation of the Preconditioned Conjugate Gradient (CG) algorithm.
// The notation is consistent with "Iterative Methods for
// Linear and Non-Linear Equations" from C.T. Kelley (1995)
// and "Iterative Methods for Sparse Linear Systems" from Y. Saad (2003).


// ----------------------------
// Constructor
// ----------------------------
template< typename VECTOR >
CgSolver< VECTOR >::CgSolver( LinearSolverParameters params,
                              LinearOperator< Vector > const & A,
                              LinearOperator< Vector > const & M )
  : KrylovSolver< VECTOR >( std::move( params ), A, M )
{
  GEOS_ERROR_IF( !m_params.isSymmetric, "Cannot use CG solver with a non-symmetric system" );
}

// ----------------------------
// Solve method
// ----------------------------
template< typename VECTOR >
void CgSolver< VECTOR >::solve( Vector const & b, Vector & x ) const
{
  GEOS_MARK_FUNCTION;

  Stopwatch watch;

  // Define residual vector
  VectorTemp r = createTempVector( b );

  // Compute initial rk =  b - Ax
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Preconditioning
  VectorTemp z = createTempVector( x );

  // Search direction
  VectorTemp p = createTempVector( z );
  VectorTemp Ap = createTempVector( z );
  p.zero();

  // Keep old value of preconditioned residual norm
  real64 tau_old = 0.0;

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

    // Update z = Mr
    m_precond.apply( r, z );

    // Compute beta
    real64 const tau = z.dot( r );
    real64 const beta = k > 0 ? tau / tau_old : 0.0;

    // Update p = z + beta*p
    p.axpby( 1.0, z, beta );

    // Compute Ap
    m_operator.apply( p, Ap );

    // compute alpha
    real64 const pAp = p.dot( Ap );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( pAp )
    real64 const alpha = tau / pAp;

    // Update x = x + alpha*p
    x.axpby( alpha, p, 1.0 );

    // Update rk = rk - alpha*Ap
    r.axpby( -alpha, Ap, 1.0 );

    // Keep the old value of tau
    tau_old = tau;
  }

  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// END_RST_NARRATIVE

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class CgSolver< TrilinosInterface::ParallelVector >;
template class CgSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class CgSolver< HypreInterface::ParallelVector >;
template class CgSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class CgSolver< PetscInterface::ParallelVector >;
template class CgSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

// ----------------------------
// Constructor
// ----------------------------
template< typename VECTOR >
UnprecCgSolver< VECTOR >::UnprecCgSolver( LinearSolverParameters params,
                                          LinearOperator< Vector > const & A )
  : KrylovSolver< VECTOR >( std::move( params ), A, A )
{
  GEOS_ERROR_IF( !m_params.isSymmetric, "Cannot use CG solver with a non-symmetric system" );
}

template< typename Vector >
void axpby( real64 alpha, const Vector & x, real64 beta, const Vector & y, Vector & res )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const localX = x.values();
  arrayView1d< real64 const > const localY = y.values();
  arrayView1d< real64 > const localRes = res.open();
  GEOS_ERROR_IF(
    localX.size() != localY.size() || localX.size() != localRes.size(),
    "Cannot use CG solver with a non-symmetric system" );

  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localX.size(), [alpha, beta, localX, localY, localRes] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localRes[ i ] = alpha * localX[ i ] + beta * localY[ i ];
  } );
  res.close();
}


template< typename Vector >
void xpby( const Vector & x, real64 beta, Vector & res )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const localX = x.values();
  arrayView1d< real64 > const localRes = res.open();
  GEOS_ERROR_IF( localX.size() != localRes.size(), "Cannot use CG solver with a non-symmetric system" );

  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localX.size(), [beta, localX, localRes] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localRes[ i ] = localX[ i ] + beta * localRes[ i ];
  } );
  res.close();
}


template< typename Vector >
real64 axpby_dot( real64 alpha, const Vector & x, real64 beta, const Vector & y, Vector & res )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const localX = x.values();
  arrayView1d< real64 const > const localY = y.values();
  arrayView1d< real64 > const localRes = res.open();
  GEOS_ERROR_IF(
    localX.size() != localY.size() || localX.size() != localRes.size(),
    "Cannot use CG solver with a non-symmetric system" );


  RAJA::ReduceSum< parallelDeviceReduce, real64 > sum( 0.0 );

  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localX.size(), [alpha, beta, localX, localY, localRes, sum] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localRes[ i ] = alpha * localX[ i ] + beta * localY[ i ];
    sum += localRes[ i ] * localRes[ i ];
  } );
  res.close();
  return static_cast< real64 >(sum.get());
}


template< typename Vector >
real64 axpby2_dot( real64 const alpha0, const Vector & x0, real64 const beta0, const Vector & y0, Vector & res0,
                   real64 const alpha1, const Vector & x1, real64 const beta1, const Vector & y1, Vector & res1 )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const localX0 = x0.values();
  arrayView1d< real64 const > const localY0 = y0.values();
  arrayView1d< real64 > const localRes0 = res0.open();

  arrayView1d< real64 const > const localX1 = x1.values();
  arrayView1d< real64 const > const localY1 = y1.values();
  arrayView1d< real64 > const localRes1 = res1.open();

  GEOS_ERROR_IF( localX0.size() != localY0.size() ||
                 localX0.size() != localRes0.size() ||
                 localX0.size() != localX1.size() ||
                 localX1.size() != localY1.size() ||
                 localX1.size() != localRes1.size(),
                 "Cannot use CG solver with a non-symmetric system" );


  RAJA::ReduceSum< parallelDeviceReduce, real64 > sum( 0.0 );

  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localX0.size(), [ alpha0, beta0, localX0, localY0, localRes0,
                                      alpha1, beta1, localX1, localY1, localRes1, sum] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localRes0[ i ] = alpha0 * localX0[ i ] + beta0 * localY0[ i ];
    localRes1[ i ] = alpha1 * localX1[ i ] + beta1 * localY1[ i ];
    sum += localRes1[ i ] * localRes1[ i ];
  } );
  res0.close();
  res1.close();
  return static_cast< real64 >(sum.get());
}


template< typename Vector >
real64 axpy2_dot( real64 const alpha0, const Vector & x0, Vector & res0,
                  real64 const alpha1, const Vector & x1, Vector & res1 )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< real64 const > const localX0 = x0.values();
  arrayView1d< real64 > const localRes0 = res0.open();

  arrayView1d< real64 const > const localX1 = x1.values();
  arrayView1d< real64 > const localRes1 = res1.open();

  GEOS_ERROR_IF( localX0.size() != localRes0.size() ||
                 localX0.size() != localX1.size() ||
                 localX1.size() != localRes1.size(),
                 "Cannot use CG solver with a non-symmetric system" );


  RAJA::ReduceSum< parallelDeviceReduce, real64 > sum( 0.0 );

  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localX0.size(), [ alpha0, localX0, localRes0,
                                      alpha1, localX1, localRes1, sum] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localRes0[ i ] = localRes0[ i ] + alpha0 * localX0[ i ];
    localRes1[ i ] = localRes1[ i ] + alpha1 * localX1[ i ];
    sum += localRes1[ i ] * localRes1[ i ];
  } );
  res0.close();
  res1.close();
  return static_cast< real64 >(sum.get());
}

template< typename Vector >
real64 dot2( const Vector & x, const Vector & y )
{
  GEOS_MARK_FUNCTION;
  RAJA::ReduceSum< parallelDeviceReduce, real64 > vsum( 0.0 );
  arrayView1d< real64 const > const localX = x.values();
  arrayView1d< real64 const > const localY = y.values();

  using POLICY = parallelDeviceAsyncPolicy< 1024 >;
  forAll< POLICY >( localX.size(), [vsum, localX, localY] GEOS_HOST_DEVICE ( localIndex const i )
  {
    vsum += localX[ i ] * localY[ i ];
  } );
  return static_cast< real64 >(vsum.get());
}

template< typename VECTOR >
void UnprecCgSolver< VECTOR >::solve( Vector const & b, Vector & x ) const
{
  GEOS_MARK_FUNCTION;

  Stopwatch watch;
  Stopwatch watch2;

  m_result.minIterTime = 1e99;

  // Define residual vector
  VectorTemp r = createTempVector( b );

  // Compute initial rk =  b - Ax
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;


  // Search direction
  VectorTemp p = createTempVector( r );
  VectorTemp Ap = createTempVector( r );
  p.zero();

  // Keep old value of preconditioned residual norm
  real64 tau_old = 0.0;

  // Initialize iteration state
  m_result.status = LinearSolverResult::Status::NotConverged;
  m_residualNorms.clear();
  watch.zero();

  integer & k = m_result.numIterations;
  real64 tau = dot2( r, r );
  for( k = 0; k <= m_params.krylov.maxIterations; ++k )
  {
    watch2.zero();
    real64 const rnorm = std::sqrt( tau );
    m_residualNorms.emplace_back( rnorm );
    logProgress();

    // Convergence check on ||rk||/||b||
    if( rnorm <= absTol && k > 0 )
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    // Compute beta
    real64 const beta = k > 0 ? tau / tau_old : 0.0;

    // Update p = r + beta*p
    // p.axpby( 1.0, r, beta );
    xpby( r, beta, p );

    // Compute Ap
    m_operator.apply( p, Ap );

    // compute alpha
    real64 const pAp = dot2( p, Ap );
    //GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( pAp )
    real64 const alpha = tau / pAp;

    // Keep the old value of tau
    tau_old = tau;

    #define GC_FUSION_LEVEL 3
#if GC_FUSION_LEVEL==3
    tau = axpy2_dot( +alpha, p, x,
                     -alpha, Ap, r );
#elif GC_FUSION_LEVEL==2
    tau = axpby2_dot( +alpha, p, 1.0, x, x,
                      -alpha, Ap, 1.0, r, r );
#elif GC_FUSION_LEVEL==1
    // Update x = x + alpha*p
    axpby( alpha, p, 1.0, x, x );

    tau = axpby_dot( -alpha, Ap, 1.0, r, r );
#else
    // Update x = x + alpha*p
    axpby( alpha, p, 1.0, x, x );

    // Update rk = rk - alpha*Ap
    axpby( -alpha, Ap, 1.0, r, r );

    tau = dot2( r, r );
#endif

    real64 iterTime = watch2.elapsedTime();
    m_result.minIterTime = LvArray::math::min( iterTime, m_result.minIterTime );

  }
  // std::cout << "iter: " << k << std::endl;
  // std::cout << "solution: \n" << x << std::endl;

  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// END_RST_NARRATIVE

// -----------------------
// Explicit Instantiations
// -----------------------

#ifdef GEOSX_USE_HYPRE
template class UnprecCgSolver< HypreInterface::ParallelVector >;
#endif

} //namespace geos
