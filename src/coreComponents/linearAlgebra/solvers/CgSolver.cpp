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

//  std::cout << "son qui\n\n\n\n" << std::endl;

  // Compute initial rk =  b - Ax
  m_operator.residual( x, b, r );

  // std::cout << "x0: \n" << x << std::endl;
  // std::cout << "b: \n" << b << std::endl;

  // std::cout << "residual: \n" << r << std::endl;

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
  watch.zero();

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
//    std::cout << "p("<<k<<"): \n" << p << std::endl;

    // Compute Ap
    m_operator.apply( p, Ap );
//    std::cout << "Ap("<<k<<"): \n" << Ap << std::endl;

    // compute alpha
    real64 const pAp = p.dot( Ap );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( pAp )
    real64 const alpha = tau / pAp;

    // Update x = x + alpha*p
    x.axpby( alpha, p, 1.0 );
//    std::cout << "x("<<k<<"+1): \n" << x << std::endl;

    // Update rk = rk - alpha*Ap
    r.axpby( -alpha, Ap, 1.0 );

//    std::cout << "r("<<k<<"+1): \n" << r << std::endl;


    // Keep the old value of tau
    tau_old = tau;
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
  arrayView1d< real64 const > const localX = x.values();
  arrayView1d< real64 const > const localY = y.values();
  arrayView1d< real64 > const localRes = res.open();
  GEOS_ERROR_IF(
    localX.size() != localY.size() || localX.size() != localRes.size(),
    "Cannot use CG solver with a non-symmetric system" );

  using POLICY = parallelDeviceAsyncPolicy<>;
  forAll< POLICY >( localX.size(), [alpha, beta, localX, localY, localRes] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localRes[ i ] = alpha * localX[ i ] + beta * localY[ i ];
  } );
  res.close();
}

template< typename Vector >
real64 dot( const Vector & x, const Vector & y )
{
  RAJA::ReduceSum< RAJA::cuda_reduce_atomic, real64 > vsum( 0.0 );
  arrayView1d< real64 const > const localX = x.values();
  arrayView1d< real64 const > const localY = y.values();

  using POLICY = parallelDeviceAsyncPolicy<>;
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

  // Define residual vector
  VectorTemp r = createTempVector( b );

//  std::cout << "son qui\n\n\n\n" << std::endl;

  // Compute initial rk =  b - Ax
  m_operator.residual( x, b, r );

  // std::cout << "x0: \n" << x << std::endl;
  // std::cout << "b: \n" << b << std::endl;

  // std::cout << "residual: \n" << r << std::endl;

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Preconditioning
  // VectorTemp z = createTempVector( x );

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
  for( k = 0; k <= m_params.krylov.maxIterations; ++k )
  {
    // real64 const rnorm = r.norm2();
    real64 const tau = dot( r, r );
    real64 const rnorm = std::sqrt( tau );
    m_residualNorms.emplace_back( rnorm );
    logProgress();

    // Convergence check on ||rk||/||b||
    if( rnorm <= absTol )
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    // Update z = Mr
    // m_precond.apply( r, z );

    // Compute beta
    // real64 const tau = r.dot( r );
    // real64 const tau = dot( r, r );
    real64 const beta = k > 0 ? tau / tau_old : 0.0;

    // Update p = r + beta*p
    // p.axpby( 1.0, r, beta );
    axpby( 1.0, r, beta, p, p );
//    std::cout << "p("<<k<<"): \n" << p << std::endl;

    // Compute Ap
    m_operator.apply( p, Ap );
//    std::cout << "Ap("<<k<<"): \n" << Ap << std::endl;

    // compute alpha
    // real64 const pAp = p.dot( Ap );
    real64 const pAp = dot( p, Ap );
    //GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( pAp )
    real64 const alpha = tau / pAp;

    // Update x = x + alpha*p
    // x.axpby( alpha, p, 1.0 );
    axpby( alpha, p, 1.0, x, x );
//    std::cout << "x("<<k<<"+1): \n" << x << std::endl;

    // Update rk = rk - alpha*Ap
    // r.axpby( -alpha, Ap, 1.0 );
    axpby( -alpha, Ap, 1.0, r, r );

//    std::cout << "r("<<k<<"+1): \n" << r << std::endl;


    // Keep the old value of tau
    tau_old = tau;
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
