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
 * @file CGsolver.cpp
 *
 */

#include "CGsolver.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
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
// Empty constructor.
template< typename VECTOR >
CGsolver< VECTOR >::CGsolver( LinearOperator< Vector > const & A,
                              LinearOperator< Vector > const & M,
                              real64 const tolerance,
                              localIndex const maxIterations,
                              integer const verbosity )
  : KrylovSolver< VECTOR >( A, M, tolerance, maxIterations, verbosity )
{}

// ----------------------------
// Destructor
// ----------------------------
template< typename VECTOR >
CGsolver< VECTOR >::~CGsolver() = default;

// ----------------------------
// Monolithic CG solver
// ----------------------------
template< typename VECTOR >
void CGsolver< VECTOR >::solve( Vector const & b,
                                Vector & x ) const

{
  // Shortcuts for operators
  LinearOperator< VECTOR > const & A = m_operator;
  LinearOperator< VECTOR > const & M = m_precond;

  // Get the global size and resize residual norm vector
  localIndex const N = m_maxIterations;
  m_residualNormVector.resize( N + 1 );

  // Get the norm of the right hand side
  real64 const normb = b.norm2();

  // Define residual vector
  VectorTemp rk( x );

  // Compute initial rk =  b - Ax
  A.residual( x, b, rk );
  m_residualNormVector[0]= rk.norm2();

  // Preconditioning
  VectorTemp zk( x );
  M.apply( rk, zk );

  // pk = zk
  VectorTemp pk( zk );
  VectorTemp Apk( zk );

  // Declare older vectors
  VectorTemp rkold( rk );
  VectorTemp zkold( zk );

  localIndex k;
  for( k = 0; k < N; k++ )
  {
    // Compute Apk
    A.apply( pk, Apk );

    // compute alpha
    real64 const alpha = rk.dot( zk ) / pk.dot( Apk );

    // Update x = x + alpha*ph
    x.axpby( alpha, pk, 1.0 );

    // Update rk = rk - alpha*Apk
    rkold.copy( rk );
    zkold.copy( zk );
    rk.axpby( -alpha, Apk, 1.0 );
    m_residualNormVector[ k+1 ] = rk.norm2();

    // Convergence check on ||rk||/||b||
    if( m_residualNormVector[ k+1 ] / normb < 1e-8 )
    {
      break;
    }

    // Update zk = Mrk
    M.apply( rk, zk );

    // Compute beta
    real64 const beta = zk.dot( rk ) / zkold.dot( rkold );

    // Update pk = pk + beta*zk
    pk.axpby( 1.0, zk, beta );
  }

  // Convergence statistics
  m_numIterations = k;
  m_convergenceFlag = k < N;
  m_residualNormVector.resize( m_numIterations + 1 );

  if( m_verbosity >= 1 )
  {
    GEOSX_LOG_RANK_0( "CG " << (k < N ? "converged" : "did not converge") << " in " << k << " iterations." );
  }
}

// END_RST_NARRATIVE

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class CGsolver< TrilinosInterface::ParallelVector >;
template class CGsolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class CGsolver< HypreInterface::ParallelVector >;
template class CGsolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class CGsolver< PetscInterface::ParallelVector >;
template class CGsolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

} //namespace geosx
