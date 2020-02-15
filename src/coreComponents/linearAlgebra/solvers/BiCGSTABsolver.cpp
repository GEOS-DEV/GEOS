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
 * @file BiCGSTABsolver.cpp
 *
 */

#include "BiCGSTABsolver.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockOperatorView.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
{

// BEGIN_RST_NARRATIVE BiCGSTABsolver.rst
// ==============================
// BiCGSTAB Solver
// ==============================
// Implementation of the BiCGSTAB algorithm.
// The notation is consistent with "Iterative Methods for
// Linear and Non-Linear Equations" from C.T. Kelley (1995)
// and "Iterative Methods for Sparse Linear Systems" from Y. Saad (2003).

// ----------------------------
// Constructor
// ----------------------------
template< typename VECTOR >
BiCGSTABsolver< VECTOR >::BiCGSTABsolver( LinearOperator< Vector > const & A,
                                          LinearOperator< Vector > const & M,
                                          real64 const tolerance,
                                          localIndex const maxIterations,
                                          integer const verbosity )
  : KrylovSolver< VECTOR >( A, M, tolerance, maxIterations, verbosity )
{

}

// ----------------------------
// Destructor
// ----------------------------
template< typename VECTOR >
BiCGSTABsolver< VECTOR >::~BiCGSTABsolver() = default;

// ----------------------------
// Monolithic BiCGSTAB solver
// ----------------------------
template< typename VECTOR >
void BiCGSTABsolver< VECTOR >::solve( Vector const & b,
                                      Vector & x ) const
{
  // Shortcuts for operators
  LinearOperator< VECTOR > const & A = m_operator;
  LinearOperator< VECTOR > const & M = m_precond;

  // Get the global size
  localIndex const N = ( m_maxIterations > 0 ) ? m_maxIterations : integer_conversion<localIndex>( x.globalSize() );

  // Get the norm of the right hand side
  real64 const normb = b.norm2();

  // Define vectors
  VectorTemp rk( x );

  // Compute initial rk
  A.residual( x, b, rk );

  // Define vectors
  VectorTemp r0_hat( rk );

  // Define scalars and initialize some
  real64 rhok_1 = rk.dot( r0_hat );
  real64 alpha = 1.0;
  real64 omegak = 1.0;

  // Define vectors and set them to 0
  VectorTemp vk( rk ); vk.zero();
  VectorTemp pk( rk ); pk.zero();
  VectorTemp y( rk );
  VectorTemp z( rk );
  VectorTemp t( rk );
  VectorTemp s( rk );
  VectorTemp h( x );
  VectorTemp q( x );

  localIndex k;
  for( k = 0; k < N; k++ )
  {
    // Compute r0_hat.rk
    real64 const rhok = rk.dot( r0_hat );

    // Compute beta
    real64 const beta = rhok / rhok_1 * alpha / omegak;

    // Keep the old value of rho
    rhok_1 = rhok;

    // Update pk = rk + beta*(pk - omega*vk)
    pk.axpy( -omegak, vk );
    pk.axpby( 1., rk, beta );

    // Uptate vk = MApk
    M.multiply( pk, y );
    A.multiply( y, vk );

    // Compute alpha
    alpha = rhok / vk.dot( r0_hat );

    // compute h = x + alpha*y
    h.copy( x );
    h.axpy( alpha, y );

    // Compute s = rk - alpha*vk
    s.copy( rk );
    s.axpy( -alpha, vk );

    // Compute z = Ms
    M.multiply( s, z );

    // Compute t = Az
    A.multiply( z, t );

    // Compute t = Mt
    M.multiply( t, q );

    // Update omega
    omegak = q.dot( z ) / q.dot( q );

    // Update x = h + omega*z
    h.axpy( omegak, z );
    x.copy( h );

    // Update rk = s - omega*t
    s.axpy( -omegak, t );
    rk = s;

    // Convergence check on ||rk||/||b||
    if( rk.norm2() / normb < 1e-8 )
    {
      break;
    }

  }

  if( m_verbosity >= 1 )
  {
    GEOSX_LOG_RANK_0( "BiCGSTAB " << (k < N ? "converged" : "did not converge") << " in " << k << " iterations." );
  }
}

// END_RST_NARRATIVE

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BiCGSTABsolver<TrilinosInterface::ParallelVector>;
template class BiCGSTABsolver<BlockVectorView<TrilinosInterface::ParallelVector>>;
#endif

#ifdef GEOSX_USE_HYPRE
template class BiCGSTABsolver<HypreInterface::ParallelVector>;
template class BiCGSTABsolver<BlockVectorView<HypreInterface::ParallelVector>>;
#endif

#ifdef GEOSX_USE_PETSC
template class BiCGSTABsolver<PetscInterface::ParallelVector>;
template class BiCGSTABsolver<BlockVectorView<PetscInterface::ParallelVector>>;
#endif

} //namespace geosx


