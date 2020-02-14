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
#include "linearAlgebra/utilities/BlockMatrixView.hpp"
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
// Empty constructor
template< typename LAI >
BiCGSTABsolver< LAI >::BiCGSTABsolver()
{}

// ----------------------------
// Monolithic BiCGSTAB solver
// ----------------------------
template< typename LAI >
void BiCGSTABsolver< LAI >::solve( typename LAI::ParallelMatrix const & A,
                                   typename LAI::ParallelVector & x,
                                   typename LAI::ParallelVector const & b,
                                   typename LAI::ParallelMatrix const & M )
{

  // Get the global size
  globalIndex N = x.globalSize();

  // Placeholder for the number of iterations
  localIndex numIt = 0;

  // Get the norm of the right hand side
  real64 normb = b.norm2();

  // Define vectors
  ParallelVector rk( x );

  // Compute initial rk
  A.residual( x, b, rk );

  // Define vectors
  ParallelVector r0_hat( rk );

  // Define scalars and initialize some
  real64 rhok, rhokminus1, alpha, beta, omegak;
  rhok = 1.;
  alpha = 1.;
  omegak = 1.;

  // Define vectors and set them to 0
  ParallelVector vk( rk );
  vk.scale( 0. );
  ParallelVector pk( rk );
  pk.scale( 0. );
  ParallelVector y( rk );
  y.scale( 0. );
  ParallelVector z( rk );
  z.scale( 0. );
  ParallelVector t( rk );
  t.scale( 0. );
  ParallelVector t_( rk );

  // Declare scalar for convergence check
  real64 convCheck;

  // Declare temp scalars for alpha and beta computations
  real64 temp1;
  real64 temp2;

  for( globalIndex k = 0 ; k < N ; k++ ) // TODO: needs a maxIter param of type localIndex
  {
    // Keep the old value of rho
    rhokminus1 = rhok;

    // Compute r0_hat.rk
    rhok = rk.dot( r0_hat );

    // Compute beta
    beta = rhok / rhokminus1 * alpha / omegak;

    // Update pk = rk + beta*(pk - omega*vk)
    pk.axpy( -omegak, vk );
    pk.axpby( 1., rk, beta );

    // Uptate vk = MApk
    M.multiply( pk, y );
    A.multiply( y, vk );

    // Compute alpha
    temp1 = vk.dot( r0_hat );
    alpha = rhok / temp1;

    // compute h = x + alpha*y
    ParallelVector h( x );
    h.axpy( alpha, y );

    // Compute s = rk - alpha*vk
    ParallelVector s( rk );
    s.axpy( -alpha, vk );

    // Compute z = Ms
    M.multiply( s, z );

    // Compute t = Az
    A.multiply( z, t );

    // Compute t = Mt
    t_.copy( t );
    M.multiply( t_, t );

    // Update omega
    temp1 = t.dot( z );
    temp2 = t.dot( t );
    omegak = temp1 / temp2;

    // Update x = h + omega*z
    h.axpy( omegak, z );
    x.copy( h );

    // Update rk = s - omega*t
    s.axpy( -omegak, t );
    rk.copy( s );

    // Convergence check on ||rk||/||b||
    convCheck = rk.norm2();
    if( convCheck / normb < 1e-8 )
    {
      numIt = k;
      break;
    }

  }

  GEOSX_LOG_RANK_0( "Native BiCGSTAB (no preconditioner) converged in " << numIt << " iterations.");

}

// ----------------------------
// Block BiCGSTAB solver
// ----------------------------
template< typename LAI >
void BiCGSTABsolver<LAI>::solve( BlockMatrixView<LAI> const & GEOSX_UNUSED_PARAM( A ),
                                 BlockVectorView<LAI> & GEOSX_UNUSED_PARAM( x ),
                                 BlockVectorView<LAI> const & GEOSX_UNUSED_PARAM( b ),
                                 BlockMatrixView<LAI> const & GEOSX_UNUSED_PARAM( M ) )

{
  GEOSX_ERROR( "Not implemented" );

  // TODO: BlockVectorView is a view that doesn't handle any vector
  //       storage.  The copy and copy constructor functions below
  //       won't work.

#if 0

  // Get the global size
  globalIndex N = x.globalSize();

  // Placeholder for the number of iterations
  globalIndex numIt = 0;

  // Get the norm of the right hand side
  real64 normb = b.norm2();

  // Define vectors
  BlockVectorView<LAI> rk( x );

  // Compute initial rk
  A.residual( x, b, rk );

  // Define vectors
  BlockVectorView<LAI> r0_hat( rk );

  // Define scalars and initialize some
  real64 rhok, rhokminus1, alpha, beta, omegak;
  rhok = 1.;
  alpha = 1.;
  omegak = 1.;

  // Define vectors and set values to 0
  BlockVectorView<LAI> vk( rk );
  vk.scale( 0. );
  BlockVectorView<LAI> pk( rk );
  pk.scale( 0. );
  BlockVectorView<LAI> y( rk );
  y.scale( 0. );
  BlockVectorView<LAI> z( rk );
  z.scale( 0. );
  BlockVectorView<LAI> t( rk );
  t.scale( 0. );
  BlockVectorView<LAI> u( rk );
  t.scale( 0. );

  // Declare scalar for convergence check
  real64 convCheck;

  // Declare temp scalars for alpha and beta computations
  real64 temp1;
  real64 temp2;

  for( globalIndex k = 0 ; k < N ; k++ ) //TODO: needs a maxIter param of size localIndex
  {
    // Keep previous value of rho
    rhokminus1 = rhok;

    // Compute r0_hat.rk
    rhok = rk.dot( r0_hat );

    // Compute beta
    beta = rhok/rhokminus1*alpha/omegak;

    // Update pk = rk + beta*(pk - omega*vk)
    pk.axpy( -omegak, vk );
    pk.axpby( 1., rk, beta );

    // Uptate vk = MApk
    M.multiply( pk, y );
    A.multiply( y, vk );

    // Compute alpha
    temp1 = vk.dot( r0_hat );
    alpha = rhok/temp1;

    // compute h = x + alpha*y
    BlockVectorView<LAI> h( x );
    h.axpy( alpha, y );

    // Compute s = rk - alpha*vk
    BlockVectorView<LAI> s( rk );
    s.axpy( -alpha, vk );

    // Compute z = Ms
    M.multiply( s, z );

    // Compute t = At
    A.multiply( z, t );

    // Compute u = Mt (TODO remove u and do it in place in t)
    M.multiply( t, u );

    // Update omega
    temp1 = u.dot( z );
    temp2 = u.dot( u );
    omegak = temp1/temp2;

    // Update x = h + omega*z
    h.axpy( omegak, z );
    x.copy( h );

    // Update rk = s - omega*t
    s.axpy( -omegak, u );
    rk.copy( s );

    // Convergence check ||rk||/||b||
    convCheck = rk.norm2();
    if( convCheck/normb < 1e-8 )
    {
      numIt = k;
      break;
    }

  }

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // verbose output (TODO verbosity manager?)
  if( rank == 1 )
    std::cout << std::endl << "Block BiCGSTAB converged in " << numIt << " iterations." << std::endl;
  return;

#endif
}

// END_RST_NARRATIVE

#ifdef GEOSX_USE_TRILINOS
template class BiCGSTABsolver<TrilinosInterface>;
#endif

#ifdef GEOSX_USE_HYPRE
//template class BiCGSTABsolver<HypreInterface>;
#endif

#ifdef GEOSX_USE_PETSC
template class BiCGSTABsolver<PetscInterface>;
#endif

} //namespace geosx
