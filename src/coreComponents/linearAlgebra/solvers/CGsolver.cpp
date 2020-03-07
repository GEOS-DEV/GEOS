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
#include "linearAlgebra/utilities/BlockMatrixView.hpp"
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
template< typename LAI >
CGsolver<LAI>::CGsolver()
{}

// ----------------------------
// Monolithic CG solver
// ----------------------------
template< typename LAI >
void CGsolver<LAI>::solve( typename LAI::ParallelMatrix const &A,
                           typename LAI::ParallelVector &x,
                           typename LAI::ParallelVector const &b,
                           typename LAI::ParallelMatrix const &M )

{

  // Get the global size
  globalIndex N = x.globalSize();

  // Placeholder for the number of iterations
  localIndex numIt = 0;

  // Get the norm of the right hand side
  real64 normb = b.norm2();

  // Define residual vector
  ParallelVector rk( x );

  // Compute initial rk =  b - Ax
  A.residual( x, b, rk );

  // Preconditioning
  ParallelVector zk( x );
  M.multiply( rk, zk );

  // pk = zk
  ParallelVector pk( zk );
  ParallelVector Apk( zk );

  // Declare alpha and beta scalars
  real64 alpha, beta;

  // Convergence check
  real64 convCheck = rk.norm2();

  // Declare temp scalar for alpha computation.
  real64 temp;

  // Declare older vectors
  ParallelVector rkold( rk );
  ParallelVector zkold( zk );

  for( globalIndex k = 0 ; k < N ; k++ ) // TODO: this needs a max_iter param of type localIndex
  {
    // Compute rkT.rk
    alpha = rk.dot( zk );

    // Compute Apk
    A.multiply( pk, Apk );

    // compute alpha
    temp = pk.dot( Apk );
    alpha = alpha/temp;

    // Update x = x + alpha*ph
    x.axpby( alpha, pk, 1.0 );

    // Update rk = rk - alpha*Apk
    rkold.copy( rk );
    zkold.copy( zk );
    rk.axpby( -alpha, Apk, 1.0 );

    // Convergence check on ||rk||/||b||
    convCheck = rk.norm2();
    if( convCheck/normb < 1e-8 )
    {
      numIt = k;
      break;
    }

    // Update zk = Mrk
    M.multiply( rk, zk );

    // Compute beta
    beta = zk.dot( rk );
    temp = zkold.dot( rkold );
    beta = beta/temp;

    // Update pk = pk + beta*zk
    pk.axpby( 1.0, zk, beta );

  }

  GEOSX_LOG_RANK_0( "Native CG (no preconditioner) converged in " << numIt << " iterations.");
  return;

}

// ----------------------------
// Block CG solver
// ----------------------------
template< typename LAI >
void CGsolver<LAI>::solve( BlockMatrixView<LAI> const & GEOSX_UNUSED_PARAM( A ),
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
  localIndex numIt = 0;

  // Get the norm of the right hand side
  real64 normb = b.norm2();

  // Define vectors
  BlockVectorView<LAI> rk( x );

  // Compute initial rk =  b - Ax
  A.residual( x, b, rk );

  // Preconditioning
  BlockVectorView<LAI> zk( x );
  M.multiply( rk, zk );

  // pk = zk
  BlockVectorView<LAI> pk( zk );
  BlockVectorView<LAI> Apk( zk );

  // Declare alpha and beta scalars
  real64 alpha, beta;

  // Convergence check
  real64 convCheck = rk.norm2();

  // Declare temp scalar for alpha computation
  real64 temp;

  // Declare older vectors
  BlockVectorView<LAI> rkold( rk );
  BlockVectorView<LAI> zkold( zk );

  for( globalIndex k = 0 ; k < N ; k++ ) // TODO: needs maxIter param of type localIndex
  {
    // Compute rkT.rk
    alpha = rk.dot( zk );

    // Compute Apk
    A.multiply( pk, Apk );

    // compute alpha
    temp = pk.dot( Apk );

    alpha = alpha/temp;

    // Update x = x + alpha*pk
    x.axpby( alpha, pk, 1.0 );

    // Update rk = rk - alpha*Apk
    rkold.copy( rk );
    zkold.copy( zk );
    rk.axpby( -alpha, Apk, 1.0 );

    // Convergence check on ||rk||/||b||
    convCheck = rk.norm2();
    if( convCheck/normb < 1e-8 )
    {
      numIt = k;
      break;
    }

    // Update zk = Mrk
    M.multiply( rk, zk );

    // Compute beta
    beta = zk.dot( rk );
    temp = zkold.dot( rkold );
    beta = beta/temp;

    // Update pk = pk + beta*zk
    pk.axpby( 1.0, zk, beta );

  }

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // verbose output (TODO verbosity manager?)
  if( rank == 1 )
    std::cout << std::endl << "Block CG converged in " << numIt << " iterations." << std::endl;
  return;

#endif
}

// END_RST_NARRATIVE

#ifdef GEOSX_USE_TRILINOS
template class CGsolver<TrilinosInterface>;
#endif

#ifdef GEOSX_USE_HYPRE
//template class CGsolver<HypreInterface>;
#endif

#ifdef GEOSX_USE_PETSC
template class CGsolver<PetscInterface>;
#endif

} //namespace geosx
