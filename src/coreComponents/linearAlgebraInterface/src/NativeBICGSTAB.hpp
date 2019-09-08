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

/**
 * @file BiCGSTABsolver.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Matthias Cremon
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_SOLVER_BICGSTAB_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_SOLVER_BICGSTAB_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"
#include "BlockMatrixView.hpp"
#include "BlockVectorView.hpp"

namespace geosx
{

/**
 * \class BiCGSTABsolver
 * \brief This class creates and provides basic support for block
 *        BiCGSTAB (templated on the LA interface).
 * \note  The notation is consistent with "Iterative Methods for
 *        Linear and Non-Linear Equations" from C.T. Kelley (1995)
 *        and "Iterative Methods for Sparse Linear Systems"
 *        from Y. Saad (2003).
 */

template< typename LAI >
class BiCGSTABsolver
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty solver object constructor.
   *
   * Create an empty solver object.
   */
  BiCGSTABsolver();

  /**
   * @brief Virtual destructor.
   */
  virtual ~BiCGSTABsolver() = default;
  //@}

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with BiCGSTAB
   * using monolithic GEOSX matrices.
   *
   * \param A system matrix.
   * \param x system solution (input = initial guess, output = solution).
   * \param b system right hand side.
   * \param M preconditioner.
   */
  void solve( ParallelMatrix const &A,
              ParallelVector &x,
              ParallelVector const &b,
              ParallelMatrix const &M );

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with BiCGSTAB
   * using block GEOSX matrices.
   *
   * \param A system block matrix.
   * \param x system block solution (input = initial guess, output = solution).
   * \param b system block right hand side.
   * \param M block preconditioner.
   */
  void solve( BlockMatrixView<LAI> const &A,
              BlockVectorView<LAI> &x,
              BlockVectorView<LAI> const &b,
              BlockMatrixView<LAI> const &M );

private:

};

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
BiCGSTABsolver<LAI>::BiCGSTABsolver()
{}

// ----------------------------
// Monolithic BiCGSTAB solver
// ----------------------------
template< typename LAI >
void BiCGSTABsolver<LAI>::solve( typename LAI::ParallelMatrix const &A,
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
    M.multiply( t, t );

    // Update omega
    temp1 = t.dot( z );
    temp2 = t.dot( t );
    omegak = temp1/temp2;

    // Update x = h + omega*z
    h.axpy( omegak, z );
    x.copy( h );

    // Update rk = s - omega*t
    s.axpy( -omegak, t );
    rk.copy( s );

    // Convergence check on ||rk||/||b||
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

  // TODO verbosity management
  if( rank == 0 )
    std::cout << "Native BiCGSTAB (no preconditioner) converged in " << numIt << " iterations." << std::endl;
  return;

}

// ----------------------------
// Block BiCGSTAB solver
// ----------------------------
template< typename LAI >
void BiCGSTABsolver<LAI>::solve( BlockMatrixView<LAI> const &A,
                                 BlockVectorView<LAI> &x,
                                 BlockVectorView<LAI> const &b,
                                 BlockMatrixView<LAI> const &M )

{
  GEOS_ERROR( "Not implemented" );

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

} // namespace GEOSX

#endif /* SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BiCGSTABsolver_HPP_ */
