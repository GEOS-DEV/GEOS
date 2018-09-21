/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file CGsolver.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Matthias Cremon
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_CGSOLVER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_CGSOLVER_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"

namespace geosx
{

/**
 * \class CGsolver
 * \brief This class creates and provides basic support for block
 *        CG (templated on the LA interface).
 * \note  The notation is consistent with "Iterative Methods for
 *        Linear and Non-Linear Equations" from C.T. Kelley (1995)
 *        and "Iterative Methods for Sparse Linear Systems"
 *        from Y. Saad (2003).
 */

template< typename LAI >
class CGsolver
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty constructor.
   *
   * Creates a solver object.
   */
  CGsolver();

  /**
   * @brief Virtual destructor.
   */
  virtual ~CGsolver() = default;
  //@}

  /**
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with CG
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
   * @brief Solve the system <tt>M^{-1}(Ax - b) = 0</tt> with CG
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

// Empty constructor
template< typename LAI >
CGsolver<LAI>::CGsolver()
{}

template< typename LAI >
void CGsolver<LAI>::solve( typename LAI::ParallelMatrix const &A,
                           typename LAI::ParallelVector &x,
                           typename LAI::ParallelVector const &b,
                           typename LAI::ParallelMatrix const &M )

{

  // Get the global size
  typename LAI::laiGID N = x.globalSize();

  // Placeholder for the number of iterations
  typename LAI::laiGID numIt = 0;

  // Define vectors
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
  real64 convCheck;
  rk.norm2( convCheck );


  for( typename LAI::laiGID k = 0 ; k < N ; k++ )
  {
    // Compute rkT.rk
    rk.dot( zk, alpha );

    // Compute Apk
    A.multiply( pk, Apk );

    // compute alpha
    real64 temp;
    pk.dot( Apk, temp );
    alpha = alpha/temp;

    // Update x = x + alpha*ph
    x.update( alpha, pk, 1.0 );

    // Update rk = rk - alpha*Apk
    ParallelVector rkold( rk );
    ParallelVector zkold( zk );
    rk.update( -alpha, Apk, 1.0 );

    // Convergence check on ||rk||_2
    rk.norm2( convCheck );

    if( convCheck < 1e-8 )
    {
      numIt = k;
      break;
    }

    // Update zk = Mrk
    M.multiply( rk, zk );

    // Compute beta
    zk.dot( rk, beta );
    zkold.dot( rkold, temp );
    beta = beta/temp;

    // Update pk = pk + beta*zk
    pk.update( 1.0, zk, beta );

    //std::cout << k << ", " << convCheck << std::endl;

  }

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // veborse output (TODO verbosity manager?)
  if( rank == 1 )
    std::cout << "CG converged in " << numIt << " iterations." << std::endl;
  return;

}

template< typename LAI >
void CGsolver<LAI>::solve( BlockMatrixView<LAI> const &A,
                           BlockVectorView<LAI> &x,
                           BlockVectorView<LAI> const &b,
                           BlockMatrixView<LAI> const &M )

{

  // Get the global size
  typename LAI::laiGID N = x.globalSize();

  // Placeholder for the number of iterations
  typename LAI::laiGID numIt = 0;

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
  real64 convCheck;
  rk.norm2( convCheck );

  for( typename LAI::laiGID k = 0 ; k < N ; k++ )
  {
    // Compute rkT.rk
    rk.dot( zk, alpha );

    // Compute Apk
    A.multiply( pk, Apk );

    // compute alpha
    real64 temp;
    pk.dot( Apk, temp );

    alpha = alpha/temp;

    // Update x = x + alpha*pk
    x.update( alpha, pk, 1.0 );

    // Update rk = rk - alpha*Apk
    BlockVectorView<LAI> rkold( rk );
    BlockVectorView<LAI> zkold( zk );
    rk.update( -alpha, Apk, 1.0 );

    // Convergence check on ||rk||_2
    rk.norm2( convCheck );

    if( convCheck < 1e-8 )
    {
      numIt = k;
      break;
    }

    // Update zk = Mrk
    M.multiply( rk, zk );

    // Compute beta
    zk.dot( rk, beta );
    zkold.dot( rkold, temp );
    beta = beta/temp;

    // Update pk = pk + beta*zk
    pk.update( 1.0, zk, beta );

    //std::cout << k << ", " << convCheck << std::endl;

  }

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  // veborse output (TODO verbosity manager?)
  if( rank == 1 )
    std::cout << "Block CG converged in " << numIt << " iterations." << std::endl;
  return;

}

} // namespace GEOSX

#endif /* SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_CGSOLVER_HPP_ */
