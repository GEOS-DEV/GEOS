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

/*
 * BiCGSTABsolver.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Matthias
 */

#ifndef SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BICGSTABSOLVER_HPP_
#define SRC_CORECOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BICGSTABSOLVER_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"

namespace geosx
{

/**
 * \class BlockLinearSolvers
 * \brief This class creates and provides basic support for block
 *        linear solvers (templated on the LA interface).
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
   * @brief Empty matrix constructor.
   *
   * Create an empty block matrix.
   */
  BiCGSTABsolver();

  /**
   * @brief Virtual destructor.
   */
  virtual ~BiCGSTABsolver() = default;
  //@}

  void solve( ParallelMatrix const &A,
              ParallelVector &x,
              ParallelVector const &b,
              ParallelMatrix const &M );

  void solve( BlockMatrixView<LAI> const &A,
              BlockVectorView<LAI> &x,
              BlockVectorView<LAI> const &b,
              BlockMatrixView<LAI> const &M );

private:

};

// Empty constructor
template< typename LAI >
BiCGSTABsolver<LAI>::BiCGSTABsolver()
{}

template< typename LAI >
void BiCGSTABsolver<LAI>::solve( typename LAI::ParallelMatrix const &A,
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

  // Compute initial rk
  A.residual( x, b, rk );

  // Define vectors
  ParallelVector r0_hat( rk );

  real64 rhok, rhokminus1, alpha, beta, omegak;
  rhok = 1.;
  alpha = 1.;
  omegak = 1.;

  // Define vectors
  ParallelVector vk( rk );
  vk.scale(0.);
  ParallelVector pk( rk );
  pk.scale(0.);
  ParallelVector y( rk );
  y.scale(0.);
  ParallelVector z( rk );
  z.scale(0.);
  ParallelVector t( rk );
  t.scale(0.);
  ParallelVector u( rk );
  u.scale(0.);

  // Declare scalar for convergence check
  real64 convCheck;

  for( typename LAI::laiGID k = 0 ; k < N ; k++ )
  {
    rhokminus1 = rhok;

    // Compute r0_hat.rk
    rk.dot( r0_hat, &rhok );

    // Compute beta
    beta = rhok/rhokminus1*alpha/omegak;

    // Update pk
    pk.update(-omegak,vk,1.);
    pk.update(1.,rk,beta);

    // Uptate vk
    M.multiply(pk, y);
    A.multiply(y, vk);

    // Compute alpha
    real64 temp1;
    vk.dot( r0_hat, &temp1 );
    alpha = rhok/temp1;

    // compute h = x + alpha*y
    ParallelVector h( x );
    h.update(alpha, y, 1.);

    // Compute s = rk -alpha*vk
    ParallelVector s( rk );
    s.update(-alpha, vk, 1.);

    // Compute z
    M.multiply(s, z);

    // Compute t
    A.multiply(z, t);

    // Compute u
    M.multiply(t, u);

    // Update omega
    real64 temp2;
    u.dot( z, &temp1 );
    u.dot( u, &temp2 );
    omegak = temp1/temp2;

    h.update(omegak, z, 1.);

    x.update(1., h, 0.);

    s.update(-omegak, t, 1.);

    rk.update(1., s, 0.);

    // Convergence check
    rk.norm2( convCheck );

    if( convCheck < 1e-8 )
    {
      numIt = k;
      break;
    }

    //std::cout << k << ", " << convCheck << std::endl;

  }

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if ( rank == 1 )
    std::cout << "BiCGSTAB converged in " << numIt << " iterations." << std::endl;
  return;

}

template< typename LAI >
void BiCGSTABsolver<LAI>::solve( BlockMatrixView<LAI> const &A,
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

  // Compute initial rk
  A.residual( x, b, rk );

  // Preconditioning
  BlockVectorView<LAI> zk( x );
  M.multiply( rk, zk );

  // pk = zk
  BlockVectorView<LAI> pk( zk );
  BlockVectorView<LAI> Apk( zk );

  real64 alpha, beta;

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

    // Update x
    x.update( alpha, pk, 1.0 );

    // Update rk
    BlockVectorView<LAI> rkold( rk );
    BlockVectorView<LAI> zkold( zk );
    rk.update( -alpha, Apk, 1.0 );

    // Convergence check
    rk.norm2( convCheck );

    if( convCheck < 1e-8 )
    {
      numIt = k;
      break;
    }

    M.multiply( rk, zk );

    zk.dot( rk, beta );
    zkold.dot( rkold, temp );
    beta = beta/temp;

    // Update pk
    pk.update( 1.0, zk, beta );

    //std::cout << k << ", " << convCheck << std::endl;

  }

  // Get the MPI rank
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  if ( rank == 1 )
    std::cout << "CG converged in " << numIt << " iterations." << std::endl;
  return;

}

} // namespace GEOSX

#endif /* SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BiCGSTABsolver_HPP_ */
