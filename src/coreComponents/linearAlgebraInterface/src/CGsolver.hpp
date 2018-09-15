/*
 * BlockLinearSolvers.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Matthias
 */

#ifndef SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_CGSOLVER_HPP_
#define SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_CGSOLVER_HPP_

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
class CGsolver
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
  CGsolver();

  /**
   * @brief Virtual destructor.
   */
  virtual ~CGsolver() = default;
  //@}

  void solve( ParallelMatrix const &A,
              ParallelVector &x,
              ParallelVector const &b,
              ParallelMatrix const &M );

  void solve( BlockMatrixView<LAI> &A,
              BlockVectorView<LAI> &x,
              BlockVectorView<LAI> &b,
              BlockMatrixView<LAI> &M );

private:

};

// Empty constructor (inlined)
template< typename LAI >
inline
CGsolver<LAI>::CGsolver()
{}

template< typename LAI >
void CGsolver<LAI>::solve( typename LAI::ParallelMatrix const &A,
                           typename LAI::ParallelVector &x,
                           typename LAI::ParallelVector const &b,
                           typename LAI::ParallelMatrix const &M )

{

  // Get the global size
  TrilinosInterface::laiGID N = x.globalSize();

  // Placeholder for the number of iterations
  TrilinosInterface::laiGID numIt = 0;

  // Define vectors
  ParallelVector rk( x );

  // Compute initial rk
  A.residual( x, b, rk );
  rk.scale( -1.0 );

  // Preconditioning
  ParallelVector zk( x );
  M.multiply( rk, zk );

  // pk = zk
  ParallelVector pk( zk );
  ParallelVector Apk( zk );

  real64 alpha, beta;

  real64 convCheck;
  rk.norm2( convCheck );

  for( TrilinosInterface::laiGID k = 0 ; k < N ; k++ )
  {
    // Compute rkT.rk
    rk.dot( zk, &alpha );

    // Compute Apk
    A.multiply( pk, Apk );

    // compute alpha
    real64 temp;
    pk.dot( Apk, &temp );
    alpha = alpha/temp;

    // Update x
    x.update( alpha, pk, 1.0 );

    // Update rk
    ParallelVector rkold( rk );
    ParallelVector zkold( zk );
    rk.update( -alpha, Apk, 1.0 );

    // Convergence check
    rk.norm2( convCheck );

    if( convCheck < 1e-8 )
    {
      numIt = k;
      break;
    }

    M.multiply( rk, zk );

    zk.dot( rk, &beta );
    zkold.dot( rkold, &temp );
    beta = beta/temp;

    // Update pk
    pk.update( 1.0, zk, beta );

    //std::cout << k << ", " << convCheck << std::endl;

  }

  std::cout << "CG converged in " << numIt << " iterations." << std::endl;
  return;

}

template< typename LAI >
void CGsolver<LAI>::solve( BlockMatrixView<LAI> &A,
                           BlockVectorView<LAI> &x,
                           BlockVectorView<LAI> &b,
                           BlockMatrixView<LAI> &M )

{

  // Get the global size
  TrilinosInterface::laiGID N = x.globalSize();

  // Placeholder for the number of iterations
  TrilinosInterface::laiGID numIt = 0;

  // Define vectors
  BlockVectorView<LAI> rk( x );

  // Compute initial rk
  A.residual( x, b, rk );
  rk.scale( -1.0 );

  // Preconditioning
  BlockVectorView<LAI> zk( x );
  M.multiply( rk, zk );

  // pk = zk
  BlockVectorView<LAI> pk( zk );
  BlockVectorView<LAI> Apk( zk );

  real64 alpha, beta;

  real64 convCheck;
  rk.norm2( convCheck );

  for( TrilinosInterface::laiGID k = 0 ; k < N ; k++ )
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

  std::cout << "CG converged in " << numIt << " iterations." << std::endl;
  return;

}

} // namespace GEOSX

#endif /* SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_CGSOLVER_HPP_ */
