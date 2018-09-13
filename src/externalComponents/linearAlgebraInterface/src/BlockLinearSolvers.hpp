/*
 * BlockLinearSolvers.hpp
 *
 *  Created on: Sep 12, 2018
 *      Author: Matthias
 */

#ifndef SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKLINEARSOLVERS_HPP_
#define SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKLINEARSOLVERS_HPP_

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
class BlockLinearSolvers
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

public:

//void CG( ParallelMatrix const &A,
//         ParallelVector &x,
//         ParallelVector const &b,
//         ParallelMatrix const &M);

private:

};

template< typename LAI >
void CG( typename LAI::ParallelMatrix const &A,
         typename LAI::ParallelVector &x,
         typename LAI::ParallelVector const &b,
         typename LAI::ParallelMatrix const &M )

{

  using ParallelVector = typename LAI::ParallelVector;

  // Get the global size
  globalIndex N = x.globalSize();

  // Placeholder for the number of iterations
  globalIndex numIt = 0;

  // Define vectors
  ParallelVector rk(x);

  // Compute initial rk
  A.residual(x,b,rk);
  rk.scale(-1.0);

  // Preconditioning
  ParallelVector zk(x);
  M.multiply(rk,zk);

  // pk = zk
  ParallelVector pk(zk);
  ParallelVector Apk(zk);

  real64 alpha, beta;

  real64 convCheck;
  rk.norm2(convCheck);

  for (globalIndex k = 0; k < N; k++)
  {
    // Compute rkT.rk
    rk.dot(zk,&alpha);

    // Compute Apk
    A.multiply(pk,Apk);

    // compute alpha
    real64 temp;
    pk.dot(Apk,&temp);
    alpha = alpha/temp;

    // Update x
    x.update(alpha,pk,1.0);

    // Update rk
    ParallelVector rkold(rk);
    ParallelVector zkold(zk);
    rk.update(-alpha,Apk,1.0);

    // Convergence check
    rk.norm2(convCheck);

    if (convCheck < 1e-6)
    {
      numIt = k;
      break;
    }

    M.multiply(rk,zk);

    zk.dot(rk,&beta);
    zkold.dot(rkold,&temp);
    beta = beta/temp;

    // Update pk
    pk.update(1.0,zk,beta);

  }

  std::cout << "CG converged in " << numIt << " iterations." << std::endl;
  return;

}

} // namespace GEOSX

#endif /* SRC_EXTERNALCOMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKLINEARSOLVERS_HPP_ */
