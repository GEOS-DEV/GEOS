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
 * @file LinearSolverParameters.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @brief Set of parameters for a linear sovler or preconditioner.
 *
 * This class holds a simple tree of linear solver options.
 * They are set to default values, but can be overwritten as needed.
 */
class LinearSolverParameters
{
public:

  integer logLevel = 0;                //!< Output level [0=none, 1=basic, 2=everything]
  string solverType = "cg";            //!< Solver type [direct, cg, gmres, bicgstab]
  string preconditionerType = "ilut";  //!< Preconditioner type [none, ilu, ilut, icc, amg]
  integer dofsPerNode = 1;             //!< Can be used to enable dense-block algorithms if available

  struct
  {
    real64 tolerance = 1e-6;      //!< Iterative solver relative residual tolerance
    integer maxIterations = 200;  //!< Maximum overall number of iterations
    integer maxRestart = 200;     //!< Maximum number of restarts (for GMRES)
    bool useAdaptiveTol = false;  //!< Whether to use an adaptive linear tolerance scheme
  }
  krylov;                         //!< Iterative (Krylov subspace) solver parameters

  struct
  {
    bool useRowScaling = false;     //!< Whether to apply row scaling to the linear system
    bool useRowColScaling = false;  //!< Whether to apply column scaling to the linear system
  }
  scaling;                          //!< Linear system scaling parameters

  struct
  {
    integer maxLevels = 20;                 //!< Maximum number of levels in the hierarchy
    string cycleType = "V";                 //!< Multigrid cycle type
    string smootherType = "gaussSeidel";    //!< Type of smoother
    string coarseType = "direct";           //!< Type of coarse solver
    integer numSweeps = 2;                  //!< Number of smoother sweeps at each level
    bool isSymmetric = true;                //!< Whether the system matrix is symmetric
    bool separateComponents = false;        //!< Whether to use separate component matrix approximation for constructing AMG hierarchy
    string nullSpaceType = "constantModes"; //!< Type of nullspace information for AMG
  }
  amg;                                      //!< Algebraic Multigrid (AMG) parameters

  struct
  {
    integer fill = 0;        //!< Level of fill for ILU(k)-type methods
    real64 threshold = 0.0;  //!< Thresholding value for ILUT-type methods
  }
  ilu;                       //!< Incomplete factorization parameters

  struct
  {
    integer overlap = 0;  //!< Overlap level between ranks
  }
  dd;                     //!< Domain decomposition parameters

  /**
   * @brief Constructor.
   */
  LinearSolverParameters() = default;

  /**
   * @brief Destructor.
   *
   */
  ~LinearSolverParameters() = default;

  /**
   * @brief Compute Einsenstat-Walker adaptive tolerance
   * @param newNewtonNorm current iteration residual norm
   * @param oldNewtonNorm previous iteration residual norm
   * @return linear solver tolerance for the next Newton iteration
   */
  static real64 eisenstatWalker( real64 const newNewtonNorm, real64 const oldNewtonNorm )
  {
    const real64 weakTol = 1e-3;
    const real64 strongTol = 1e-8;
    const real64 exponent = 2.0;
    const real64 gamma = 0.9;

    real64 normRatio = newNewtonNorm / oldNewtonNorm;
    if( normRatio > 1 ) normRatio = 1;

    real64 newKrylovTol = gamma*std::pow( normRatio, exponent );
    real64 altKrylovTol = gamma*std::pow( oldNewtonNorm, exponent );

    real64 krylovTol = std::max( newKrylovTol, altKrylovTol );
    krylovTol = std::min( krylovTol, weakTol );
    krylovTol = std::max( krylovTol, strongTol );

    return krylovTol;
  };
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
