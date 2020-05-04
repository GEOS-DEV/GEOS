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

#include "dataRepository/Group.hpp"

namespace geosx
{

/**
 * @class LinearSolverParameters
 * This class holds a simple tree of linear solver options.  They are set to
 * default values, but can be overwritten as needed.
 */

class LinearSolverParameters
{
public:

  integer logLevel = 0;                //!< Output level [0=none, 1=basic, 2=everything]
  string solverType = "direct";        //!< Solver type [direct, cg, gmres, bicgstab, preconditioner]
  string preconditionerType = "ilut";  //!< Preconditioner type [none, iluk, ilut, icc, amg]
  integer dofsPerNode = 1;             //!< Can be used to enable dense-block algorithms if available

  struct
  {
    real64 tolerance = 1e-6;
    integer maxIterations = 200;
    integer maxRestart = 200;
    integer useAdaptiveTol = false;
  }
  krylov;

  struct
  {
    integer useRowScaling = false;
    integer useRowColScaling = false; // not currently used
  }
  scaling;

  struct
  {
    integer maxLevels = 20;
    string cycleType = "V";
    string smootherType = "gaussSeidel";
    string coarseType = "direct";
    integer numSweeps = 2;
    real64 aggregationThreshold = 0.0;
    integer isSymmetric = true;
    integer separateComponents = false;
    string nullSpaceType = "constantModes";
  }
  amg;

  struct
  {
    integer fill = 0;
    real64 threshold = 0.0;
  }
  ilu;

  struct
  {
    integer overlap = 0;
  }
  dd;

  /**
   * @brief Constructor.
   */
  LinearSolverParameters() = default;

  /**
   * @brief Destructor.
   */
  ~LinearSolverParameters() = default;

  /**
   * @brief Einsenstat-Walker adaptive tolerance
   */
  static real64 eisenstatWalker( real64 newNewtonNorm, real64 oldNewtonNorm )
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

/**
 * @class LinearSolverParametersGroup
 * This class is a derived version of LinearSolverParameters with
 * dataRepository::Group capabilities (to allow for XML input)
 */
class LinearSolverParametersGroup : public LinearSolverParameters, public dataRepository::Group
{
public:

  LinearSolverParametersGroup() = delete;

  LinearSolverParametersGroup( std::string const & name, Group * const parent );

  LinearSolverParametersGroup( LinearSolverParametersGroup && ) = default;

  virtual ~LinearSolverParametersGroup() override = default;

  static string CatalogName() { return "LinearSolverParameters"; }

  virtual void PostProcessInput() override;

  // note: Only a subset of frequently used parameters should be exposed to users.
  //       Many advanced parameters can be set by the physicSolver itself
  //       (e.g. a solver will already know how many dofsPerNode it has).
  //       In typical usage the user really should just choose between
  //       direct and krylov, and set the solver tolerance for the latter.
  //       The physicsSolver be set up with ``optimal`` strategies by default.

  struct viewKeysStruct
  {
    static constexpr auto solverTypeString         = "solverType";
    static constexpr auto preconditionerTypeString = "preconditionerType";

    static constexpr auto krylovTolString         = "krylovTol";
    static constexpr auto krylovAdaptiveTolString = "krylovAdaptiveTol";
    static constexpr auto krylovMaxIterString     = "krylovMaxIter";

    static constexpr auto amgNumSweepsString   = "amgNumSweeps";
    static constexpr auto amgSmootherString    = "amgSmootherType";
    static constexpr auto amgCoarseString      = "amgCoarseSolver";
    static constexpr auto amgAggregationString = "amgAggregationThreshold";

    static constexpr auto iluFillString      = "iluFill";
    static constexpr auto iluThresholdString = "iluThreshold";
  } viewKeys;

  //struct groupKeysStruct
  //{} groupKeys;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
