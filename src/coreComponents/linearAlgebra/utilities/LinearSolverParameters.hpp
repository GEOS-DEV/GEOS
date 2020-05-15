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
 * @brief Set of parameters for a linear solver or preconditioner.
 *
 * This class holds a simple tree of linear solver options.
 * They are set to default values, but can be overwritten as needed.
 */
class LinearSolverParameters
{
public:

  integer logLevel = 0;                ///< Output level [0=none, 1=basic, 2=everything]
  string solverType = "direct";        ///< Solver type [direct, cg, gmres, bicgstab, preconditioner]
  string preconditionerType = "iluk";  ///< Preconditioner type [none, iluk, ilut, amg, mgr, block]
  integer dofsPerNode = 1;             ///< Dofs per node (or support location) for non-scalar problems

  /// Krylov-method parameters
  struct Krylov
  {
    real64 relTolerance = 1e-6;       ///< Relative convergence tolerance for iterative solvers
    integer maxIterations = 200;      ///< Max iterations before declaring convergence failure
    integer maxRestart = 200;         ///< Max number of vectors in Krylov basis before restarting
    integer useAdaptiveTol = false;   ///< Use Eisenstat-Walker adaptive tolerance
    real64 weakestTol = 1e-3;         ///< Weakest allowed tolerance when using adaptive method
  }
  krylov;                             ///< Krylov-method parameter struct

  /// Matrix-scaling parameters
  struct Scaling
  {
    integer useRowScaling = false;      ///< Apply row scaling
    integer useRowColScaling = false;   ///< Apply row and column scaling (not yet implemented)
  }
  scaling;                              ///< Matrix-scaling parameter struct

  /// Algebraic multigrid parameters
  struct AMG
  {
    integer maxLevels = 20;                  ///< Maximum number of coarsening levels
    string cycleType = "V";                  ///< AMG cycle type
    string smootherType = "gaussSeidel";     ///< Smoother type
    string coarseType = "direct";            ///< Coarse-level solver/smoother
    integer numSweeps = 2;                   ///< Number of smoother sweeps
    string preOrPostSmoothing = "both";      ///< Pre and/or post smoothing [pre,post,both]
    real64 threshold = 0.0;                  ///< Threshold for "strong connections" (for classical and
                                             ///< smoothed-aggregation AMG)
    integer isSymmetric = true;              ///< Identify if matrix is symmetric
    integer separateComponents = false;      ///< Apply a separate component filter before AMG construction
    string nullSpaceType = "constantModes";  ///< Null space type [constantModes,rigidBodyModes]
  }
  amg;                                      //!< Algebraic Multigrid (AMG) parameters

  /// Incomplete factorization parameters
  struct ILU
  {
    integer fill = 0;        ///< Fill level
    real64 threshold = 0.0;  ///< Dropping threshold
  }
  ilu;                       ///< Incomplete factorization parameter struct

  /// Domain decomposition parameters
  struct DD
  {
    integer overlap = 0;   ///< Ghost overlap
  }
  dd;                      ///< Domain decomposition parameter struct

  /// Constructor.
  LinearSolverParameters() = default;

  /// Destructor.
  ~LinearSolverParameters() = default;

  /**
   * @brief Eisenstat-Walker adaptive tolerance
   *
   * This method enables an inexact-Newton method is which the linear solver
   * tolerance is chosen based on the nonlinear solver convergence behavior.
   * In early Newton iterations, the search direction is usually imprecise, and
   * therefore a weak linear convergence tolerance can be chosen to minimize
   * computational cost.  As the search gets closer to the true solution, however,
   * more stringent linear tolerances are necessary to maintain quadratic convergence
   * behavior.
   *
   * The user can set the weakest tolerance allowed, with a default of 1e-3.
   * Even weaker values (e.g. 1e-2,1e-1) can be used for further speedup, but may
   * occasionally cause convergence problems.  Use this parameter with caution.  The
   * most stringent tolerance is hardcoded to 1e-8, which is sufficient for
   * most problems.
   *
   * See Eisenstat, S.C. and Walker, H.F., 1996. Choosing the forcing terms in an
   * inexact Newton method. SIAM Journal on Scientific Computing, 17(1), pp.16-32.
   *
   * @param newNewtonNorm Residual norm at current iteration
   * @param oldNewtonNorm Residual norm at previous iteration
   * @param weakestTol Weakest tolerance allowed (default 1e-3).
   * @return Adaptive tolerance recommendation
   */
  static real64 eisenstatWalker( real64 const newNewtonNorm,
                                 real64 const oldNewtonNorm,
                                 real64 const weakestTol )
  {
    real64 const strongestTol = 1e-8;
    real64 const exponent = 2.0;
    real64 const gamma = 0.9;

    real64 normRatio = newNewtonNorm / oldNewtonNorm;
    if( normRatio > 1 ) normRatio = 1;

    real64 newKrylovTol = gamma*std::pow( normRatio, exponent );
    real64 altKrylovTol = gamma*std::pow( oldNewtonNorm, exponent );

    real64 krylovTol = std::max( newKrylovTol, altKrylovTol );
    krylovTol = std::min( krylovTol, weakestTol );
    krylovTol = std::max( krylovTol, strongestTol );

    return krylovTol;
  };
};

/**
 * @brief Linear solver parameters with Group capabilities
 *
 * This class is a derived version of LinearSolverParameters with
 * dataRepository::Group capabilities to allow for XML input.
 *
 * Developer note: This class exposes LAI solver settings to external users.
 * As a general philosophy, only a subset of frequently tuned parameters should
 * be exposed.  Many advanced parameters can be set by the PhysicsSolver itself
 * since it has knowledge of the underlying problem (e.g. isSymmetric = true,
 * dofsPerNode = 3, etc.).  While we want to enable power users to tune the
 * solvers, most users prefer a short list of options with good default settings
 */
class LinearSolverParametersGroup : public LinearSolverParameters, public dataRepository::Group
{
public:

  LinearSolverParametersGroup() = delete;

  /// Constructor
  LinearSolverParametersGroup( std::string const & name, Group * const parent );

  /// Copy constructor
  LinearSolverParametersGroup( LinearSolverParametersGroup && ) = default;

  /// Destructor
  virtual ~LinearSolverParametersGroup() override = default;

  /// Catalog name
  static string CatalogName() { return "LinearSolverParameters"; }

  /// Postprocessing of input
  virtual void PostProcessInput() override;

  /// Keys appearing in XML
  struct viewKeysStruct
  {
    static constexpr auto solverTypeString         = "solverType";         ///< Solver type key
    static constexpr auto preconditionerTypeString = "preconditionerType"; ///< Preconditioner type key

    static constexpr auto krylovMaxIterString     = "krylovMaxIter";     ///< Krylov max iterations key
    static constexpr auto krylovTolString         = "krylovTol";         ///< Krylov tolerance key
    static constexpr auto krylovAdaptiveTolString = "krylovAdaptiveTol"; ///< Krylov adaptive tolerance key
    static constexpr auto krylovWeakTolString     = "krylovWeakestTol";  ///< Krylov weakest tolerance key

    static constexpr auto amgNumSweepsString = "amgNumSweeps";             ///< AMG number of sweeps key
    static constexpr auto amgSmootherString  = "amgSmootherType";          ///< AMG smoother type key
    static constexpr auto amgCoarseString    = "amgCoarseSolver";          ///< AMG coarse solver key
    static constexpr auto amgThresholdString = "amgThreshold";             ///< AMG threshold key

    static constexpr auto iluFillString      = "iluFill";       ///< ILU fill key
    static constexpr auto iluThresholdString = "iluThreshold";  ///< ILU threshold key
  } viewKeys;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
