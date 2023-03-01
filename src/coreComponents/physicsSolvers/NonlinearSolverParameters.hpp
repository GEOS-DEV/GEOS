/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_
#define GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "physicsSolvers/SolverBaseKernels.hpp"

namespace geosx
{

/**
 * @class Holds parameters and status for execution of nonlinear solution schemes.
 */
class NonlinearSolverParameters : public dataRepository::Group
{
public:

  /**
   * @brief Default constructor.
   */
  NonlinearSolverParameters() = delete;

  /**
   * @brief Constructor
   * @param[in] name The name of the new instantiation of this Group.
   * @param[in] parent A pointer to the parent of this Group.
   */
  NonlinearSolverParameters( string const & name,
                             Group * const parent );

  /**
   * @brief Default destructor
   */
  virtual ~NonlinearSolverParameters() = default;

  /**
   * @brief Default Move Constructor
   * @param The source object of the move.
   */
  NonlinearSolverParameters( NonlinearSolverParameters && ) = default;

  /**
   * @brief The name of this object in the catalog.
   * @return A string containing the name of this object in the catalog.
   * The CatalogName is the string that will result in the creation of a new
   * NonlinearSolverParameters object when calling
   * Group::getCatalog()::Allocate().
   */
  static string catalogName() { return "NonlinearSolverParameters"; }

  virtual void postProcessInput() override;

  struct viewKeysStruct
  {
    static constexpr char const * lineSearchActionString()        { return "lineSearchAction"; }
    static constexpr char const * lineSearchMaxCutsString()       { return "lineSearchMaxCuts"; }
    static constexpr char const * lineSearchCutFactorString()     { return "lineSearchCutFactor"; }
    static constexpr char const * lineSearchInterpolationTypeString() { return "lineSearchInterpolationType"; }

    static constexpr char const * normTypeString()                { return "normType"; }
    static constexpr char const * newtonTolString()               { return "newtonTol"; }
    static constexpr char const * newtonMaxIterString()           { return "newtonMaxIter"; }
    static constexpr char const * newtonMinIterString()           { return "newtonMinIter"; }
    static constexpr char const * newtonNumIterationsString()     { return "newtonNumberOfIterations"; }
    static constexpr char const * newtonSplitOperMaxIterString()  { return "newtonSplitOperMaxIter"; }

    static constexpr char const * allowNonConvergedString()       { return "allowNonConverged"; }
    static constexpr char const * timeStepDecreaseIterLimString() { return "timeStepDecreaseIterLimit"; }
    static constexpr char const * timeStepIncreaseIterLimString() { return "timeStepIncreaseIterLimit"; }
    static constexpr char const * timeStepDecreaseFactorString()  { return "timeStepDecreaseFactor"; }
    static constexpr char const * timeStepIncreaseFactorString()  { return "timeStepIncreaseFactor"; }

    static constexpr char const * maxSubStepsString()             { return "maxSubSteps"; }
    static constexpr char const * maxTimeStepCutsString()         { return "maxTimeStepCuts"; }
    static constexpr char const * minNumNewtonIterationsString()  { return "minNumberOfNewtonIterations"; }
    static constexpr char const * timeStepCutFactorString()       { return "timeStepCutFactor"; }
    static constexpr char const * maxAllowedResidualNormString()  { return "maxAllowedResidualNorm"; }

    static constexpr char const * numConfigurationAttemptsString()    { return "numConfigurationAttempts"; }
    static constexpr char const * maxNumConfigurationAttemptsString() { return "maxNumConfigurationAttempts"; }

    static constexpr char const * couplingTypeString()     { return "couplingType"; }
    static constexpr char const * subcyclingOptionString() { return "subcycling"; }
  } viewKeys;


  /**
   * @brief Calculates the upper limit for the number of iterations to allow a
   * decrease to the next time step.
   * @return The scaled value of the limit (m_timeStepDecreaseIterLimit * m_maxIterNewton)
   */
  integer timeStepDecreaseIterLimit() const
  {
    return std::ceil( m_timeStepDecreaseIterLimit * m_maxIterNewton );
  }


  /**
   * @brief Calculates the lower limit for the number of iterations to force an
   * increase to the next time step.
   * @return The scaled value of the limit (m_timeStepIncreaseIterLimit * m_maxIterNewton)
   */
  integer timeStepIncreaseIterLimit() const
  {
    return std::ceil( m_timeStepIncreaseIterLimit * m_maxIterNewton );
  }

  /**
   * @brief Getter for the factor used to decrease the time step size
   * @return the factor used to decrease the time step size
   */
  real64 timeStepDecreaseFactor() const
  {
    return m_timeStepDecreaseFactor;
  }

  /**
   * @brief Getter for the factor used to increase the time step size
   * @return the factor used to increase the time step size
   */
  real64 timeStepIncreaseFactor() const
  {
    return m_timeStepIncreaseFactor;
  }

  /**
   * @brief Getter for the norm type used to check convergence in the flow/well solvers
   * @return the norm type
   */
  solverBaseKernels::NormType normType() const
  {
    return m_normType;
  }

  /**
   * @brief Indicates the handling of line search in a Newton loop.
   */
  enum class LineSearchAction : integer
  {
    None,    ///< Do not use line search
    Attempt, ///< Use line search. Allow exit from line search without achieving smaller residual than starting residual.
    Require, ///< Use line search. If smaller residual than starting residual is not achieved, cut time step.
  };

  /**
   * @brief Indicates the handling of line each interpolation strategy.
   */
  enum class LineSearchInterpolationType : integer
  {
    Linear,    ///< linear decrease of line search scaling factor.
    Parabolic, ///< use parabolic interpolation to define line search scaling factor.
  };

  /**
   * @brief Coupling type.
   */
  enum class CouplingType : integer
  {
    FullyImplicit,      ///< Fully-implicit coupling
    Sequential ///< Sequential coupling
  };


  /// Flag to apply a line search.
  LineSearchAction m_lineSearchAction;

  /// Flag to pick the type of linesearch
  LineSearchInterpolationType m_lineSearchInterpType;

  /// The maximum number of line search cuts to attempt.
  integer m_lineSearchMaxCuts;

  /// The reduction factor for each line search cut.
  real64 m_lineSearchCutFactor;

  /// Norm used to check the nonlinear loop convergence
  solverBaseKernels::NormType m_normType;

  /// The tolerance for the nonlinear convergence check.
  real64 m_newtonTol;

  /// The maximum number of nonlinear iterations that are allowed.
  integer m_maxIterNewton;

  /// The minimum number of nonlinear iterations that may be applied.
  integer m_minIterNewton;

  /// The number of nonlinear iterations that have been exectued.
  integer m_numNewtonIterations;

  /// The maximum value of residual norm that we allow (otherwise, we cut the time step)
  real64 m_maxAllowedResidualNorm;

  /// Flag to allow for a non-converged nonlinear solution and continue with the problem.
  integer m_allowNonConverged;

  /// Fraction of the max Newton iterations above which the solver asks for the time-step to be decreased for the next timeStep
  real64 m_timeStepDecreaseIterLimit;

  /// Fraction of the max Newton iterations below which the solver asks for the time-step to be increased for the next timeStep
  real64 m_timeStepIncreaseIterLimit;

  /// Factor used to decrease the time step size
  real64 m_timeStepDecreaseFactor;

  /// Factor used to increase the time step size
  real64 m_timeStepIncreaseFactor;

  /// Maximum number of time sub-steps allowed for the solver
  integer m_maxSubSteps;

  /// Max number of time step cuts
  integer m_maxTimeStepCuts;

  /// Factor by which the time step will be cut if a timestep cut is required.
  real64 m_timeStepCutFactor;

  /// number of times that the time-step had to be cut
  integer m_numTimeStepAttempts;

  /// number of times that the configuration had to be changed
  integer m_numConfigurationAttempts;

  /// Max number of times that the configuration can be changed
  integer m_maxNumConfigurationAttempts;

  /// Type of coupling
  CouplingType m_couplingType;

  integer m_subcyclingOption;
};

ENUM_STRINGS( NonlinearSolverParameters::LineSearchAction,
              "None",
              "Attempt",
              "Require" );

ENUM_STRINGS( NonlinearSolverParameters::LineSearchInterpolationType,
              "Linear",
              "Parabolic" );

ENUM_STRINGS( NonlinearSolverParameters::CouplingType,
              "FullyImplicit",
              "Sequential" );



} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_ */
