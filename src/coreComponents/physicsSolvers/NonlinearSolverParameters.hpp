/*
 * NonlinearSolverParameters.hpp
 *
 *  Created on: Dec 14, 2019
 *      Author: settgast
 */

#ifndef GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_
#define GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{

class NonlinearSolverParameters : public dataRepository::Group
{
public:

  NonlinearSolverParameters() = delete;

  NonlinearSolverParameters( std::string const & name,
                             Group * const parent );
  virtual ~NonlinearSolverParameters() = default;

  NonlinearSolverParameters(NonlinearSolverParameters &&) = default;

  static string CatalogName() { return "NonlinearSolverParameters"; }

  virtual void PostProcessInput() override;

  struct viewKeysStruct
  {
    static constexpr auto lineSearchActionString        = "lineSearchAction";
    static constexpr auto lineSearchMaxCutsString       = "lineSearchMaxCuts";
    static constexpr auto lineSearchCutFactorString     = "lineSearchCutFactor";

    static constexpr auto newtonTolString               = "newtonTol";
    static constexpr auto newtonMaxIterString           = "newtonMaxIter";
    static constexpr auto newtonMinIterString           = "newtonMinIter";
    static constexpr auto newtonNumIterationsString     = "newtonNumberOfIterations";
    static constexpr auto newtonSplitOperMaxIterString  = "newtonSplitOperMaxIter";

    static constexpr auto allowNonConvergedString       = "allowNonConverged";
    static constexpr auto dtCutIterLimString            = "dtCutIterLimit";
    static constexpr auto dtIncIterLimString            = "dtIncIterLimit";
    static constexpr auto maxSubStepsString             = "maxSubSteps";
    static constexpr auto maxTimeStepCutsString         = "maxTimeStepCuts";
    static constexpr auto minNumNewtonIterationsString  = "minNumberOfNewtonIterations";
    static constexpr auto timeStepCutFactorString       = "timestepCutFactor";

  } viewKeys;


    real64  dtCutIterLimit() const              { return m_dtCutIterLimit * m_maxIterNewton; }
    real64  dtIncIterLimit() const              { return m_dtIncIterLimit * m_maxIterNewton; }

  integer m_lineSearchAction;
  integer m_lineSearchMaxCuts;
  real64 m_lineSearchCutFactor;

  real64 m_newtonTol;
  integer m_maxIterNewton;
  integer m_minIterNewton;
  integer m_numNewtonIterations;

  integer m_allowNonConverged;
  real64  m_dtCutIterLimit;
  real64  m_dtIncIterLimit;
  integer m_maxSubSteps;
  integer m_maxTimeStepCuts;
  real64  m_timeStepCutFactor;
  integer m_numdtAttempts; // number of times that the time-step had to be cut.


};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_ */
