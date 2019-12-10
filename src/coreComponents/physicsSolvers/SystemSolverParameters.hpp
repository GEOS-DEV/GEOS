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
 * @file SystemSolverParameters.hpp
 */

#ifndef GEOSX_SYSTEMSOLVERPARAMETERS_HPP_
#define GEOSX_SYSTEMSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{

class SystemSolverParameters : public dataRepository::Group
{
public:
  SystemSolverParameters() = delete;

  SystemSolverParameters( std::string const & name,
                          Group * const parent );

  SystemSolverParameters(SystemSolverParameters &&) = default;

  virtual ~SystemSolverParameters() override = default;

  static string CatalogName() { return "SystemSolverParameters"; }

  virtual void PostProcessInput() override;

  struct viewKeysStruct
  {
    static constexpr auto solverTypeString          = "solverType";
    static constexpr auto krylovTolString           = "krylovTol";
    static constexpr auto useAdaptiveKrylovString   = "useAdaptiveKrylovTol";
    static constexpr auto numKrylovIterString       = "numKrylovIter";
    static constexpr auto kspaceString              = "kspace";
    static constexpr auto ilut_fillString           = "ilut_fill";
    static constexpr auto ilut_dropString           = "ilut_drop";
    static constexpr auto useMLPrecondString        = "useMLPrecond";
    static constexpr auto useInnerSolverString      = "useInnerSolver";
    static constexpr auto scalingOptionString       = "scalingOption";
    static constexpr auto useBicgstabString         = "useBicgstab";
    static constexpr auto useDirectSolverString     = "useDirectSolver";
    static constexpr auto krylovResidualInitString  = "krylovResidualInit";
    static constexpr auto krylovResidualFinalString = "krylovResidualFinal";
    static constexpr auto useNewtonSolveString      = "useNewtonSolve";
    static constexpr auto newtonTolString           = "newtonTol";
    static constexpr auto maxIterNewtonString       = "maxIterNewton";
    static constexpr auto minIterNewtonString       = "minIterNewton";
    static constexpr auto numNewtonIterationsString = "numberOfNewtonIterations";
    static constexpr auto maxTimeStepCutsString     = "maxTimeStepCuts";
    static constexpr auto timeStepCutFactorString   = "timestepCutFactor";
    static constexpr auto maxLineSearchCutsString   = "maxLineSearchCuts";
    static constexpr auto lineSearchCutFactorString = "lineSearchCutFactor";
    static constexpr auto allowNonConvergedString   = "allowNonConverged";
    static constexpr auto maxSubStepsString         = "maxSubSteps";
    static constexpr auto dtCutIterLimString        = "dtCutIterLimit";
    static constexpr auto dtIncIterLimString        = "dtIncIterLimit";

  } viewKeys;

  struct groupKeysStruct
  {} groupKeys;

  string  solverType() const                  { return m_solverType; }
  real64 krylovTol() const                    { return m_krylovTol; }
  integer useAdaptiveKrylovTol() const        { return m_useAdaptiveKrylovTol; }
  integer  numKrylovIter() const              { return m_numKrylovIter; }
  integer  kspace() const                     { return m_kspace; }
  real64 ilut_fill() const                    { return m_ilut_fill; }
  real64 ilut_drop() const                    { return m_ilut_drop; }
  integer   useMLPrecond() const              { return m_useMLPrecond; }
  integer   useInnerSolver() const            { return m_useInnerSolver; }
  integer  scalingOption() const              { return m_scalingOption; }
  integer   useBicgstab() const               { return m_useBicgstab; }
  integer   useDirectSolver() const           { return m_useDirectSolver; }
  real64 krylovResidualInit() const           { return m_krylovResidualInit; }
  real64 krylovResidualFinal() const          { return m_krylovResidualFinal; }
  integer   useNewtonSolve() const            { return m_useNewtonSolve; }
  real64 newtonTol() const                    { return m_newtonTol; }
  integer  maxIterNewton() const              { return m_maxIterNewton; }
  integer  minIterNewton() const              { return m_minIterNewton; }
  integer & minIterNewton()                   { return m_minIterNewton; }
  integer const & numNewtonIterations() const { return m_numNewtonIterations; }
  integer & numNewtonIterations()             { return m_numNewtonIterations; }

  integer maxTimeStepCuts() const             { return m_maxTimeStepCuts; }
  real64  timeStepCutFactor() const           { return m_timeStepCutFactor; }
  integer maxLineSearchCuts() const           { return m_maxLineSearchCuts; }
  real64  lineSearchCutFactor() const         { return m_lineSearchCutFactor; }
  integer allowNonConverged() const           { return m_allowNonConverged; }
  integer maxSubSteps() const                 { return m_maxSubSteps; }
  integer & numdtAttempts()                   { return m_numdtAttempts; }
  real64  dtCutIterLimit() const              { return m_dtCutIterLimit * m_maxIterNewton; }
  real64  dtIncIterLimit() const              { return m_dtIncIterLimit * m_maxIterNewton; }

  string  m_solverType;
  real64  m_krylovTol;
  integer m_useAdaptiveKrylovTol;
  integer m_numKrylovIter;
  integer m_kspace;
  real64  m_ilut_fill;
  real64  m_ilut_drop;
  integer m_useMLPrecond;
  integer m_useInnerSolver;
  integer m_scalingOption;
  integer m_useBicgstab;
  integer m_useDirectSolver;
  real64  m_krylovResidualInit;
  real64  m_krylovResidualFinal;
  real64  m_krylovAuxTime;
  real64  m_krylovSetupTime;
  real64  m_krylovSolveTime;
  integer m_useNewtonSolve;
  real64  m_newtonTol;
  integer m_maxIterNewton;
  integer m_minIterNewton;
  integer m_numNewtonIterations;
  
  integer m_maxTimeStepCuts;
  real64  m_timeStepCutFactor;
  integer m_maxLineSearchCuts;
  real64  m_lineSearchCutFactor;
  integer m_allowNonConverged;
  integer m_maxSubSteps;
  integer m_maxIters = 1000;
  integer m_numdtAttempts; // number of times that the time-step had to be cut.
  real64  m_dtCutIterLimit;
  real64  m_dtIncIterLimit;
};

} /* namespace geosx */

#endif /*GEOSX_SYSTEMSOLVERPARAMETERS_HPP_*/
