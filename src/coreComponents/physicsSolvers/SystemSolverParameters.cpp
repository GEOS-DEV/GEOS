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
 * @file SystemSolverParameters.cpp
 *
 */

#include "SystemSolverParameters.hpp"

namespace geosx
{
using namespace dataRepository;
SystemSolverParameters::SystemSolverParameters( std::string const & name,
                                                Group * const parent ):
  Group(name,parent),
  m_solverType("Klu"),
  m_krylovTol(),
  m_numKrylovIter(),
  m_kspace(),
  m_ilut_fill(3.0),
  m_ilut_drop(),
  m_useMLPrecond(),
  m_useInnerSolver(),
  m_scalingOption(),
  m_useBicgstab(),
  m_useDirectSolver(),
  m_krylovResidualInit(),
  m_krylovResidualFinal(),
  m_useNewtonSolve(),
  m_newtonTol(),
  m_maxIterNewton(),
  m_numNewtonIterations()
{
  setInputFlags(InputFlags::OPTIONAL);
  
  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper(viewKeysStruct::solverTypeString, &m_solverType, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::krylovTolString, &m_krylovTol, false )->
    setApplyDefaultValue(1.0e-6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Desired tolerance for Krylov solve");

  registerWrapper(viewKeysStruct::useAdaptiveKrylovString, &m_useAdaptiveKrylovTol, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Enable Eisenstat-Walker adaptive Krylov tolerance");

  registerWrapper(viewKeysStruct::numKrylovIterString, &m_numKrylovIter, false )->
    setApplyDefaultValue(100)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum number of Krylov Iterations");

  registerWrapper(viewKeysStruct::kspaceString, &m_kspace, false )->
    setApplyDefaultValue(300)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::ilut_fillString, &m_ilut_fill, false )->
    setApplyDefaultValue(3.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::ilut_dropString, &m_ilut_drop, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::useMLPrecondString, &m_useMLPrecond, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::useInnerSolverString, &m_useInnerSolver, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::scalingOptionString, &m_scalingOption, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::useBicgstabString, &m_useBicgstab, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::useDirectSolverString, &m_useDirectSolver, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::useNewtonSolveString, &m_useNewtonSolve, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::newtonTolString, &m_newtonTol, false )->
    setApplyDefaultValue(1.0e-6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::maxIterNewtonString, &m_maxIterNewton, false )->
    setApplyDefaultValue(5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum number of Newton iterations");

  registerWrapper(viewKeysStruct::minIterNewtonString, &m_minIterNewton, false )->
      setApplyDefaultValue(1)->
      setInputFlag(InputFlags::OPTIONAL)->
      setDescription("Minimum number of Newton iterations.");

  registerWrapper( viewKeysStruct::maxTimeStepCutsString, &m_maxTimeStepCuts, false )->
    setApplyDefaultValue(2)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Max number of time step cuts");

  registerWrapper( viewKeysStruct::timeStepCutFactorString, &m_timeStepCutFactor, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Time step cut factor");

  registerWrapper( viewKeysStruct::maxLineSearchCutsString, &m_maxLineSearchCuts, false )->
    setApplyDefaultValue(4)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Max number of line search cuts");

  registerWrapper( viewKeysStruct::lineSearchCutFactorString, &m_lineSearchCutFactor, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Line search cut factor");

  registerWrapper( viewKeysStruct::allowNonConvergedString, &m_allowNonConverged, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Allow non-converged solution to be accepted");

  registerWrapper( viewKeysStruct::maxSubStepsString, &m_maxSubSteps, false )->
    setApplyDefaultValue(10)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum number of time sub-steps allowed for the solver");

  registerWrapper(viewKeysStruct::krylovResidualInitString, &m_krylovResidualInit, false )->
    setApplyDefaultValue(0)->
    setDescription("Initial Krylov solver residual.");

  registerWrapper(viewKeysStruct::krylovResidualFinalString, &m_krylovResidualFinal, false )->
    setApplyDefaultValue(0)->
    setDescription("Final Krylov solver residual.");

  registerWrapper(viewKeysStruct::numNewtonIterationsString, &m_numNewtonIterations, false )->
    setApplyDefaultValue(0)->
    setDescription("Number of Newton's iterations.");

  registerWrapper(viewKeysStruct::dtCutIterLimString, &m_dtCutIterLimit, false )->
    setApplyDefaultValue(0.7)->
  setInputFlag(InputFlags::OPTIONAL)->
  setDescription("Fraction of the Max Newton iterations above which the solver asks for the time-step to be cut for the next dt.");

  registerWrapper(viewKeysStruct::dtIncIterLimString, &m_dtIncIterLimit, false )->
      setApplyDefaultValue(0.4)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Fraction of the Max Newton iterations below which the solver asks for the time-step to be doubled for the next dt.");

}

void SystemSolverParameters::PostProcessInput()
{
  if (m_dtCutIterLimit <= m_dtIncIterLimit)
  {
    GEOSX_ERROR(" dtIncIterLimit should be smaller than dtCutIterLimit!!" );
  }
}

REGISTER_CATALOG_ENTRY( Group, SystemSolverParameters, std::string const &, Group * const )

} /* namespace geosx */


