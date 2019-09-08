/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * SystemSolverParameters.cpp
 *
 *  Created on: Sep 12, 2017
 *      Author: settgast
 */

#include "SystemSolverParameters.hpp"

namespace geosx
{
using namespace dataRepository;
SystemSolverParameters::SystemSolverParameters( std::string const & name,
                                                Group * const parent ):
  Group(name,parent),
  m_verbose(0),
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
  m_KrylovResidualInit(),
  m_KrylovResidualFinal(),
  m_useNewtonSolve(),
  m_newtonTol(),
  m_maxIterNewton(),
  m_numNewtonIterations()
{
  setInputFlags(InputFlags::OPTIONAL);
  
  setRestartFlags(RestartFlags::NO_WRITE);
  registerWrapper(viewKeysStruct::verbosityString, &m_verbose, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("verbosity level");

  registerWrapper(viewKeysStruct::solverTypeString, &m_solverType, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("");

  registerWrapper(viewKeysStruct::krylovTolString, &m_krylovTol, false )->
    setApplyDefaultValue(1.0e-6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Allowable tolerance for krylov solve");

  registerWrapper(viewKeysStruct::numKrylovIterString, &m_numKrylovIter, false )->
    setApplyDefaultValue(100)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum number of Krylov Iterations");

  registerWrapper(viewKeysStruct::kspaceString, &m_kspace, false )->
    setApplyDefaultValue(0)->
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

  registerWrapper(viewKeysStruct::KrylovResidualInitString, &m_KrylovResidualInit, false )->
    setApplyDefaultValue(0)->
    setDescription("verbosity level");

  registerWrapper(viewKeysStruct::KrylovResidualFinalString, &m_KrylovResidualFinal, false )->
    setApplyDefaultValue(0)->
    setDescription("verbosity level");

  registerWrapper(viewKeysStruct::numNewtonIterationsString, &m_numNewtonIterations, false )->
    setApplyDefaultValue(0)->
    setDescription("verbosity level");


}


REGISTER_CATALOG_ENTRY( Group, SystemSolverParameters, std::string const &, Group * const )

} /* namespace geosx */


