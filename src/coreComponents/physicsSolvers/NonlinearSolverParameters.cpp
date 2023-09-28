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

#include "NonlinearSolverParameters.hpp"

namespace geos
{

using namespace dataRepository;

NonlinearSolverParameters::NonlinearSolverParameters( string const & name,
                                                      Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper( viewKeysStruct::lineSearchActionString(), &m_lineSearchAction ).
    setApplyDefaultValue( LineSearchAction::Attempt ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "How the line search is to be used. Options are: \n "
                    "* None    - Do not use line search.\n"
                    "* Attempt - Use line search. Allow exit from line search without achieving smaller residual than starting residual.\n"
                    "* Require - Use line search. If smaller residual than starting resdual is not achieved, cut time step." );

  registerWrapper( viewKeysStruct::lineSearchInterpolationTypeString(), &m_lineSearchInterpType ).
    setApplyDefaultValue( LineSearchInterpolationType::Linear ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Strategy to cut the solution update during the line search. Options are: \n "
                    "* Linear\n"
                    "* Parabolic" );

  registerWrapper( viewKeysStruct::lineSearchMaxCutsString(), &m_lineSearchMaxCuts ).
    setApplyDefaultValue( 4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of line search cuts." );

  registerWrapper( viewKeysStruct::lineSearchCutFactorString(), &m_lineSearchCutFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Line search cut factor. For instance, a value of 0.5 will result in the effective application of"
                    " the last solution by a factor of (0.5, 0.25, 0.125, ...)" );

  registerWrapper( viewKeysStruct::normTypeString(), &m_normType ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( solverBaseKernels::NormType::Linf ).
    setDescription( "Norm used by the flow solver to check nonlinear convergence. "
                    "Valid options:\n* " + EnumStrings< solverBaseKernels::NormType >::concat( "\n* " ) );

  registerWrapper( viewKeysStruct::newtonTolString(), &m_newtonTol ).
    setApplyDefaultValue( 1.0e-6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The required tolerance in order to exit the Newton iteration loop." );

  registerWrapper( viewKeysStruct::newtonMaxIterString(), &m_maxIterNewton ).
    setApplyDefaultValue( 5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of iterations that are allowed in a Newton loop." );

  registerWrapper( viewKeysStruct::newtonMinIterString(), &m_minIterNewton ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Minimum number of iterations that are required before exiting the Newton loop." );

  registerWrapper( viewKeysStruct::newtonNumIterationsString(), &m_numNewtonIterations ).
    setApplyDefaultValue( 0 ).
    setDescription( "Number of Newton's iterations." );

  registerWrapper( viewKeysStruct::maxAllowedResidualNormString(), &m_maxAllowedResidualNorm ).
    setApplyDefaultValue( 1e9 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum value of residual norm that is allowed in a Newton loop" );


  registerWrapper( viewKeysStruct::allowNonConvergedString(), &m_allowNonConverged ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Allow non-converged solution to be accepted. "
                    "(i.e. exit from the Newton loop without achieving the desired tolerance)" );

  registerWrapper( viewKeysStruct::timeStepDecreaseIterLimString(), &m_timeStepDecreaseIterLimit ).
    setApplyDefaultValue( 0.7 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fraction of the max Newton iterations above which the solver asks for the time-step to be decreased for the next time step." );

  registerWrapper( viewKeysStruct::timeStepIncreaseIterLimString(), &m_timeStepIncreaseIterLimit ).
    setApplyDefaultValue( 0.4 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fraction of the max Newton iterations below which the solver asks for the time-step to be increased for the next time step." );

  registerWrapper( viewKeysStruct::timeStepDecreaseFactorString(), &m_timeStepDecreaseFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which the time step is decreased when the number of Newton iterations is large." );

  registerWrapper( viewKeysStruct::timeStepIncreaseFactorString(), &m_timeStepIncreaseFactor ).
    setApplyDefaultValue( 2.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which the time step is increased when the number of Newton iterations is small." );

  registerWrapper( viewKeysStruct::timeStepCutFactorString(), &m_timeStepCutFactor ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Factor by which the time step will be cut if a timestep cut is required." );

  registerWrapper( viewKeysStruct::maxTimeStepCutsString(), &m_maxTimeStepCuts ).
    setApplyDefaultValue( 2 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max number of time step cuts" );

  registerWrapper( viewKeysStruct::maxSubStepsString(), &m_maxSubSteps ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of time sub-steps allowed for the solver" );

  registerWrapper( viewKeysStruct::maxNumConfigurationAttemptsString(), &m_maxNumConfigurationAttempts ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Max number of times that the configuration can be changed" );

  /// GEOS mainly uses FIM coupling so let's define FIM as the default.
  registerWrapper( viewKeysStruct::couplingTypeString(), &m_couplingType ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( CouplingType::FullyImplicit ).
    setDescription( "Type of coupling. "
                    "Valid options:\n* " + EnumStrings< CouplingType >::concat( "\n* " ) );

  registerWrapper( viewKeysStruct::sequentialConvergenceCriterionString(), &m_sequentialConvergenceCriterion ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( SequentialConvergenceCriterion::ResidualNorm ).
    setDescription( "Criterion used to check outer-loop convergence in sequential schemes. "
                    "Valid options:\n* " + EnumStrings< SequentialConvergenceCriterion >::concat( "\n* " ) );

  registerWrapper( viewKeysStruct::subcyclingOptionString(), &m_subcyclingOption ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to decide whether to iterate between sequentially coupled solvers or not." );
}

void NonlinearSolverParameters::postProcessInput()
{
  GEOS_ERROR_IF_LE_MSG( m_timeStepDecreaseIterLimit, m_timeStepIncreaseIterLimit,
                        getWrapperDataContext( viewKeysStruct::timeStepIncreaseIterLimString() ) <<
                        ": should be smaller than " << viewKeysStruct::timeStepDecreaseIterLimString() );
}



REGISTER_CATALOG_ENTRY( Group, NonlinearSolverParameters, string const &, Group * const )

} /* namespace geos */
