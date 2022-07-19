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

/*
 * @file SolverStatistics.cpp
 */

#include "SolverStatistics.hpp"

namespace geosx
{

using namespace dataRepository;

SolverStatistics::SolverStatistics( string const & name, Group * const parent )
  : Group( name, parent ),
  m_currentNumNonlinearIterations( 0 ),
  m_currentNumLinearIterations( 0 )
{
  registerWrapper( viewKeyStruct::numTimeStepsString(), &m_numTimeSteps ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Number of time steps" );

  registerWrapper( viewKeyStruct::numTimeStepCutsString(), &m_numTimeStepCuts ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Number of time step cuts" );

  registerWrapper( viewKeyStruct::numSuccessfulNonlinearIterationsString(), &m_numSuccessfulNonlinearIterations ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Cumulative number of successful nonlinear iterations" );

  registerWrapper( viewKeyStruct::numSuccessfulLinearIterationsString(), &m_numSuccessfulLinearIterations ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Cumulative number of successful linear iterations" );

  registerWrapper( viewKeyStruct::numFailedNonlinearIterationsString(), &m_numFailedNonlinearIterations ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Cumulative number of failed nonlinear iterations" );

  registerWrapper( viewKeyStruct::numFailedLinearIterationsString(), &m_numFailedLinearIterations ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Cumulative number of failed linear iterations" );
}

void SolverStatistics::initializeTimeStepStatistics()
{
  // the time step begins, we reset the individual-timestep counters
  m_currentNumNonlinearIterations = 0;
  m_currentNumLinearIterations = 0;
}

void SolverStatistics::logNonlinearIteration( integer const numLinearIterations )
{
  // we have just performed a Newton iteration, so we increment the individual-timestep counters
  m_currentNumNonlinearIterations++;
  m_currentNumLinearIterations += numLinearIterations;
}

void SolverStatistics::logNonlinearIteration()
{
  // we have just performed an outer iteration, so we increment the individual-timestep counter (number of outer iteration)
  m_currentNumNonlinearIterations++;
}

void SolverStatistics::logTimeStepCut()
{
  // we have just cut the time step, so we increment the cumulative counters for failed timesteps
  m_numFailedNonlinearIterations += m_currentNumNonlinearIterations;
  m_numFailedLinearIterations += m_currentNumLinearIterations;
  m_numTimeStepCuts++;

  // we are going to restart the timestep from the previous converged time step, so we have to re-initialize the statistics
  initializeTimeStepStatistics();
}

void SolverStatistics::saveTimeStepStatistics()
{
  // the timestep has converged, so we increment the cumulative counters for successful timesteps
  m_numSuccessfulNonlinearIterations += m_currentNumNonlinearIterations;
  m_numSuccessfulLinearIterations += m_currentNumLinearIterations;
  m_numTimeSteps++;
}

void SolverStatistics::outputStatistics() const
{
  bool const isExplicitScheme = m_numSuccessfulNonlinearIterations == 0 && m_numFailedNonlinearIterations == 0;
  bool const isOuterLoopSolver = m_numSuccessfulLinearIterations == 0 && m_numFailedLinearIterations == 0;

  GEOSX_LOG_RANK_0( getParent().getName() << ", number of time steps: " << m_numTimeSteps );
  if( !isExplicitScheme )
  {
    GEOSX_LOG_RANK_0( getParent().getName() << ", number of successful nonlinear iterations: " << m_numSuccessfulNonlinearIterations );
    if( !isOuterLoopSolver ) // don't print for the outer iterations in sequential schemes
    {
      GEOSX_LOG_RANK_0( getParent().getName() << ", number of successful linear iterations: " << m_numSuccessfulLinearIterations ); //
    }

    GEOSX_LOG_RANK_0( getParent().getName() << ", number of time step cuts: " << m_numTimeStepCuts );
    GEOSX_LOG_RANK_0( getParent().getName() << ", number of failed nonlinear iterations: " << m_numFailedNonlinearIterations );
    if( !isOuterLoopSolver ) // don't print for the outer iterations in sequential schemes
    {
      GEOSX_LOG_RANK_0( getParent().getName() << ", number of failed linear iterations: " << m_numFailedLinearIterations );
    }
  }
}
} // namespace geosx
